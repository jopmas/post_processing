"""
Zero-offset acoustic propagation with Devito, as an alternative to tremelique.

``DevitoZeroOffset`` solves the same variable-density acoustic wave equation
as ``tremelique.Acoustic``,

    1 / (rho v^2) d2u/dt2 = div(1 / rho grad u),

with the finite-difference compiler Devito (Louboutin et al., 2019,
https://www.devitoproject.org). Devito generates and compiles C code for the
stencil, with OpenMP threads, so any space order can be used; here 8th order in
space (staggered first derivatives for the density term), 2nd order in time.

It has the same interface as the ``ZeroOffset`` class of
``seismic_section.py`` (exploding-reflector source ``source_scale`` times the
wavelet, injection of traces at the receivers for RTM, traces of the receivers
and last wave field kept in memory), so the two propagators are
interchangeable. The differences are the discretization (8th instead of 4th
order in space) and the absorbing boundaries: a damping term ``sigma du/dt``
with a quadratic profile in the padding zone (left, right and bottom), whose
strength is set for a reflection of ``reflection`` at normal incidence
(Cerjan et al., 1985; Sochacki et al., 1987), instead of tremelique's
multiplicative taper.

On macOS, Apple's clang needs to be told where ``omp.h`` and ``libomp`` are:
if ``CFLAGS``/``LDFLAGS`` are not set, they are set to the ``include`` and
``lib`` folders of the active Python environment (conda's ``llvm-openmp``).
"""

import os
import sys
from pathlib import Path

if sys.platform == "darwin" and "CFLAGS" not in os.environ:
    _prefix = Path(sys.prefix)
    if (_prefix / "include" / "omp.h").exists():
        os.environ["CFLAGS"] = f"-Xclang -fopenmp -I{_prefix / 'include'}"
        os.environ.setdefault(
            "LDFLAGS", f"-L{_prefix / 'lib'} -lomp -Wl,-rpath,{_prefix / 'lib'}"
        )
        os.environ.setdefault("DEVITO_LANGUAGE", "openmp")
os.environ.setdefault("DEVITO_LOGGING", "WARNING")

import numpy as np  # noqa: E402
import sympy  # noqa: E402
from devito import (  # noqa: E402
    Eq,
    Function,
    Grid,
    Operator,
    SparseTimeFunction,
    TimeFunction,
    div,
    grad,
    solve,
)


def ricker(f_cut):
    """Causal Ricker wavelet with maximum frequency ``f_cut`` (as tremelique's)."""
    sqrt_pi = np.sqrt(np.pi)
    fc = f_cut / (3 * sqrt_pi)

    def wavelet(time):
        t = np.pi * fc * (np.asarray(time) - 2 * sqrt_pi / f_cut)
        return (2 * np.pi * t**2 - 1) * np.exp(-np.pi * t**2)

    return wavelet


def critical_dt(spacing, vmax, space_order, safety=0.8):
    """
    Largest stable time step of the leapfrog scheme with staggered derivatives.

    The staggered first derivative with weights ``c_k`` (on both sides) has a
    largest eigenvalue ``sum|c_k| / h`` (at the Nyquist wavenumber), so the 2D
    operator has ``2 (sum|c_k| / h)^2`` and the leapfrog scheme is stable for
    ``dt <= sqrt(2) h / (vmax sum|c_k|)`` (0.55 h / vmax for 8th order).
    """
    half = space_order // 2
    points = [k - 0.5 for k in range(-half + 1, half + 1)]
    weights = sympy.finite_diff_weights(1, points, 0)[1][-1]
    c = float(sum(abs(w) for w in weights))
    return safety * np.sqrt(2) * spacing / (vmax * c)


def _set(function, values):
    """Fill a Function and its halo (edge values), so stencils near the edges are smooth."""
    halo = [tuple(h) for h in function._size_halo]
    function.data_with_halo[:] = np.pad(values, halo, mode="edge")


class DevitoZeroOffset:
    """
    Acoustic zero-offset simulation (exploding reflector or RTM) with Devito.

    Parameters
    ----------
    velocity, density : 2D arrays (nz, nx)
        Model, including the padding of the absorbing zone (bottom, left and
        right, ``padding`` cells).
    spacing : float
        Grid spacing [m] (the same in x and z).
    receivers : tuple of arrays
        Rows and columns of the receivers.
    absorb : 2D array
        Factor applied to the wave field at every step (1 where nothing is
        absorbed; used to absorb the water above the receivers).
    f_cut : float
        Maximum frequency of the Ricker wavelet [Hz].
    padding : int
        Width of the absorbing zone [cells].
    dt : float, optional
        Time step [s]; by default the stable one.
    space_order : int
        Order of the spatial derivatives.
    reflection : float
        Theoretical reflection coefficient of the absorbing zone.
    """

    def __init__(self, velocity, density, spacing, receivers, absorb, f_cut, padding,
                 dt=None, space_order=8, reflection=1e-4):  # fmt: skip
        self.velocity = np.asarray(velocity, dtype=np.float32)
        self.density = np.asarray(density, dtype=np.float32)
        self.shape = self.velocity.shape
        self.spacing = float(spacing)
        self.rec_z, self.rec_x = (np.asarray(r) for r in receivers)
        self.absorb = np.asarray(absorb, dtype=np.float32)
        self.padding = padding
        self.space_order = space_order
        self.reflection = reflection
        self.wavelet = ricker(f_cut)
        self.dt = dt or critical_dt(self.spacing, self.velocity.max(), space_order)
        self.source_scale = None  # exploding reflector: 2D field x wavelet
        self.injection = None  # RTM: traces injected at the receivers
        self.traces = None
        self.last = None

    def damping(self):
        """Damping coefficient sigma [1/s]: quadratic in the padding zone."""
        nz, nx = self.shape
        pad = self.padding
        width = pad * self.spacing
        sigma_max = 3 * self.velocity.max() * np.log(1 / self.reflection) / (2 * width)
        dist = np.zeros(self.shape)
        cols = np.arange(nx)
        rows = np.arange(nz)
        dist = np.maximum(dist, np.clip(pad - cols, 0, None)[None, :])
        dist = np.maximum(dist, np.clip(cols - (nx - 1 - pad), 0, None)[None, :])
        dist = np.maximum(dist, np.clip(rows - (nz - 1 - pad), 0, None)[:, None])
        return (sigma_max * (dist / pad) ** 2).astype(np.float32)

    def run(self, nsteps):
        """Run ``nsteps`` time steps; keeps ``traces`` (nsteps, nrec) and ``last``."""
        nz, nx = self.shape
        h = self.spacing
        so = self.space_order
        grid = Grid(shape=(nz, nx), extent=((nz - 1) * h, (nx - 1) * h),
                    dtype=np.float32)  # fmt: skip
        time = grid.time_dim
        u = TimeFunction(name="u", grid=grid, time_order=2, space_order=so)
        fields = {}
        for name, values in [("m", 1 / self.velocity**2), ("rho", self.density),
                             ("b", 1 / self.density), ("sigma", self.damping()),
                             ("absorb", self.absorb)]:  # fmt: skip
            fields[name] = Function(name=name, grid=grid, space_order=so)
            _set(fields[name], values.astype(np.float32))
        m, rho, b, sigma, absorb = (fields[k] for k in ("m", "rho", "b", "sigma", "absorb"))

        pde = m * u.dt2 + m * sigma * u.dt - rho * div(b * grad(u, shift=0.5), shift=-0.5)
        # As in tremelique, the factor is applied to u(t+dt) and again to u(t)
        # after the step, so u(t-dt) enters the update multiplied twice: the
        # asymmetry damps du/dt, which is what absorbs the waves
        update = absorb * solve(pde, u.forward).subs(u.backward, absorb * u.backward)

        coords = np.stack([self.rec_z * h, self.rec_x * h], axis=1).astype(np.float32)
        npoint = len(coords)
        equations = []
        if self.source_scale is not None:
            src = Function(name="src", grid=grid, space_order=0)
            src.data[:] = self.source_scale
            wav = Function(name="wav", dimensions=(time,), shape=(nsteps + 1,),
                           dtype=np.float32)  # fmt: skip
            wav.data[:] = self.wavelet(np.arange(nsteps + 1) * self.dt)
            update = update + src * wav[time]
        equations.append(Eq(u.forward, update))
        if self.injection is not None:
            inj = SparseTimeFunction(name="inj", grid=grid, npoint=npoint, nt=nsteps + 1,
                                     coordinates=coords)  # fmt: skip
            n = min(len(self.injection), nsteps + 1)
            inj.data[:n] = self.injection[:n]
            equations += inj.inject(field=u.forward, expr=inj)
        rec = SparseTimeFunction(name="rec", grid=grid, npoint=npoint, nt=nsteps + 1,
                                 coordinates=coords)  # fmt: skip
        equations += rec.interpolate(expr=u.forward)

        op = Operator(equations, subs=grid.spacing_map, name="ZeroOffset")
        op.apply(time_m=0, time_M=nsteps - 1, dt=np.float32(self.dt))
        self.traces = np.array(rec.data[:nsteps])
        self.last = np.array(u.data[nsteps % 3])
        return self
