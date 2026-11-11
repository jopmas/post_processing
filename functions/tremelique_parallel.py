"""
Multithreaded version of the acoustic simulation of tremelique.

The kernels below are the ones of ``tremelique._acoustic`` (``timestep_esg``,
``nonreflexive_bc``) and ``tremelique._utils`` (``apply_damping``), with the
loops over the grid split between threads with ``numba.prange``
(``parallel=True``). Each cell is updated with exactly the same operations, so
the results are identical to the serial ones; only the work is shared between
the cores. The number of threads is that of numba (all the cores by default,
``numba.set_num_threads`` or the ``NUMBA_NUM_THREADS`` environment variable
change it).

``ParallelAcoustic`` is a drop-in replacement for ``tremelique.Acoustic``
(same arguments, HDF5 cache, plots and animations).

The kernels are adapted from tremelique (https://www.fatiando.org/tremelique),
Copyright (c) 2025 The Tremelique Developers, distributed under the terms of
the BSD 3-Clause License.
"""

import numba
import numpy as np
import tremelique as tr
from numba import prange


@numba.jit(nopython=True, parallel=True)
def timestep_esg(u_tp1, u_t, u_tm1, x1, x2, z1, z2, dt, dx, dz, vel, dens):
    """
    One time step of the equivalent staggered grid method (DiBartolo et al., 2012).

    Same as ``tremelique._acoustic.timestep_esg``, with the rows of the grid
    split between threads.
    """
    dt2 = dt**2
    dx2 = dx**2
    dz2 = dz**2
    for i in prange(z1, z2):
        for j in range(x1, x2):
            zderiv = (1.125 / dz2) * (
                0.5
                * (1 / dens[i + 1, j] + 1 / dens[i, j])
                * (
                    1.125 * (u_t[i + 1, j] - u_t[i, j])
                    - (u_t[i + 2, j] - u_t[i - 1, j]) / 24.0
                )
                - 0.5
                * (1 / dens[i, j] + 1 / dens[i - 1, j])
                * (
                    1.125 * (u_t[i, j] - u_t[i - 1, j])
                    - (u_t[i + 1, j] - u_t[i - 2, j]) / 24.0
                )
            ) - (1 / (24 * dz2)) * (
                0.5
                * (1 / dens[i + 2, j] + 1 / dens[i + 1, j])
                * (
                    1.125 * (u_t[i + 2, j] - u_t[i + 1, j])
                    - (u_t[i + 3, j] - u_t[i, j]) / 24.0
                )
                - 0.5
                * (1 / dens[i - 1, j] + 1 / dens[i - 2, j])
                * (
                    1.125 * (u_t[i - 1, j] - u_t[i - 2, j])
                    - (u_t[i, j] - u_t[i - 3, j]) / 24.0
                )
            )
            xderiv = (1.125 / dx2) * (
                0.5
                * (1 / dens[i, j + 1] + 1 / dens[i, j])
                * (
                    1.125 * (u_t[i, j + 1] - u_t[i, j])
                    - (u_t[i, j + 2] - u_t[i, j - 1]) / 24.0
                )
                - 0.5
                * (1 / dens[i, j] + 1 / dens[i, j - 1])
                * (
                    1.125 * (u_t[i, j] - u_t[i, j - 1])
                    - (u_t[i, j + 1] - u_t[i, j - 2]) / 24.0
                )
            ) - (1.0 / (24.0 * dx2)) * (
                0.5
                * (1 / dens[i, j + 2] + 1 / dens[i, j + 1])
                * (
                    1.125 * (u_t[i, j + 2] - u_t[i, j + 1])
                    - (u_t[i, j + 3] - u_t[i, j]) / 24.0
                )
                - 0.5
                * (1 / dens[i, j - 1] + 1 / dens[i, j - 2])
                * (
                    1.125 * (u_t[i, j - 1] - u_t[i, j - 2])
                    - (u_t[i, j] - u_t[i, j - 3]) / 24.0
                )
            )
            u_tp1[i, j] = (
                2 * u_t[i, j]
                - u_tm1[i, j]
                + dt2 * dens[i, j] * vel[i, j] ** 2 * (xderiv + zderiv)
            )


@numba.jit(nopython=True, parallel=True)
def apply_damping(array, nx, nz, pad, decay):
    """
    Damping of the padding regions (left, right and bottom).

    Same as ``tremelique._utils.apply_damping``, split between threads.
    """
    for i in prange(nz):
        for j in range(pad):
            array[i, j] *= np.exp(-((decay * (pad - j)) ** 2))
    for i in prange(nz):
        for j in range(nx - pad, nx):
            array[i, j] *= np.exp(-((decay * (j - nx + pad)) ** 2))
    for i in prange(nz - pad, nz):
        for j in range(nx):
            array[i, j] *= np.exp(-((decay * (i - nz + pad)) ** 2))


@numba.jit(nopython=True, parallel=True)
def nonreflexive_bc(u_tp1, u_t, nx, nz, dt, dx, dz, mu, dens):
    """
    Non-reflexive boundary conditions.

    Same as ``tremelique._acoustic.nonreflexive_bc``, split between threads.
    The four borders are done in the same order as in tremelique, so the
    corners get the same values.
    """
    # Left
    for i in prange(nz):
        for j in range(3):
            u_tp1[i, j] = (
                u_t[i, j]
                + dt * np.sqrt(mu[i, j] / dens[i, j]) * (u_t[i, j + 1] - u_t[i, j]) / dx
            )
    # Right
    for i in prange(nz):
        for j in range(nx - 3, nx):
            u_tp1[i, j] = (
                u_t[i, j]
                - dt * np.sqrt(mu[i, j] / dens[i, j]) * (u_t[i, j] - u_t[i, j - 1]) / dx
            )
    # Bottom
    for j in prange(nx):
        for i in range(nz - 3, nz):
            u_tp1[i, j] = (
                u_t[i, j]
                - dt * np.sqrt(mu[i, j] / dens[i, j]) * (u_t[i, j] - u_t[i - 1, j]) / dz
            )
    # Top
    for j in prange(nx):
        u_tp1[2, j] = u_tp1[3, j]
        u_tp1[1, j] = u_tp1[2, j]
        u_tp1[0, j] = u_tp1[1, j]


@numba.jit(nopython=True, parallel=True)
def multiply_inplace(a, b):
    """``a *= b`` for 2D arrays, split between threads."""
    for i in prange(a.shape[0]):
        for j in range(a.shape[1]):
            a[i, j] *= b[i, j]


@numba.jit(nopython=True, parallel=True)
def add_scaled_inplace(a, b, scale):
    """
    ``a += scale * b`` for 2D arrays, split between threads.

    The sum is done in float64 and rounded to the type of ``a``, as NumPy does
    for ``a += b * scale`` with a float32 ``a`` and a float64 ``scale``, so the
    result is the same as with NumPy.
    """
    for i in prange(a.shape[0]):
        for j in range(a.shape[1]):
            a[i, j] = np.float64(a[i, j]) + np.float64(b[i, j]) * np.float64(scale)


class ParallelAcoustic(tr.Acoustic):
    """
    ``tremelique.Acoustic`` with the time step split between threads.

    Same arguments, methods and results as ``tremelique.Acoustic``; only
    ``_timestep`` changes, to use the multithreaded kernels of this module.
    """

    def _timestep(self, u, tm1, t, tp1, iteration):
        nz, nx = self.shape
        timestep_esg(
            u[tp1], u[t], u[tm1], 3, nx - 3, 3, nz - 3,
            self.dt, self.dx, self.dz, self.velocity, self.density,
        )  # fmt: skip
        apply_damping(u[t], nx, nz, self.padding, self.taper)
        nonreflexive_bc(
            u[tp1], u[t], nx, nz, self.dt, self.dx, self.dz, self.velocity, self.density
        )
        apply_damping(u[tp1], nx, nz, self.padding, self.taper)
        for pos, src in self.sources:
            i, j = pos
            scale = -self.density[i, j] * (self.velocity[i, j] * self.dt) ** 2
            u[tp1, i, j] += scale * src(iteration * self.dt)
