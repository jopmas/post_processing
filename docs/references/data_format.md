# Data format

## Mesh fields (`meshes2netcdf.jl`)

One file per field, named after the variable:

| File | Variable(s) | Units |
| --- | --- | --- |
| `density.nc` | `density` | kg/m3 |
| `temperature.nc` | `temperature` | oC |
| `pressure.nc` | `pressure` | Pa |
| `strain.nc` | `strain` | - |
| `strain_rate.nc` | `strain_rate` | 1/s |
| `viscosity.nc` | `viscosity` | Pa.s |
| `heat.nc` | `heat` | W/kg |
| `velocity.nc` | `vx`, `vz` | m/s |
| `Phi.nc`, `dPhi.nc`, `X_depletion.nc` | `Phi` (melt), `dPhi` (incremental melt), `X_depletion` | - |
| `surface.nc` | `surface` | m |
| `lithology.nc` | `lithology` | Mandyoc layer code |

- Dims `(time, z, x)`: `x` from 0 to Lx, `z` from -Lz to 0 [m] (positive upwards), so the first row is the bottom
  of the model. `time` in Myr and the coordinate `step`.
- There are no `nx`, `nz`, `lx`, `lz` attributes: use `sizes` and the coordinates.

  ```python
  ds = xr.open_dataset("strain.nc")
  Nx, Nz = ds.sizes["x"], ds.sizes["z"]
  Lx = float(ds.x[-1] - ds.x[0]); Lz = float(ds.z[-1] - ds.z[0])
  ```

- The values in the air are kept (the old Julia scripts set them to zero). Mask them with the density where
  needed: `np.where(density < 1200, 0, strain)`.
- `surface.nc` has dims `(time, x)` with more points than the mesh; it is the z of the surface in the coordinates
  of the model (-40 km at the sea level for 40 km of air).
- `lithology.nc` is 5 times finer: `Nxl = 5*(Nx-1)+1`, `Nzl = 5*(Nz-1)+1`. The codes are the order of the layers in
  `interfaces.txt` (0 = bottom layer), for example:

  | Model | Codes |
  | --- | --- |
  | double keel | 0 asthenosphere, 1 lower craton, 2 upper craton, 3 lithospheric mantle, 4 lower crust 1, 5 seed, 6 lower crust 2, 7 upper crust, 8 decolement, 9 sediments, 10 air |
  | rifting | 0 asthenosphere, 1 lithospheric mantle, 2 seed, 3 lithospheric mantle, 4 lower crust, 5 upper crust, 6 air |
  | accordion (AR) | as rifting, with 6 basaltic crust (2900 kg/m3) and 7 air |

## Particles (`particles2netcdf.jl`)

`particles_trajectories.nc`:

- `x`, `z` [m] and `lithology` (layer code) with dims `(time, id)`; `id` is the particle id (sorted).
- `-999` in `x`/`z` and `-1` in `lithology` where the particle does not exist in a step.
- Only the particles that exist in the reference step (attribute `reference_timestep`), inside `SELECT_X`,
  `SELECT_Z` and `SELECT_LITHOLOGY` of the Julia script.

## Tracked particles (`_track_xzPT_all_steps.nc`)

Written by the tracking scripts and read by the frame scripts and `functions/mandyocIO.py`:

| Variable | Meaning |
| --- | --- |
| `xtrack`, `ztrack` | Positions [m], dim `index`. |
| `ptrack`, `ttrack` | Pressure [MPa] and temperature [oC], dim `index`. |
| `step`, `time` | Steps and times [Myr] of the tracking. |
| `ntracked` | Number of tracked particles `n`. |
| `particles_layers` | Layer code of each tracked particle. |

The `index` dimension is `(time, particle)` flattened and reversed. To read it:

```python
track = xr.open_dataset("_track_xzPT_all_steps.nc")
n = int(track.ntracked)
T = track.ttrack.values[::-1].reshape(-1, n)  # (time, particle)
layers = track.particles_layers.values[::-1]
```
