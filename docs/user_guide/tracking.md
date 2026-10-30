# Particle tracking

`track_particles_nc_dk.py` (double keel models) and `track_particles_nc_rft.py` (rifting models):

1. select particles in the final step of the tracking: particles above `h_air + search_thickness` (depth below the
   top of the model), between `x_begin` and `x_end`, that are not air (and not asthenosphere, unless
   `take_asthenosphere = True`);
2. read the positions of these particles in every step from `particles_trajectories.nc`;
3. interpolate (bilinear) the pressure and the temperature of `pressure.nc` and `temperature.nc` at the positions;
4. write `_track_xzPT_all_steps.nc` (see [Data format](../references/data_format.md)), one text file per step
   (`track_xzPT_step_<step>.txt`, zipped in `<model>.zip` and then removed) and, for double keel models,
   `infos_track_<model>.txt`.

Input files: `temperature.nc`, `pressure.nc` and `particles_trajectories.nc`, created with the final step of the
tracking as reference step (see [Workflow](workflow.md)).

## Settings

### `track_particles_nc_dk.py`

| Variable | Meaning |
| --- | --- |
| `take_specific_time`, `instant_to_take` | Track until `instant_to_take` [Myr] instead of the last step. |
| `take_asthenosphere` | Also select asthenosphere particles. |
| `h_air` | Thickness of the air layer [m]. |
| `search_thickness` | Select particles down to this depth below the surface [m]. |
| `x_begin`, `x_end` | Horizontal limits of the selection [m]. |
| `double_keel` | Layer codes of the double keel model (11 layers) or of the homogeneous model (9 layers). |

### `track_particles_nc_rft.py`

The settings depend on the name of the model: `take_asthenosphere = True` for models with `He` in the name, and
`instant_to_take` (40 Myr) is chosen from the extension rate (`W1`/`S1`, `W2`/`S2`, `W3`/`S3`). The selection uses
`search_thickness = 12` km, `x_begin = 600` km and `x_end = 1600` km, and the layer codes of the rifting models
(asthenosphere 0, lithospheric mantle 1 and 3, seed 2, lower crust 4, upper crust 5, air 6).

## Notes

- `particles_trajectories.nc` has `-999` where a particle does not exist: the positions, P and T are `NaN` there.
- The positions are read before the `pymp` block, since NetCDF/HDF5 files must not be read by processes forked
  after opening them.
