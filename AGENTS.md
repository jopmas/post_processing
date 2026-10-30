# post_processing

Scripts and notebooks to post-process Mandyoc models from the NetCDF files created by tapioca
(`src/postrunning/meshes2netcdf.jl` and `src/postrunning/particles2netcdf.jl`, https://github.com/ggciag/tapioca).

## Layout

- `track_particles_nc_dk.py`, `track_particles_nc_rft.py`: select particles in the final step (double keel and
  rifting models), read their trajectories from `particles_trajectories.nc`, interpolate P and T from
  `pressure.nc`/`temperature.nc` and write `_track_xzPT_all_steps.nc`.
- `plot_PTT_frames_dk.py`, `plot_PTT_frames_rft.py`: frames of the lithology with the P-T-t paths of
  `_track_xzPT_all_steps.nc` (temperature-coded for double keel, depth-coded for rifting), videos and gifs.
- `frames_generator.py`: frames of the lithology, strain, isotherms and melt (accordion and rifting models).
- `Geometries.ipynb`: a field at chosen instants (plotted with `plot_property`), crustal thickness and width of the
  rifted margins. It has no functions of its own: they are in `functions/mandyocIO.py`.
- `functions/mandyocIO.py`: functions shared by the scripts and the notebook: lithology colors
  (`LITHOLOGY_COLORS`), readers of the tapioca files (`read_tapioca_datasets`, `read_params`, `read_times_steps`,
  `read_lithology`), crust measurements (`measure_crust_from_lithology`, `measure_crustal_thickness`,
  `measure_margins_width`), `plot_lithology_legend`, `add_colorbar`, `plot_property`, the P-T-t path functions and
  helpers (`find_nearest`, `change_dataset`, `trim_axs`).
- `docs/`: user guide and references (Markdown).

## Data conventions (tapioca files)

- Mesh fields (`<var>.nc`): data with dims `(time, z, x)`, `x` from 0 to Lx and `z` from -Lz to 0 [m], so the first
  row is the bottom. The step is the coordinate `step` and the time is in Myr. There are no `nx`/`lx` attributes:
  get `Nx`, `Nz`, `Lx`, `Lz` from `sizes` and the coordinates.
- File names are the variable names: `Phi.nc` (melt), `dPhi.nc` (incremental melt), `X_depletion.nc`, `heat.nc`,
  `velocity.nc` (`vx`, `vz`). Only `surface.nc` has its own `x` (more points than the mesh): interpolate it to the
  mesh x.
- tapioca keeps the raw values in the air. Where a figure needs the air empty (strain shading, melt), set the
  values to zero where `density < 1200` kg/m3, as the old Julia files did.
- `lithology.nc` has the Mandyoc layer codes (order of the layers in `interfaces.txt`, 0 = bottom layer) on a mesh
  5 times finer: `Nxl = 5*(Nx-1)+1`, `Nzl = 5*(Nz-1)+1`, and `lithology[:, ::5]` is on the mesh x. The lithology is
  plotted from `lithology.nc`, never from the density. `plot_property` and the notebook use `LITHOLOGY_COLORS` of
  `functions/mandyocIO.py`; the frame scripts build a `lithology_colors` dict from their layer code variables. Plot with
  `ListedColormap(colors)`, `vmin=-0.5`, `vmax=len(colors)-0.5`. The codes depend on the model type
  (double keel: 11 layers; rifting: 7; accordion AR: 8, with a basaltic layer of 2900 kg/m3).
- `particles_trajectories.nc`: `x`, `z`, `lithology` with dims `(time, id)`, `-999` (and `-1` for the lithology)
  where the particle does not exist in a step. It only has the particles that exist in its reference step
  (attribute `reference_timestep`): the tracking scripts select particles in the final step, so the file should be
  created with the final step as reference.
- `_track_xzPT_all_steps.nc` is read by the frame scripts and by `functions/mandyocIO.py`: keep its structure
  (`xtrack`, `ztrack`, `ptrack` [MPa], `ttrack` [oC] flattened as (time, particle) and reversed, `step`, `time`,
  `ntracked`, `particles_layers`).

## Conventions

- Keep the structure and the line of reasoning of the scripts: settings at the top (flags, region, layer codes,
  colors), the same variable names, `pymp` for the frames and the tracking. Prefer small changes over rewrites.
- The scripts run in the folder of the scenario (`model_path = os.getcwd()`) and import `functions/mandyocIO.py`
  from the folder of the script (`path_to_functions`).
- Legends and colorbars use `plot_lithology_legend` and `add_colorbar`: their size is fixed in inches (legend
  0.74 x 0.66 in, colorbar 1.8 x 0.11 in), only the position is relative to the axis.
- With `pymp` (fork), read the HDF5/NetCDF data needed by the parallel part before it, or open the files inside it.
- The tracking scripts zip `track_xzPT_step_*.txt` with `zip -m` (the files are removed after they are in the zip).
- `track_particles_nc_rft.py` sets `take_asthenosphere` and `instant_to_take` from the model name (`He`, `W1`/`S1`,
  ...). `frames_generator.py` has `oceanic_crust` (accordion AR models) and `plot_PTT_frames_dk.py` has
  `double_keel` to choose the layer codes.
- `Geometries.ipynb` measures the crust counting the cells of `crust_codes` in each column of `lithology.nc`
  (thickness = number of cells x dz/5), as `thickness_oceanic_crust` in `frames_generator_litho_melt_thickness.py`.
- Code, comments, docstrings, commits and docs are written in English.

## Testing

There is no test suite. Test a change by running the scripts on a scenario converted by tapioca with a few steps
(link the `.txt` folders of a real scenario into a temporary folder with only some `time/time_*.txt` files and run
`meshes2netcdf.jl` there), never writing in the folders of the user's scenarios. Compare the outputs with the
previous version (the `_track_xzPT_all_steps.nc` and the frames) and look at the figures.

## Commands

```bash
conda activate tapioca                                                   # environment of tapioca
julia -t 14 ~/opt/tapioca/src/postrunning/meshes2netcdf.jl ${PWD}       # in the folder of the scenario
julia -t 14 ~/opt/tapioca/src/postrunning/particles2netcdf.jl ${PWD} <reference_step>
python ~/opt/post_processing/track_particles_nc_dk.py                   # in the folder of the scenario
python ~/opt/post_processing/plot_PTT_frames_dk.py
```

Mandyoc runs are done on the group's cluster, Hypatia.
