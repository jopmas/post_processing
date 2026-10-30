# post_processing

Scripts and notebooks to post-process [Mandyoc](https://github.com/ggciag/mandyoc) models from the NetCDF files
created by [tapIOca](https://github.com/ggciag/tapioca): particle tracking (P-T-t paths), frames and videos of the
lithology, and measurements of the crust (crustal thickness and width of the rifted margins).

## About

The scripts read **only** the `.nc` files created by the Julia scripts of tapIOca
(`src/postrunning/meshes2netcdf.jl` and `src/postrunning/particles2netcdf.jl`): `temperature.nc`, `strain.nc`,
`lithology.nc`, `particles_trajectories.nc`, ... The old files (`_strain.nc`, `_output_strain.nc`, `_steps.nc`) are
not read. See [docs/references/data_format.md](docs/references/data_format.md) for the structure of these files.

| Script | What it does |
| --- | --- |
| `track_particles_nc_dk.py` | Selects particles of double keel models and interpolates P and T along their trajectories (`_track_xzPT_all_steps.nc`). |
| `track_particles_nc_rft.py` | The same for rifting models. |
| `plot_PTT_frames_dk.py` | Frames of the lithology with the temperature-coded P-T-t paths (double keel models). |
| `plot_PTT_frames_rft.py` | Frames of the lithology with the depth-coded P-T-t paths (rifting models). |
| `frames_generator.py` | Frames of the lithology, strain, isotherms and melt (accordion and rifting models). |
| `Geometries.ipynb` | Lithology (or other field) at chosen instants, crustal thickness and width of the margins. |
| `functions/mandyocIO.py` | Plotting functions used by the scripts (lithology legend, colorbars, P-T-t paths, `plot_property`). |

## How to install

The scripts use the `tapioca` conda environment (python, xarray, h5netcdf, matplotlib, scipy, pymp, psutil, ffmpeg):

```
git clone https://github.com/ggciag/tapioca.git
cd tapioca
conda env create -f environment.yml -y
conda activate tapioca
pip install -e .
```

Then clone this repository next to tapioca, e.g. in `~/opt/post_processing`. See
[docs/user_guide/install.md](docs/user_guide/install.md).

## How to use

1. Convert the outputs of Mandyoc to NetCDF with tapioca, in the folder of the scenario:

   ```
   julia -t 14 ~/opt/tapioca/src/postrunning/meshes2netcdf.jl ${PWD}
   julia -t 14 ~/opt/tapioca/src/postrunning/particles2netcdf.jl ${PWD} <reference_step>
   ```

2. Run the scripts **in the folder of the scenario** (they read and write there):

   ```
   python ~/opt/post_processing/track_particles_nc_dk.py
   python ~/opt/post_processing/plot_PTT_frames_dk.py
   ```

Each script has a block of settings at the top (region of the particles, layer codes, colors, which frames and
videos to make). See the [user guide](docs/index.md):

- [Workflow](docs/user_guide/workflow.md)
- [Particle tracking](docs/user_guide/tracking.md)
- [Frames and videos](docs/user_guide/frames.md)
- [Geometries notebook](docs/user_guide/geometries.md)
- [Reference of `functions/mandyocIO.py`](docs/references/mandyocIO.md)

## Contributing

Report bugs and suggest improvements in the issues. Keep the structure of the scripts (settings at the top, the
same variable names) when changing them, and test the changes on a scenario converted by tapioca (see
[AGENTS.md](AGENTS.md)).

## License

This is free software: you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in [LICENSE](LICENSE).
