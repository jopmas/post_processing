# post_processing

**Scripts and notebooks to post-process [Mandyoc](https://github.com/ggciag/mandyoc) models from the NetCDF files
created by [tapIOca](https://github.com/ggciag/tapioca).**

## About

The scripts read only the `.nc` files created by the Julia scripts of tapIOca (`meshes2netcdf.jl` and
`particles2netcdf.jl`): particle tracking (P-T-t paths), frames and videos of the lithology, and measurements of
the crust (crustal thickness and width of the rifted margins).

## User guide

- [How to install](user_guide/install.md)
- [Workflow](user_guide/workflow.md): from the outputs of Mandyoc to the figures.
- [Pipeline](user_guide/pipeline.md): `pipeline.sh`, from the zip of a scenario to the snapshots and the video.
- [Particle tracking](user_guide/tracking.md): `track_particles_nc_dk.py` and `track_particles_nc_rft.py`.
- [Frames and videos](user_guide/frames.md): `plot_PTT_frames_dk.py`, `plot_PTT_frames_rft.py` and
  `frames_generator.py`.
- [Geometries notebook](user_guide/geometries.md): fields at chosen instants, crustal thickness and width of the
  margins.

## References

- [Data format](references/data_format.md): structure of the tapioca `.nc` files and of `_track_xzPT_all_steps.nc`.
- [`functions/mandyocIO.py`](references/mandyocIO.md): overview of the functions.
- [API](references/api.md): docstrings of the functions of `functions/mandyocIO.py`.

## Building this documentation

With the `tapioca` environment (it has Sphinx, myst-parser and sphinx_rtd_theme), in the `docs` folder:

```
make html
```

The pages are written in `docs/build/html` (open `docs/build/html/index.html`).

## License

This is free software: you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in
[LICENSE](https://github.com/jopmas/post_processing/blob/main/LICENSE).

```{toctree}
:maxdepth: 2
:hidden:

user_guide/index
references/index
```
