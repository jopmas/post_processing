# post_processing

**Scripts and notebooks to post-process [Mandyoc](https://github.com/ggciag/mandyoc) models from the NetCDF files
created by [tapIOca](https://github.com/ggciag/tapioca).**

## User guide

- [How to install](user_guide/install.md)
- [Workflow](user_guide/workflow.md): from the outputs of Mandyoc to the figures.
- [Particle tracking](user_guide/tracking.md): `track_particles_nc_dk.py` and `track_particles_nc_rft.py`.
- [Frames and videos](user_guide/frames.md): `plot_PTT_frames_dk.py`, `plot_PTT_frames_rft.py` and
  `frames_generator.py`.
- [Geometries notebook](user_guide/geometries.md): fields at chosen instants, crustal thickness and width of the
  margins.

## References

- [Data format](references/data_format.md): structure of the tapioca `.nc` files and of `_track_xzPT_all_steps.nc`.
- [`functions/mandyocIO.py`](references/mandyocIO.md): plotting functions.
