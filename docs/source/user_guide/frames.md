# Frames and videos

The frame scripts write the frames in `_output/`, then (if `make_videos`, `make_gifs`, `zip_files` are `True`) the
videos (`ffmpeg`), the gifs and zip files of the frames, videos and gifs. With `zip_files = True` the frames are
removed after they are zipped.

The lithology is plotted from `lithology.nc` with one color per Mandyoc layer code (the order of the layers in
`interfaces.txt`, 0 = bottom layer). Each script builds a `lithology_colors` dict from its layer code variables:
change it when the layers of the model are different. The strain is shaded with `log10(strain)` (zero in the air,
where the density is lower than 1200 kg/m3). The lithology/strain legend and the colorbars have a fixed size in
inches (see [`plot_lithology_legend` and `add_colorbar`](../references/mandyocIO.md)).

## `plot_PTT_frames_dk.py`

Double keel models. Each frame has the lithology, the strain, the isotherms and the tracked particles, and the
P-T-t paths over the metamorphic facies. The particles are colored by their maximum temperature and pressure
(wedge or subducted, below or above 900 oC); particles that never reach 400 oC are not plotted.

Input files: `strain.nc`, `lithology.nc`, `temperature.nc`, `density.nc` and `_track_xzPT_all_steps.nc` (from
`track_particles_nc_dk.py`).

| Variable | Meaning |
| --- | --- |
| `double_keel` | Layer codes and colors of the double keel model (11 layers) or of the homogeneous model (9). |
| `start`, `end`, `step` | Time indices of the frames (`step = 4`: one frame every 4 steps). |
| `strain_vmin`, `strain_vmax`, `strain_alpha` | Strain shading (log10). |
| `plot_colorbar` | Plot the lithology/strain legend. |

Frames: `_output/<model>_lithology_and_PTt_temperature_coded_<step>.png`.

## `plot_PTT_frames_rft.py`

Rifting models. Each frame has a zoom of the lithology with the tracked particles, the P-T-t paths and the whole
model below.

Input files: `strain.nc`, `temperature.nc`, `density.nc`, `lithology.nc` (and `Phi.nc`, `dPhi.nc` with
`plot_melt`) and `_track_xzPT_all_steps.nc` (from `track_particles_nc_rft.py`).

| Variable | Meaning |
| --- | --- |
| `plot_depth_coded` | Particles colored by their initial depth in the lithospheric mantle (upper, intermediate, lower). |
| `plot_three_particles` | Only the coldest, the hottest and an intermediate particle. |
| `plot_melt` | Melt (`Phi`) and incremental melt (`dPhi`). |
| `properties` | Fields to plot (`lithology` by default; other fields use `plot_property`). |
| `xlims`, `ylims` | Limits of the zoom [km]. |

Frames: `_output/<model>_<plot_type>_<property>_and_PTt_<step>.png`.

## `frames_generator.py`

Accordion, rifting and double keel models (used by the [pipeline](pipeline.md) for the video). Each frame has the lithology, the strain, the isotherms and, with `plot_melt`, the
melt.

Input files: `strain.nc`, `lithology.nc`, `temperature.nc`, `density.nc` (and `Phi.nc`, `dPhi.nc`).

| Variable | Meaning |
| --- | --- |
| `lithology_colors` | `'auto'` (from the number of layers in `lithology.nc`: 11 double keel, 7 rifting, 8 accordion), `'double_keel'`, `'rifting'`, `'accordion'` or a list with one color per layer code. |
| `plot_melt` | Melt and incremental melt. |
| `fig_format` | Format of the frames (`jpeg`). |

Frames: `_output/<model>_lithology_<step>.jpeg`.
