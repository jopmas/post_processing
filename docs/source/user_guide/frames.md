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
| `plot_depth_coded` | Particles colored by their initial depth in the lithospheric mantle: upper (30 km below the crust), intermediate (the next 30 km) and lower (down to the base of the lithosphere). |
| `hcrust` | Thickness of the crust (40 km if the model name has `4` and `0`, else 35 km). |
| `h_lithosphere` | Thickness of the lithosphere (crust and lithospheric mantle): from `lithosphere: ... km` in `infos_*.txt` (`get_lithosphere_thickness`) or, without it, 120 km. |
| `plot_three_particles` | Only the coldest, the hottest and an intermediate particle. |
| `plot_melt` | Melt (`Phi`) and incremental melt (`dPhi`). |
| `properties` | Fields to plot (`lithology` by default; other fields use `plot_property`). |
| `xlims`, `ylims` | Limits of the zoom [km]. |

Frames: `_output/<model>_<plot_type>_<property>_and_PTt_<step>.png`.

## `plot_PTT_frames_rft2.py`

As `plot_PTT_frames_rft.py`, but only with the particles exhumed with the lithospheric mantle. In the last step of
the tracking (`lithology.nc`), the exhumed lithospheric mantle of a column is the mantle right below the air or
below an upper crust thinner than `max_thin_crust_elements` (pockets of air inside this crust are part of it, but do
not count in its thickness). The columns where it is thicker than `min_exhumed_elements` are kept, and a particle is
kept if its final position is in this mantle.

| Variable | Meaning |
| --- | --- |
| `min_exhumed_elements` | Minimum thickness of the exhumed lithospheric mantle [elements of the mesh, dz]. |
| `max_thin_crust_elements` | Upper crust on top of the mantle thinner than this is ignored [elements of the mesh, dz]; 0: the mantle must be right below the air. |
| `thin_crust_codes` | Layers of this thin crust (`[upper_crust_code]`). |
| `exhumed_mantle_codes` | Layers of the lithospheric mantle (with the seed). |

Frames: `_output/<model>_<plot_type>_exhumed_<property>_and_PTt_<step>.png`.

## `frames_generator.py`

Accordion, rifting and double keel models (used by the [pipeline](pipeline.md) for the video). Each frame has the lithology, the strain, the isotherms and, if the model has melt, the
melt.

Input files: `strain.nc`, `lithology.nc`, `temperature.nc`, `density.nc` (and `Phi.nc`, `dPhi.nc`; `surface.nc` for the properties other than the lithology).

| Variable | Meaning |
| --- | --- |
| `lithology_colors` | `'auto'` (any number of layers, see `get_lithology_colors`), `'double_keel'`, `'rifting'`, `'accordion'`, `'svg_accordion'` or a list with one color per layer code. |
| `sediment_layers` | Scenarios with surface processes: each deposited sediment layer in a color of `SEDIMENT_COLORS`, in cycle, to show the stratigraphy (`True`, default), or all of them in the color of the sediments (`False`). Water and basalt get their own colors with `lithology_colors = 'auto'`, and the legend shows the layers of the scenario. |
| `properties` | Properties to plot (default `['lithology']`). `'lithology'` is plotted by the script; the others (`'density'`, `'radiogenic_heat'`, `'pressure'`, `'strain'`, `'strain_rate'`, `'temperature'`, `'temperature_anomaly'`, `'surface'`, `'viscosity'`) with `plot_property`, from `<property>.nc` (`heat.nc` for `radiogenic_heat`). One set of frames, video and gif per property. |
| `plot_isotherms` | Isotherms (`isotherms`). |
| `topo_from_density` | Air/crust interface of the properties plotted with `plot_property`: `False` (default) from `surface.nc` (from `density.nc` if there is no `surface.nc`), `True` from `density.nc`. |
| `plot_melt` | Melt and incremental melt: `'auto'` (default) only if the model has `Phi.nc` and `dPhi.nc` and melt in the last step, `True` or `False`. |
| `xlims`, `ylims` | Limits of the frames [km] (default: the whole model; `ylims` with 0 at the initial surface and 40 at the top of the air). |
| `fig_format` | Format of the frames (`jpeg`). |

Options of the command line, for the settings that change from a run to another. Without them, the script uses
the settings above:

| Option | Setting |
| --- | --- |
| `--prop PROP [PROP ...]` | `properties` |
| `--plot-melt`, `--no-plot-melt` | `plot_melt` (`True`, `False`) |
| `--xlims XMIN XMAX` | `xlims` [km] |
| `--ylims ZMIN ZMAX` | `ylims` [km] |

```bash
python ~/opt/post_processing/frames_generator.py                                   # settings of the script
python ~/opt/post_processing/frames_generator.py --prop lithology temperature --no-plot-melt
python ~/opt/post_processing/frames_generator.py --xlims 600 1200 --ylims -150 40  # zoom
python ~/opt/post_processing/frames_generator.py --help
```

Frames: `_output/<model>_<property>[_MeltFrac]_<step>.jpeg`.
