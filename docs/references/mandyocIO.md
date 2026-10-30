# `functions/mandyocIO.py`

Plotting functions used by the scripts and the notebook. The fields come from the tapioca `.nc` files (see
[Data format](data_format.md)).

## Legend and colorbars

Both have a fixed size in inches, so they look the same whatever the size of the figure. Only the position is
relative to the axis: `position` is the point (x, y) in fraction of `ax` where the corner `loc` of the inset is
placed.

### `plot_lithology_legend(ax, colors, vmin=-0.5, vmax=0.7, alpha=0.2, position=(0.98, 0.08), loc='lower right', width=0.74, height=0.66, fsize=10)`

Lithology/strain legend: one band per color (`colors` from top to bottom) shaded by the strain from `vmin` to
`vmax` (log10) with `alpha`, which must be the values used to plot the strain. Returns the axis of the box.

```python
plot_lithology_legend(ax, colors=[color_uc, color_lc, color_lit, color_ast], vmin=-0.5, vmax=0.9,
                      position=(0.98, 0.20), loc='lower right')
```

### `add_colorbar(fig, ax, im, label, position=(0.02, 0.08), loc='lower left', length=1.8, thickness=0.11, orientation='horizontal', fsize=10, fmt=None, ticks=None)`

Colorbar of `im` (output of `imshow`, `contourf`, ...) inside `ax`. Returns the colorbar.

## Fields

### `plot_property(dataset, prop, xlims, ylims, model_path, fig, ax, ...)`

Plots a field of one time step (`dataset` with the variables needed, e.g. merged from `strain.nc`, `density.nc`,
`temperature.nc`), with isotherms and the topography filled in white:

- `prop`: `lithology` (density + strain shading), `density`, `temperature`, `temperature_anomaly`, `pressure`,
  `strain`, `strain_rate`, `viscosity`, `radiogenic_heat` or `surface`.
- `plot_melt`: melt from `Phi` and `dPhi` in `dataset`.
- `topo_from_density`: topography from the density field (`True`) or from `surface.nc` in `model_path`.
- `plot_particles`: particles of the step from `particles_trajectories.nc` in `model_path`.
- `bbox_to_anchor`: the legend of `lithology` has its lower right corner at the lower right corner of this box
  (x, y, width, height in fraction of `ax`).

The colorbars use `add_colorbar` and the legend uses `plot_lithology_legend`.

## Tracked particles

They read `_track_xzPT_all_steps.nc` (`trackdataset`) and plot the time index `i`:

- `plot_tracked_particles(trackdataset, ax, i, ...)` and `plot_ptt_paths_three_particles(trackdataset, ax, i, current_time, ...)`:
  the coldest, the hottest and an intermediate particle of the lithospheric mantle.
- `plot_tracked_particles_depth_coded(trackdataset, ax, i, hcrust, ...)` and
  `plot_ptt_paths_depth_coded_frame(trackdataset, ax, i, current_time, hcrust, ...)`: particles colored by their
  initial depth in the lithospheric mantle (upper, intermediate, lower).

## Helpers

- `find_nearest(array, value)`: index of the element of `array` nearest to `value`.
- `change_dataset(properties, datasets)`: the fields to read to plot `properties` (e.g. `lithology` needs
  `strain` and `density`).
- `_extract_interface(z, Z, Nx, field_datai, value_to_search)`: depth where a field (e.g. the density) crosses a
  value, in each column.
