# `functions/mandyocIO.py`

Plotting functions used by the scripts and the notebook. The fields come from the tapioca `.nc` files (see
[Data format](data_format.md)).

## Colors of the lithology

`COLOR_AIR`, `COLOR_SED`, `COLOR_DEC`, `COLOR_BAS`, `COLOR_UC`, `COLOR_LC`, `COLOR_LIT`, `COLOR_MLIT_UC`,
`COLOR_MLIT_LC`, `COLOR_AST` are the colors of the lithologies, and `LITHOLOGY_COLORS` has one color per Mandyoc
layer code (order of the layers in `interfaces.txt`, 0 = bottom layer) of each model type:

| Key | Codes |
| --- | --- |
| `'double_keel'` | asthenosphere, lower craton, upper craton, lithospheric mantle, lower crust 1, seed, lower crust 2, upper crust, decolement, sediments, air |
| `'rifting'` | asthenosphere, lithospheric mantle, seed, lithospheric mantle, lower crust, upper crust, air |
| `'accordion'` | asthenosphere, lithospheric mantle (3 layers), lower crust, upper crust, basaltic crust, air |

## Reading the tapioca files

- `read_tapioca_datasets(model_path, variables)`: merged dataset of the `.nc` files of `variables` (`melt`,
  `incremental_melt` and `radiogenic_heat` are read from `Phi.nc`, `dPhi.nc` and `heat.nc`).
- `read_params(model_path)`: `Nx, Nz, Lx, Lz` from `temperature.nc`.
- `read_times_steps(model_path)`: times [Myr] and steps.
- `read_lithology(model_path, step)`: layer codes of a step, `(Nzl, Nxl)` array (5 times finer mesh, first row at
  the bottom).

## Crust

- `measure_crust_from_lithology(model_path, step, crust_codes)`: top, base and thickness [km] of the crust in each
  column of the mesh, counting the cells of `crust_codes` in `lithology.nc` (each cell has `dz/5`).
- `measure_crustal_thickness(model_path, step, Nx, Nz, Lx, Lz, x_begin, x_end, crust_codes)`: maximum thickness
  between `x_begin` and `x_end`.
- `measure_margins_width(model_path, step, Nx, Nz, Lx, Lz, xl_begin, xl_end, xr_begin, xr_end, crust_codes)`:
  width of the rifted margins (each one ends where the crust is thinner than 5 km).

See [Geometries notebook](../user_guide/geometries.md).

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

- `prop`: `lithology` (`lithology.nc` + strain shading), `density`, `temperature`, `temperature_anomaly`, `pressure`,
  `strain`, `strain_rate`, `viscosity`, `radiogenic_heat` or `surface`.
- `plot_melt`: melt from `Phi` and `dPhi` in `dataset`.
- `topo_from_density`: topography from the density field (`True`) or from `surface.nc` in `model_path`.
- `plot_particles`: particles of the step from `particles_trajectories.nc` in `model_path`.
- `lithology_colors`: colors of the layer codes for `lithology` (`lithology.nc` in `model_path`): a key of
  `LITHOLOGY_COLORS` or a list with one color per code. The legend shows the different colors, without the air.
- `label_isotherms`: write the temperature of the isotherms.
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
- `trim_axs(axs, N)`: the first `N` axes of a grid of subplots as a 1D list (the others are removed).
- `change_dataset(properties, datasets)`: the fields to read to plot `properties` (e.g. `lithology` needs
  `strain` and `density`).
- `_extract_interface(z, Z, Nx, field_datai, value_to_search)`: depth where a field (e.g. the density) crosses a
  value, in each column.
