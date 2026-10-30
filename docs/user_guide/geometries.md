# Geometries notebook

`Geometries.ipynb` plots a Mandyoc field (lithology, density, temperature, strain, ...) at chosen instants, one
panel per instant, and measures the crust: the maximum crustal thickness in a region and the width of the rifted
margins. It reads `temperature.nc`, `density.nc`, `strain.nc`, `lithology.nc` and, if it exists, `surface.nc`.

## Settings (cell "Set scenario path")

| Variable | Meaning |
| --- | --- |
| `scenario` | Path of the scenario, without `/Users/<user>` (added by `machine_path`). |
| `prop` | Field to plot (`lithology` by default). |
| `crust_codes` | Mandyoc layer codes of the crust (order of the layers in `interfaces.txt`, 0 = bottom). |
| `instants` | Instants to plot [Myr]; the nearest steps are used. |
| `plot_thickness` | Measure the crustal thickness and the width of the margins. |

`crust_codes` depends on the model:

- double keel: `[4, 5, 6, 7, 8, 9]` (lower crust 1, seed, lower crust 2, upper crust, decolement, sediments);
- accordion (AR) and rifting: `[4, 5]` (lower crust, upper crust), or `[4, 5, 6]` to include the basaltic crust of
  the accordion models.

The instants and the regions used to measure the thickness and the margins (`x_begin`, `x_end`, `xl_begin`,
`xl_end`, `xr_begin`, `xr_end`) are set in the plotting loop and must be adapted to each model.

## Crustal thickness

The crust is measured from `lithology.nc`, counting the cells of `crust_codes` in each column, as the thickness of
the oceanic crust (`thickness_oceanic_crust`) in `frames_generator_litho_melt_thickness.py`:

```
thickness = (number of crust cells in the column) * dz / 5
```

since the lithology mesh is 5 times finer than the mesh of the other fields (one column of each 5 is on the x of
the mesh). The top and the base of the crust are the highest and the lowest crust cells of the column.

- `measure_crust_from_lithology(fpath, step, crust_codes)`: top, base and thickness of the crust in every column.
- `measure_crustal_thickness(...)`: maximum thickness between `x_begin` and `x_end`.
- `measure_margins_width(...)`: each margin ends where the crust is thinner than 5 km (the first column in
  `[xl_begin, xl_end]` and the last one in `[xr_begin, xr_end]`).

The margins are only found when the crust is thinner than 5 km in these regions (e.g. after the rifting of the
accordion models); otherwise `measure_margins_width` raises an `IndexError`.

The last column of `lithology.nc` (x = Lx) has no crust (an artifact of the border): avoid it in the regions.
