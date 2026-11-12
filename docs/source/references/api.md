# API

Docstrings of the functions of `functions/mandyocIO.py`. See the [overview](mandyocIO.md) for how they are used by
the scripts and the notebook.

---

## Colors of the lithology

`LITHOLOGY_COLORS` has one color per Mandyoc layer code (order of the layers in `interfaces.txt`, 0 = bottom layer)
of the double keel (`'double_keel'`), rifting (`'rifting'`), accordion (`'accordion'`) and SVG accordion
(`'svg_accordion'`) models, made of the
colors `COLOR_AIR`, `COLOR_SED`, `COLOR_DEC`, `COLOR_BAS`, `COLOR_UC`, `COLOR_LC`, `COLOR_LIT`, `COLOR_MLIT_UC`,
`COLOR_MLIT_LC` and `COLOR_AST`. `get_lithology_colors` returns them or, with `'auto'`, chooses the color of each layer
from its name (`LAYER_NAME_COLORS`), for any number of layers.

```{eval-rst}
.. currentmodule:: functions.mandyocIO

.. autosummary::
   :toctree: ../generated/

    get_lithology_colors
    get_layer_roles
    layer_roles
```

For scenarios with sediments or a water layer, `'auto'` uses the roles of the layers (`get_layer_roles`) and
`ROLE_COLORS`, with the sediment layers in bands of `SEDIMENT_COLORS`.

---

## Reading the tapioca files

```{eval-rst}
.. currentmodule:: functions.mandyocIO

.. autosummary::
   :toctree: ../generated/

    read_tapioca_datasets
    read_params
    read_times_steps
    read_lithology
```

---

## Crust

```{eval-rst}
.. currentmodule:: functions.mandyocIO

.. autosummary::
   :toctree: ../generated/

    measure_crust_from_lithology
    measure_crustal_thickness
    measure_margins_width
    find_far_field
```

---

## Plotting

```{eval-rst}
.. currentmodule:: functions.mandyocIO

.. autosummary::
   :toctree: ../generated/

    plot_property
    plot_lithology_legend
    add_colorbar
```

---

## Tracked particles

```{eval-rst}
.. currentmodule:: functions.mandyocIO

.. autosummary::
   :toctree: ../generated/

    plot_tracked_particles
    plot_tracked_particles_depth_coded
    plot_ptt_paths_depth_coded_frame
    plot_ptt_paths_three_particles
```

---

## Helpers

```{eval-rst}
.. currentmodule:: functions.mandyocIO

.. autosummary::
   :toctree: ../generated/

    find_nearest
    change_dataset
    trim_axs
```
