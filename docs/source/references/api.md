# API

Docstrings of the functions of `functions/mandyocIO.py`. See the [overview](mandyocIO.md) for how they are used by
the scripts and the notebook.

---

## Colors of the lithology

`LITHOLOGY_COLORS` has one color per Mandyoc layer code (order of the layers in `interfaces.txt`, 0 = bottom layer)
of the double keel (`'double_keel'`), rifting (`'rifting'`) and accordion (`'accordion'`) models, made of the
colors `COLOR_AIR`, `COLOR_SED`, `COLOR_DEC`, `COLOR_BAS`, `COLOR_UC`, `COLOR_LC`, `COLOR_LIT`, `COLOR_MLIT_UC`,
`COLOR_MLIT_LC` and `COLOR_AST`.

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
