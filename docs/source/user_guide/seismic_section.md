# Synthetic seismic section

`seismic_section.py` (and the same code as a notebook, `seismic_section.ipynb`) builds a synthetic seismic
section of a rifted margin from a Mandyoc scenario. It works with scenarios with and without sediments
(surface processes), extracted basaltic magma and a water layer. It converts the lithology, temperature and strain
into density and P-wave velocity, computes a high-frequency convolutional
section of a wide window, and simulates the wave propagation in a zoom with
[tremelique](https://www.fatiando.org/tremelique), followed by a reverse-time migration (RTM). The notebook explains
each step and the reason for each choice in its text cells.

It needs [tremelique](https://www.fatiando.org/tremelique) installed in the `tapioca` environment.

## How to run

The script runs in the folder of the scenario, like the other scripts:

```
cd /path/to/scenario
python ~/opt/post_processing/seismic_section.py
```

The notebook reads the scenario from `DATA` (cell "Parameters"). It reads `lithology.nc`, `temperature.nc`,
`strain.nc`, `surface.nc`, `interfaces.txt` and `param.txt`, and `SP_vc.txt` and `infos_*.txt` when the scenario
has them; a lithology on a fine mesh (e.g. 200 m) and temperature and strain on a 1 km mesh are needed to resolve
the faults and the serpentinization. Set `TIME` and the windows (`X0`, `X1`, `XZ0`, `XZ1`, `CANDIDATES`) for each
scenario.

## Layer roles, sea level, sediments and basalt

The seismic properties depend on what each Mandyoc layer is, and the layers change between scenarios. With
`LAYER_ROLES = "auto"` the role of each layer code (`asthenosphere`, `lithospheric mantle`, `lower crust`,
`upper crust`, `basalt`, `sediments`, `water`, `air`) is found, in this order of priority, from:

1. the names of the layers in `infos_*.txt`, when the scenario has it;
2. the sediment codes of `SP_vc.txt` (one code per time interval of the surface processes);
3. the densities of `interfaces.txt`: < 500 kg/m³ air, < 1200 water, < 2550 sediments, ≥ 3200 mantle (code 0 is
   the asthenosphere), 2750–3200 lower crust or, **above an upper-crust layer, basalt** (extracted magma sits on
   top of the crust), otherwise upper crust.

The roles are printed; if they are wrong, give `LAYER_ROLES` by hand, e.g.
`{"asthenosphere": [0], "lithospheric mantle": [1, 2, 3], "lower crust": [4], "upper crust": [5],
"basalt": [6], "sediments": list(range(7, 40)), "water": [40], "air": [41]}`.

With `SEA_LEVEL = "auto"`, sea level is the top of the water layer at `TIME` when the scenario has one (the whole
model may move up or down during the run), otherwise Mandyoc's `sea_level` above the mean initial surface. Water
fills everything between the seafloor (`surface.nc`) and sea level; land above sea level is cut off.

- **Sediments**: Mandyoc gives them a constant density, but their seismic properties come from compaction. The
  porosity decreases with burial (Athy's law, `SEDIMENTS["phi0"]`, `["decay"]`), Vp follows Wyllie's time average of
  water and matrix and the density follows the porosity, so Vp grows from ~2.2 km/s at the seafloor to ~4.5 km/s at
  depth. Each deposited layer (each sediment code) gets its own surface porosity and matrix velocity, which makes
  the layered reflections, onlaps and wedges. Faults are not added in the sediments (`FAULTS["in_sediments"]`).
- **Basalt** (extracted magma): Vp `VP0["basalt"]` with the temperature correction, density from `interfaces.txt`,
  short-scale heterogeneity.
- **Top of the basement**: the first crystalline rock (crust, basalt or mantle) below the sediments;
  serpentinization is measured from it.

Slow sediments lower the frequency of the tremelique zoom (the slowest rock sets the grid budget); a smaller `DXZ`
compensates, at a cost that grows with `1/DXZ³`.

## Steps

Both notebooks start with a flowchart of these steps, drawn by `plot_flowchart` of
`functions/seismic_section_flowchart.py`.


1. **Load** one time step (`TIME`) of the Mandyoc fields, cropped around the window; find the roles of the layers
   and sea level, and reference the depths to sea level.
2. **Fine grid**: each role is smoothed and the largest one wins, so the contacts stay sharp without the
   staircase of the Mandyoc mesh; the age (code) of each sediment layer is interpolated the same way; water fills
   everything between the seafloor (`surface.nc`) and sea level.
3. **Properties**: density `rho0 (1 - alpha T)` and Vp per role with a temperature correction for the crust,
   basalt and mantle; compaction and one impedance per deposited layer for the sediments; plus
   - serpentinization of the cold (< `SERP_TMAX`), shallow (< `SERP_DEPTH` below the basement) mantle;
   - small-scale heterogeneity inside each layer (von Kármán random media, laminated in the lower crust), which
     does not change the geometry;
   - faults from the localized strain (strain minus its regional mean), in the brittle domain: a damage zone
     (lower Vp) and a fabric laminated along the shear band.
4. **Convolutional section** of the whole window (`X0`–`X1`) at 25 Hz: cheap, but 1D (no diffractions).
5. **Zoom** (`XZ0`–`XZ1`, down to `ZZMAX`) with `tremelique.Acoustic`: exploding-reflector modelling (zero-offset
   section from a single simulation) with receivers on the seafloor, and RTM with the same solver. The sections
   are shown in two-way time with AGC, in grey scale.

## Settings

| Variable | Meaning |
| --- | --- |
| `TIME` | Model time [Myr]; the nearest output is used. |
| `X0`, `X1`, `ZMAX`, `DX` | Window, maximum depth and grid spacing of the convolutional section [m]. |
| `LAYER_ROLES` | Role of each Mandyoc layer code: `"auto"` or a dict `{role: [codes]}`. |
| `SEA_LEVEL` | `"auto"` (top of the water layer, or initial surface + `sea_level`) or a z in Mandyoc coordinates [m]. |
| `VP0` | Vp at 0 °C of the crust, basalt and mantle [m/s]. |
| `SEDIMENTS` | Compaction (surface porosity, decay length), matrix velocity, grain density and the spread between layers. |
| `SERP_DEPTH`, `SERP_TMAX` | Depth below the basement and maximum temperature of the serpentinization. |
| `HETEROGENEITY`, `SERP_HETEROGENEITY` | Strength and correlation lengths of the random heterogeneity. |
| `FAULTS`, `BRITTLE_T` | Faults from the strain: high-pass length, damage, fabric, brittle temperatures. |
| `CANDIDATES` | Windows `(x start, x end, TWT start)` of the 4 s zooms of the convolutional section; a TWT start of `None` starts just above the seafloor. |
| `XZ0`, `XZ1`, `ZZMAX`, `DXZ` | Window, maximum depth and grid spacing of the tremelique zoom [m]. |
| `ZOOM_T0`, `ZOOM_TLEN` | TWT window of the zoom figures [s]; `ZOOM_T0 = None` starts just above the seafloor. |

The run time is in the two tremelique simulations of the zoom: it grows with the size of the zoom and with
`1/DXZ³` (minutes for 70 km × 32 km, about 40 min for 150 km × 36 km with `DXZ = 20` m).

## Outputs

Saved in the `_output` folder of the scenario:

| File | Content |
| --- | --- |
| `seismic_section_model.png` | Lithology, serpentinization, faults and Vp of the whole window. |
| `seismic_section_convolutional.png` | Convolutional section of the whole window. |
| `seismic_section_convolutional_windows.png` | 4 s zooms of the convolutional section (`CANDIDATES`). |
| `seismic_section_zoom_<window>_sections.png` | Lithology, unmigrated, migrated and convolutional sections of the zoom, in TWT. |
| `seismic_section_zoom_<window>_final.png` | Lithology of the whole window (depth) and Vp, RTM and convolutional sections of the zoom. |
| `seismic_section_zoom_<window>.nc` | Results of the zoom (see below). |

## Results of the zoom (`seismic_section_zoom_<window>.nc`)

The results of the tremelique simulations are saved in a NetCDF file, readable with `xarray.open_dataset`:

| Variable | Dims | Content |
| --- | --- | --- |
| `image` | `(depth, x)` | RTM image in depth (below sea level) [m]. |
| `unmigrated` | `(twt, x)` | Unmigrated zero-offset section, in the TWT window of the figures. |
| `migrated` | `(twt, x)` | RTM section in TWT. |
| `convolutional` | `(twt_conv, x_conv)` | Convolutional section of the zoom. |

The attributes keep the time step `dt`, the maximum frequency `F_CUT`, the scenario, the model time and
`parameters`, a JSON with all the parameters that change the results. Before running the simulations, the
notebook looks for this file: if it exists and its `parameters` are the current ones, the results are loaded
instead (the notebook then runs in a few minutes); otherwise the simulations run and the file is written again.
Set `RERUN = True` to force the simulations.

## Multithreaded version

`seismic_section_parallel.py` and `seismic_section_parallel.ipynb` are the same pipeline, with the tremelique time
step split between the cores of the computer. In tremelique the time step of `Acoustic` (`timestep_esg`,
`apply_damping` and `nonreflexive_bc`) is compiled with numba but runs on one core; the time loop cannot be
parallelised (each step needs the two previous ones), but within a step the cells are updated independently.
`functions/tremelique_parallel.py` has copies of these kernels with `parallel=True` and the loops over the rows
split between threads (`numba.prange`). Each cell gets the same operations, so the results are identical to the
serial version, bit by bit.

| | Serial | Multithreaded (12 threads) |
| --- | --- | --- |
| Zoom 750–900 km × 36 km (`DXZ = 20` m): exploding reflector + RTM | ~36 min | ~6 min |
| Whole notebook | ~38 min | ~9 min |

The number of threads is numba's (all the cores by default; change it with `numba.set_num_threads` or the
`NUMBA_NUM_THREADS` environment variable). The outputs have the prefix `seismic_section_parallel_` (figures and
`seismic_section_parallel_zoom_<window>.nc`), so they do not overwrite those of the serial version.

`functions/tremelique_parallel.py` also has `ParallelAcoustic`, a drop-in replacement for `tremelique.Acoustic`
(same arguments, HDF5 cache, `snapshot` and `animate`) for other uses of tremelique:

```python
import sys
import tremelique
sys.path.append("/path/to/post_processing")
from functions.tremelique_parallel import ParallelAcoustic

sim = ParallelAcoustic(velocity, density, spacing)
sim.add_point_source((0, nx // 2), tremelique.RickerWavelet(1, 30))
sim.run(1000)
```

## Limitations

- The properties of the sediments come from a compaction law, and the contrasts between layers are random.
- Acoustic (P waves only) and zero offset: no converted waves, multiples and amplitudes versus offset are not
  correct.
- The heterogeneity and the fault fabric are statistical: they give the texture of real data, not specific
  reflectors.
- The frequency of the zoom (~10 Hz peak with `DXZ = 20` m) is lower than in field data.
