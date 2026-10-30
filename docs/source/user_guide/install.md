# How to install

The scripts use the conda environment of tapioca, which has every package they need: python, numpy, xarray,
h5netcdf, scipy, matplotlib, pandas, astropy, pymp, psutil, jupyterlab and ffmpeg (for the videos and gifs). The
frame and tracking scripts also call `zip`.

```
git clone https://github.com/ggciag/tapioca.git ~/opt/tapioca
cd ~/opt/tapioca
conda env create -f environment.yml -y
conda activate tapioca
pip install -e .
```

The NetCDF files are created by the Julia scripts of tapioca (`src/postrunning/`). Install Julia and the packages
they need once:

```
julia ~/opt/tapioca/src/postrunning/install_packages.jl
```

Then clone this repository:

```
git clone <url of this repository> ~/opt/post_processing
```

The scripts import `functions/mandyocIO.py` from their own folder, so they can be run from any folder:

```
cd /path/to/scenario
python ~/opt/post_processing/plot_PTT_frames_dk.py
```

`Geometries.ipynb` imports `functions/mandyocIO.py` from `~/opt/post_processing` (or from the folder where the
notebook is run).
