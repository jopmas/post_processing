########################################################################
# BEFORE RUNNING THIS SCRIPT:                                          #
#  - CREATE THE .nc FILES WITH tapioca's meshes2netcdf.jl              #
#  - CREATE _gz_evolution.nc WITH gz_evolution.py                      #
#  - RUN THIS SCRIPT IN THE SAME FOLDER WHERE THE .nc FILES ARE        #
#                                                                      #
# Frames of the g_z profile, the topography and the lithology, strain  #
# and isotherms (as in frames_generator.py), and the video. The frames #
# and the video are in _output, named <model>_gz_lithology*.           #
########################################################################

import gc
import numpy as np
import os
import subprocess
import sys
import xarray as xr
import multiprocessing #needed to run pymp in macOS
multiprocessing.set_start_method('fork') #needed to run pymp in macOS arm architectures
import pymp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

matplotlib.use('agg')

sys.path.append(os.path.dirname(os.path.abspath(__file__))) # post_processing folder, with functions/mandyocIO.py
from functions.mandyocIO import plot_lithology_legend, get_lithology_colors

model_path = os.getcwd() # Get local file
model_name = model_path.split('/')[-1]
output_path = '_output'

if not os.path.isdir(output_path):
    os.makedirs(output_path)

###########################################################################
# Settings                                                                #
###########################################################################

xlims = None #[km]. None: the whole model
ylims = [-150, 40] #[km], limits of the lithology. Whole model: None

gz_step = 50.0 #mGal, the limits of g_z are multiples of it around the minimum and maximum of the whole evolution
topo_step = 1.0 #km, the same for the topography

isotherms = [500, 600, 700, 800, 900, 1300] #oC

stride = 1 #uses one of each stride time steps of _gz_evolution.nc

make_videos = True
fps = 12 #frames per second of the video
zip_files = True #moves the frames of this script to <model>_gz_imgs.zip

fig_format = 'jpeg'
dpi = 300

lithology_colors = 'auto' #see get_lithology_colors of functions/mandyocIO.py

cr = 255.
color_uc = (228./cr,156./cr,124./cr)
color_lc = (240./cr,209./cr,188./cr)
color_lit = (155./cr,194./cr,155./cr)
color_ast = (207./cr,226./cr,205./cr)

###########################################################################
# Reading the mesh and g_z                                                #
###########################################################################

#tapioca files: data (time, z, x), x from 0 to Lx and z from -Lz to 0 [m]
with xr.open_dataset(f"{model_path}/strain.nc", engine="h5netcdf") as mesh:
    Nx = int(mesh.sizes["x"])
    Nz = int(mesh.sizes["z"])
    Lx = float(mesh.x.values[-1] - mesh.x.values[0])
    Lz = float(mesh.z.values[-1] - mesh.z.values[0])
    steps_mesh = mesh.step.values

x = np.linspace(0, Lx/1000.0, Nx)
z = np.linspace(-Lz/1000.0, 0, Nz)

extent = (0, Lx/1000, -Lz/1000+40, 40)

if(xlims is None):
    xlims = [0, Lx/1000]
if(ylims is None):
    ylims = [-Lz/1000+40, 40]

#Reading _gz_evolution.nc before the parallel part
with xr.open_dataset(f"{model_path}/_gz_evolution.nc", engine="h5netcdf") as gz_dataset:
    gz = gz_dataset.gz.values
    x_gz = gz_dataset.x.values/1000.0
    times_gz = gz_dataset.time.values
    steps_gz = gz_dataset.step.values
    surface = gz_dataset.surface.values/1000.0
    surface_sunk = gz_dataset.surface_sunk.values/1000.0
    sea_level = float(gz_dataset.attrs["sea_level"])/1000.0
    obs_height = float(gz_dataset.attrs["obs_height"])
    water = bool(gz_dataset.attrs["water"])

#The topography is relative to the sea level
topo_ylims = [np.floor((surface_sunk.min() - sea_level)/topo_step)*topo_step, np.ceil((surface.max() - sea_level)/topo_step)*topo_step]

#The same limits in all frames
gz_ylims = [np.floor(gz.min()/gz_step)*gz_step, np.ceil(gz.max()/gz_step)*gz_step]

#index of each step of _gz_evolution.nc inside the mesh files
idx_mesh = np.searchsorted(steps_mesh, steps_gz)
if(np.any(idx_mesh >= steps_mesh.size) or not np.array_equal(steps_mesh[np.minimum(idx_mesh, steps_mesh.size-1)], steps_gz)):
    raise ValueError("_gz_evolution.nc has steps that are not in the mesh files.")

idx_frames = np.arange(0, steps_gz.size, stride)
n_frames = idx_frames.size

colors = get_lithology_colors(model_path, lithology_colors)
cmap = ListedColormap(colors)
np.seterr(divide='ignore')

print(f"Model name: {model_name}")
print(f"Frames: {n_frames}, g_z limits: {gz_ylims} mGal")
print("Generating frames...")
sys.stdout.flush() #otherwise, when the output goes to a file, each forked process writes the buffered text again

###########################################################################
# Frames                                                                  #
###########################################################################

with pymp.Parallel() as p:
    #Opening the files inside the parallel region: each process has its own handles
    lithology = xr.open_dataset("lithology.nc", cache=False, engine="h5netcdf")["lithology"]
    strain = xr.open_dataset("strain.nc", cache=False, engine="h5netcdf")["strain"]
    temperature = xr.open_dataset("temperature.nc", cache=False, engine="h5netcdf")["temperature"]
    density = xr.open_dataset("density.nc", cache=False, engine="h5netcdf")["density"] #used to remove the values in the air

    for frame in p.range(0, n_frames):
        k = idx_frames[frame] #index in _gz_evolution.nc
        i = idx_mesh[k] #index in the mesh files

        fig, (ax_gz, ax_topo, ax_litho) = plt.subplots(3, 1, figsize=(12, 7.5), sharex=True, height_ratios=[1, 0.7, 1.4], constrained_layout=True)

        fsize = 14

        #g_z profile
        ax_gz.axhline(0, color='0.6', lw=0.8, ls='--')
        ax_gz.plot(x_gz, gz[k], color='k', lw=1.5)
        ax_gz.set_ylim(gz_ylims)
        ax_gz.set_ylabel(r'$g_z$ [mGal]', fontsize=fsize)
        ax_gz.tick_params(axis='both', labelsize=fsize)
        ax_gz.grid(True, alpha=0.3)
        ax_gz.text(0.01, 1.05, f'{model_name}', fontsize=fsize, transform=ax_gz.transAxes)
        ax_gz.text(0.99, 1.05, f'Time = {times_gz[k]:.2f} Myr', fontsize=fsize, ha='right', transform=ax_gz.transAxes)
        ax_gz.text(0.01, 0.06, f'{obs_height:.0f} m above the water or the rock', fontsize=fsize-3, transform=ax_gz.transAxes)

        #Topography
        ax_topo.axhline(0, color='0.6', lw=0.8, ls='--')
        ax_topo.plot(x_gz, surface[k] - sea_level, color='tab:blue', lw=1.5, label='Model')
        if(water):
            ax_topo.plot(x_gz, surface_sunk[k] - sea_level, color='k', lw=1.0, ls=':', label='With the water load')
            ax_topo.legend(loc='lower right', fontsize=fsize-4, frameon=False)
        ax_topo.set_ylim(topo_ylims)
        ax_topo.set_ylabel('Topography [km]', fontsize=fsize)
        ax_topo.tick_params(axis='both', labelsize=fsize)
        ax_topo.grid(True, alpha=0.3)

        #Lithology, strain and isotherms, as in frames_generator.py
        air = density.isel(time=i).values < 1200.0 #tapioca keeps the values in the air: set them to zero there

        ax_litho.imshow(lithology.isel(time=i).values[::-1, :], aspect='auto', extent=extent, cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, alpha=1.0)
        ax_litho.imshow(np.log10(np.where(air, 0.0, strain.isel(time=i).values)[::-1, :]), extent=extent, cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2, aspect='auto')
        ax_litho.contour(x, z+40, temperature.isel(time=i).values, levels=isotherms, colors='r', linewidths=1.0)
        ax_litho.set_xlim(xlims)
        ax_litho.set_ylim(ylims)

        plot_lithology_legend(ax_litho, colors=[color_uc, color_lc, color_lit, color_ast], vmin=-0.5, vmax=0.9, alpha=0.2,
                              position=(0.98, 0.25), loc='lower right')

        ax_litho.set_xlabel('Distance [km]', fontsize=fsize)
        ax_litho.set_ylabel('Depth [km]', fontsize=fsize)
        ax_litho.tick_params(axis='both', labelsize=fsize)

        figname = f"{model_name}_gz_lithology_{str(int(steps_gz[k])).zfill(6)}.{fig_format}"
        fig.savefig(f"{output_path}/{figname}", dpi=dpi)
        plt.close(fig)

        del fig, ax_gz, ax_topo, ax_litho, air
        gc.collect()

        p.print(f"Thread {p.thread_num}: time {times_gz[k]:.2f} Myr, step {steps_gz[k]} ({frame+1}/{n_frames})")

print("Done!")

###########################################################################
# Video                                                                   #
###########################################################################

basename = f"{model_name}_gz_lithology"

if(make_videos):
    print("Generating video...")
    videoname = f"{model_path}/{output_path}/{basename}"
    comand = f"ffmpeg -r {fps} -f image2 -pattern_type glob -i \"{videoname}_*.{fig_format}\" -vcodec libx264 -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -r 24 -y -an -crf 25 -pix_fmt yuv420p \"{videoname}.mp4\""
    #Without capture_output, the error message of ffmpeg is shown in the terminal
    subprocess.run(comand, shell=True, check=True)
    print(f"\t{output_path}/{basename}.mp4 written")

if(zip_files):
    #-m moves the files into the zip: each frame is deleted only after it is in the zip
    print("Zipping frames...")
    subprocess.run(f"zip -qm {model_name}_gz_imgs.zip {basename}_*.{fig_format}", shell=True, check=True, cwd=output_path)
    print(f"\t{output_path}/{model_name}_gz_imgs.zip written")
