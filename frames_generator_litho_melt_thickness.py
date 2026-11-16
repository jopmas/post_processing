#!/usr/bin/env python
# coding: utf-8

#Frames, video and gif of a Mandyoc model converted by tapioca (run in the folder of the scenario) with two panels:
#  - top: thickness of the oceanic crust (basalt layers) along x;
#  - bottom: lithology shaded by the strain, isotherms and melt, as in frames_generator.py (same colors, strain legend
#    and legend of the layers below the figure).
#e.g. python frames_generator_litho_melt_thickness.py --xlims 1600 2200 --ylims -150 40

import os
import gc
import sys
import argparse
import multiprocessing #needed to run pymp in mac
import psutil
multiprocessing.set_start_method('fork') #needed to run pymp in mac
# multiprocessing.set_start_method('spawn')
import pymp
import subprocess
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap
import matplotlib
import matplotlib.pyplot as plt


matplotlib.use('agg')

path_to_functions = os.path.dirname(os.path.abspath(__file__)) # post_processing folder, with functions/mandyocIO.py
sys.path.append(os.path.abspath(path_to_functions))

if '' in sys.path:
    sys.path.remove('')
from functions.mandyocIO import plot_strain_legend, get_lithology_colors, get_layer_roles, COLOR_SED, get_layers_legend, layers_legend_layout, plot_layers_legend

####################################################################################################################################
#Options of the command line: the settings below are used for the options that are not given
parser = argparse.ArgumentParser(description="Frames, video and gif of the lithology and of the thickness of the oceanic crust of a Mandyoc model converted by tapioca (run in the folder of the scenario).")
parser.add_argument('--plot-melt', action=argparse.BooleanOptionalAction, default=None,
                    help="plot (--plot-melt) or not (--no-plot-melt) the melt (default: plot_melt of the script, 'auto': only if the model has melt)")
parser.add_argument('--xlims', nargs=2, type=float, metavar=('XMIN', 'XMAX'), help="horizontal limits [km] (default: xlims of the script)")
parser.add_argument('--ylims', nargs=2, type=float, metavar=('ZMIN', 'ZMAX'), help="vertical limits [km], 0 at the initial surface and 40 at the top of the air (default: ylims of the script)")
args = parser.parse_args()

model_path = os.getcwd() # Get local file
model_name = model_path.split('/')[-1]
output_path = '_output'
print(f"Model name: {model_name}\n")
print(f"Model path: {model_path}\n")
print(f"Output path: {output_path}\n")

if not os.path.isdir(output_path):
    os.makedirs(output_path)

plot_isotherms = True
# plot_isotherms = False
plot_melt = 'auto' #only if the model has melt (Phi.nc and dPhi.nc, with melt in the last step), see model_has_melt
# plot_melt = True
# plot_melt = False

def model_has_melt(model_path):
    '''True if the model has Phi.nc and dPhi.nc (tapioca) and melt (Phi > 0 out of the air) in the last step'''
    if not (os.path.isfile(f"{model_path}/Phi.nc") and os.path.isfile(f"{model_path}/dPhi.nc")):
        return False
    with xr.open_dataset(f"{model_path}/Phi.nc") as melt_file, xr.open_dataset(f"{model_path}/density.nc") as density_file:
        melt_last = melt_file["Phi"].isel(time=-1).values
        air = density_file["density"].isel(time=-1).values < 1200.0 #tapioca keeps the values in the air
    return bool(np.nanmax(np.where(air, 0.0, melt_last)) > 0.0)

if(args.plot_melt is not None):
    plot_melt = args.plot_melt
if(plot_melt == 'auto'):
    plot_melt = model_has_melt(model_path)
    print(f"Melt: {'found, plotting it' if plot_melt else 'not found (no Phi.nc/dPhi.nc or no melt), not plotting it'}")

#########################################
# Get domain informations               #
#########################################

#tapioca files: data (time, z, x), x from 0 to Lx and z from -Lz to 0 [m]
with xr.open_dataset(f"{model_path}/strain.nc") as mesh:
    Nx = int(mesh.sizes["x"])
    Nz = int(mesh.sizes["z"])
    Lx = float(mesh.x.values[-1] - mesh.x.values[0])
    Lz = float(mesh.z.values[-1] - mesh.z.values[0])
    steps = np.array(mesh.step.values)
    times = np.array(mesh.time.values)

x = np.linspace(0, Lx/1000.0, Nx)
z = np.linspace(-Lz/1000.0, 0, Nz)
xx, zz  = np.meshgrid(x, z)
print(f'nx, nz, times: {Nx}, {Nz}, {len(times)}')

#lithology.nc may have a mesh different from the other files (e.g. finer): its own x and dz are used for the thickness
with xr.open_dataset(f"{model_path}/lithology.nc") as mesh_lithology:
    x_lithology = (mesh_lithology.x.values - mesh_lithology.x.values[0])/1000.0 #km
    dz_lithology = float(abs(mesh_lithology.z.values[1] - mesh_lithology.z.values[0]))/1000.0 #km
    dx_lithology = float(abs(mesh_lithology.x.values[1] - mesh_lithology.x.values[0]))/1000.0 #km

############################################################################################################################
# Plotting
start = 0
end = len(times)
step = 1

make_videos = True
# make_videos = False

make_gifs = True
# make_gifs = False

zip_files = True
# zip_files = False

#lithology.nc from tapioca has the Mandyoc layer codes (order of the layers in interfaces.txt, 0 = bottom layer).
#Each code gets the color of its layer (LITHOLOGY_COLORS of functions/mandyocIO.py).
lithology_colors = 'auto' #any number of layers: colors from the names of the layers (infos_*.txt), see get_lithology_colors
# lithology_colors = 'accordion' #models with a basaltic layer (rho = 2900 kg/m3) below the air (e.g. AR models)
# lithology_colors = 'rifting' #rifting models (e.g. Ho35S1)
# lithology_colors = 'double_keel'
# lithology_colors = 'svg_accordion' #SVG accordion scenario of tapioca (10 layers)

colors = get_lithology_colors(model_path, lithology_colors)

#Scenarios with surface processes: each deposited sediment layer gets a color of SEDIMENT_COLORS (functions/mandyocIO.py),
#in cycle, to show the stratigraphy (sediment_layers = True), or all of them the same color (False). No effect without sediments.
sediment_layers = True
# sediment_layers = False
roles = get_layer_roles(model_path)
sediment_codes = sorted(roles.get('sediments', []))
if(not sediment_layers):
    colors = [COLOR_SED if code in sediment_codes else color for code, color in enumerate(colors)]

#Codes of the oceanic crust, whose thickness is plotted in the top panel: the basalt layers (get_layer_roles), or a list of codes
oceanic_crust_codes = roles.get('basalt', [])
# oceanic_crust_codes = [4]
if(len(oceanic_crust_codes) == 0):
    print("Warning: no basalt layer found (get_layer_roles): the thickness of the oceanic crust is zero. Set oceanic_crust_codes.")
print(f"Oceanic crust codes: {oceanic_crust_codes}")

#The lithology from the particles is noisy column by column: the thickness is smoothed by a running mean of this width
#(0: no smoothing)
thickness_window = 1.0 #km

#Size of the figure [inches]: thickness panel, lithology panel and legend of the layers below them, as in frames_generator.py
#(numerical_setup.png of tapioca: title, scale factor and rheology and density of each layer)
fig_width, thickness_height, axes_height = 12, 1.5, 3
layers_legend = get_layers_legend(model_path, colors, thickness=False)
layers_legend_ncol, layers_legend_height = layers_legend_layout(layers_legend, fig_width)

#Use vmin=-0.5 and vmax=len(colors)-0.5 so the code k gets colors[k]
cmap = ListedColormap(colors)
np.seterr(divide='ignore')
fig_format = 'jpeg'

chunk_size = 4

isotherms = [500, 600, 700, 800, 900, 1300] #oC

xlims = [1600, 2200]
ylims = [-150, 40]
# xlims = [0, Lx/1000]
# ylims = [-Lz/1000+40, 0+40]
ylims_thickness = [0, 10] #km
if(args.xlims is not None):
    xlims = list(args.xlims)
if(args.ylims is not None):
    ylims = list(args.ylims)
print(f"xlims: {xlims} km, ylims: {ylims} km")

id_thickness = 0
id_model = 1

prefix = f"{model_name}_lithology_{'MeltThickness' if plot_melt else 'Thickness'}" #frames, video and gif

def plot_lithology_frame(axs, i):
    '''Lithology of lithology.nc shaded by the strain, isotherms and melt of the time step i, as in frames_generator.py'''
    air = density.isel(time=i).values < 1200.0 #tapioca keeps the values in the air: set them to zero there

    axs.imshow(lithology.isel(time=i).values[::-1, :], aspect='auto', extent=(0, Lx/1000, -Lz/1000+40, 40), cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, alpha=1.0)
    axs.imshow(np.log10(np.where(air, 0.0, strain.isel(time=i).values)[::-1, :]), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2, aspect='auto')
    if(plot_isotherms):
        axs.contour(x, z+40, temperature.isel(time=i).values, levels=isotherms, colors='r', linewidths=1.0)
    axs.set_ylim(ylims)
    axs.set_xlim(xlims)

    #strain legend: a bar with the shading, with fixed size [inches], whatever the size of the figure
    #(the layers are in the legend below the axes, see plot_layers_legend)
    plot_strain_legend(axs, vmin=-0.5, vmax=0.9, alpha=0.2, position=(0.98, 0.2), loc='lower right')

    if(plot_melt):
        incremental_melt_i = np.where(air, 0.0, incremental_melt.isel(time=i).values)
        melt_i = np.where(air, 0.0, melt.isel(time=i).values)
        color_incremental_melt = 'xkcd:bright pink'
        color_depleted_mantle = 'xkcd:purple'

        incremental_melt_i[incremental_melt_i == 0] = np.nan # Set zero values to NaN to avoid plotting them
        melt_i[melt_i == 0] = np.nan # Set zero values to NaN to avoid plotting them

        axs.contourf(xx, zz+40, incremental_melt_i, colors=color_incremental_melt, alpha=0.4, zorder=30)
        axs.contourf(xx, zz+40, melt_i, colors=color_depleted_mantle, alpha=0.4, zorder=20)
        #melt legend: inside the bottom left corner of the axis; the second label is placed to the right of the first one
        text_fsize = 11
        box = dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1.5)
        melt_label = axs.text(0.01, 0.04, r'Melt Fraction ($\partial \phi / \partial t$)', color='xkcd:bright pink',
                              fontsize=text_fsize, ha='left', va='bottom', bbox=box, transform=axs.transAxes, zorder=60)
        axs.annotate(r'Depleted Mantle ($\phi$)', xy=(1, 0.5), xycoords=melt_label, xytext=(8, 0),
                     textcoords='offset points', ha='left', va='center', color='xkcd:bright purple', fontsize=text_fsize,
                     bbox=box, zorder=60)

def plot_thickness_frame(axs, i):
    '''Thickness of the oceanic crust (oceanic_crust_codes) along x of the time step i'''
    thickness_oceanic_crust = np.isin(lithology.isel(time=i).values, oceanic_crust_codes).sum(axis=0) * dz_lithology #km
    n_window = int(round(thickness_window/dx_lithology))
    if(n_window > 1):
        thickness_oceanic_crust = np.convolve(thickness_oceanic_crust, np.ones(n_window)/n_window, mode='same')
    axs.plot(x_lithology, thickness_oceanic_crust, color='xkcd:dark grey blue')
    axs.set_ylim(ylims_thickness)
    axs.grid()

print("Generating frames...")

with pymp.Parallel() as p:
    for chunk_start in p.range(start, end, chunk_size):
        lithology = xr.open_dataset("lithology.nc", cache=False, engine="h5netcdf")["lithology"]
        strain = xr.open_dataset("strain.nc", cache=False, engine="h5netcdf")["strain"]
        temperature = xr.open_dataset("temperature.nc", cache=False, engine="h5netcdf")["temperature"]
        density = xr.open_dataset("density.nc", cache=False, engine="h5netcdf")["density"] #used to remove the values in the air
        if(plot_melt):
            incremental_melt = xr.open_dataset("dPhi.nc", cache=False, engine="h5netcdf")["dPhi"]
            melt = xr.open_dataset("Phi.nc", cache=False, engine="h5netcdf")["Phi"]

        for i in range(chunk_start, min(chunk_start + chunk_size, end)):
            steps_model = steps[i]
            current_time = float(times[i])

            fig, axs = plt.subplots(2, 1, figsize=(fig_width, thickness_height + axes_height + layers_legend_height),
                                    constrained_layout=True, sharex=True,
                                    gridspec_kw={'height_ratios': [thickness_height, axes_height]})

            #model name and time above the top panel
            axs[0].text(0.01, 1.07, f'{model_name}', fontsize=14, zorder=52, transform=axs[0].transAxes)
            axs[0].text(0.5, 1.07, f'Time = {current_time:.2f} Myr', fontsize=14, zorder=52, transform=axs[0].transAxes)

            plot_thickness_frame(axs[id_thickness], i)
            plot_lithology_frame(axs[id_model], i)
            plot_layers_legend(fig, layers_legend, layers_legend_ncol)

            # Setting plot details
            fsize = 14
            axs[-1].set_xlabel('Distance [km]', fontsize=fsize)
            axs[id_model].set_ylabel('Depth [km]', fontsize=fsize)
            axs[id_thickness].set_ylabel(r'Thick$_{\mathrm{oc}}$ [km]', fontsize=12)
            for ax in axs:
                ax.tick_params(axis='both', labelsize=fsize)

            figname = f"{prefix}_{str(int(steps_model)).zfill(6)}.{fig_format}"

            proc = psutil.Process(os.getpid())
            # print(f"before - worker ID:{os.getpid()}; Memory Usage: {proc.memory_info().rss/1024**3:.2f} GB")
            fig.savefig(f"_output/{figname}", dpi=300)
            plt.close(fig)

            del fig
            del axs
            gc.collect()

        del lithology
        del strain
        del temperature
        del density
        if(plot_melt):
            del incremental_melt
            del melt
        gc.collect()

print("Done!")

##############################################################################################################################################################################
videoname = f'{model_path}/_output/{prefix}'

if(make_videos):
    print("Generating videos...")

    fps = 24
    if os.path.isfile(f"{videoname}.mp4"):
        os.remove(f"{videoname}.mp4")
        print("\tRemoving previous video.")

    comand = f"ffmpeg -r {fps} -f image2 -s 1920x1080 -pattern_type glob -i \"{videoname}_[0-9]*.{fig_format}\" -vcodec libx264 -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -r 24 -y -an -crf 25 -pix_fmt yuv420p {videoname}.mp4"
    result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
    print("\tDone!")

##########################################################################################################################################################################
# Converting videos to gifs
# ss: skip seconds; t: duration time of the output; i: inputs format
# vf: filtergraph (fps: frames per second; scale: resize, e.g. 720 wide; lanczos: scaling algorithm;
#     palettegen and paletteuse: custom palette; split: everything done in one command)
# loop: 0: infinite; -1: no looping; n >= 0: n+1 loops

if(make_gifs):
    print("Converting videos to gifs...")
    if os.path.isfile(f"{videoname}.gif"):
        os.remove(f"{videoname}.gif")
        print("\tRemoving previous gif.")

    comand = f"ffmpeg -ss 0 -t 15 -i '{videoname}.mp4' -vf \"fps=30,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse\" -loop 0 {videoname}.gif"
    result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
    print("\tDone!")

##########################################################################################################################################################################

if(zip_files):
    #zip plots, videos and gifs
    print('Zipping figures, videos and gifs...')
    outputs_path = f'{model_path}/_output/'
    os.chdir(outputs_path)
    subprocess.run(f"zip {prefix}_imgs.zip {prefix}_[0-9]*.{fig_format}", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"zip {prefix}_videos.zip {prefix}.mp4", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"zip {prefix}_gifs.zip {prefix}.gif", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"rm {prefix}_[0-9]*.{fig_format}", shell=True, check=True, capture_output=True, text=True)
    print('Zipping complete!')
    os.chdir(f'{model_path}')
