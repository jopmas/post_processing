#!/usr/bin/env python
# coding: utf-8

import os
import gc
import sys
import multiprocessing #needed to run pymp in mac
import psutil
multiprocessing.set_start_method('fork') #needed to run pymp in mac 
# multiprocessing.set_start_method('spawn')
import pymp
import subprocess
import numpy as np
import xarray as xr
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


matplotlib.use('agg')

path = os.getcwd().split('/')
machine_path = f'/{path[1]}/{path[2]}' #cat the /home/user/ or /Users/user from system using path

path_to_functions = os.path.dirname(os.path.abspath(__file__)) # post_processing folder, with functions/mandyocIO.py
sys.path.append(os.path.abspath(path_to_functions))

if '' in sys.path:
    sys.path.remove('')
from functions.mandyocIO import change_dataset, plot_property, plot_lithology_legend, get_lithology_colors, read_tapioca_datasets

####################################################################################################################################
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
# plot_melt = True
plot_melt = False
plot_particles=False


if("4" in model_name and "0" in model_name):
    hcrust = 40.0e3 #m
else:
    hcrust = 35.0e3 #m

if(plot_isotherms or plot_melt):
    clean_plot=False
else:
    clean_plot = True

datasets = [#Properties from mandyoc. Comment/uncomment to select properties of the dataset
            # 'density',
            # 'radiogenic_heat',
            # 'pressure',
            'strain',
            # 'strain_rate',### Read ascii outputs and save them as xarray.Datasets,
            # 'surface',
            'temperature',
            # 'viscosity'
            ]# Read data and convert them to xarray.Dataset

properties = [#Properties from mandyoc. Comment/uncomment to select which ones you would like to plot
            #  'density',
            #  'radiogenic_heat',
            #  'lithology', #default: lithology.nc shaded by the strain
            #  'pressure',
            #  'strain',
            #  'strain_rate',
            #  'temperature',
            #  'temperature_anomaly',
            #  'surface',
             'viscosity'
             ]
#'lithology' is plotted by this script; the other properties are plotted with plot_property (functions/mandyocIO.py),
#reading <property>.nc (heat.nc for radiogenic_heat). One set of frames, video and gif per property.

#######################################################
# Read ascii outputs and save them as xarray.Datasets #
#######################################################

new_datasets = change_dataset(properties, datasets)
# print(new_datasets)
to_remove = []
remove_density=False
if ('density' not in properties): #used to plot air/curst interface
        properties.append('density')
        new_datasets = change_dataset(properties, datasets)
        to_remove.append('density')
        # remove_density=True

# if ('surface' not in properties): #used to plot air/curst interface
#         properties.append('surface')
#         new_datasets = change_dataset(properties, datasets)
#         to_remove.append('surface')
        # remove_density=True

if (plot_isotherms): #add datasets needed to plot isotherms
    if ('temperature' not in new_datasets):
        properties.append('temperature')
        new_datasets = change_dataset(properties, datasets)
        to_remove.append('temperature')

# print(f"newdatasets: {new_datasets}")

if (plot_melt): #add datasets needed to plot melt fraction
    if ('melt' not in properties):
        properties.append('melt')
    if ('incremental_melt' not in properties):
        properties.append('incremental_melt')
    new_datasets = change_dataset(properties, datasets)

    #removing the auxiliary datasets to not plot
    to_remove.append('melt')
    to_remove.append('incremental_melt')

if(clean_plot): #a clean plot
    new_datasets = change_dataset(properties, datasets)

for item in to_remove:
    properties.remove(item)
    
# dataset = read_datasets(model_path, new_datasets)
# print(dataset.keys())
# Normalize velocity values
# if ("velocity_x" and "velocity_z") in dataset.data_vars:
#     v_max = np.max((dataset.velocity_x**2 + dataset.velocity_z**2)**(0.5))    
#     dataset.velocity_x[:] = dataset.velocity_x[:] / v_max
#     dataset.velocity_z[:] = dataset.velocity_z[:] / v_max

# if ('lithology' in properties):
# lithology_dataset = xr.open_dataset(f"{model_path}/_lithology.nc")
# strain_dataset = xr.open_dataset(f"{model_path}/_strain.nc")
# temperature_dataset = xr.open_dataset(f"{model_path}/_temperature.nc")

# if(plot_melt):
#     Phi_dataset = xr.open_dataset(f"{model_path}/_melt.nc")
#     dPhi_dataset = xr.open_dataset(f"{model_path}/_incremental_melt.nc")

#########################################
# Get domain and particles informations #
#########################################

#tapioca files: data (time, z, x), x from 0 to Lx and z from -Lz to 0 [m]
with xr.open_dataset(f"{model_path}/strain.nc") as mesh:
    Nx = int(mesh.sizes["x"])
    Nz = int(mesh.sizes["z"])
    Lx = float(mesh.x.values[-1] - mesh.x.values[0])
    Lz = float(mesh.z.values[-1] - mesh.z.values[0])

x = np.linspace(0, Lx/1000.0, Nx)
z = np.linspace(-Lz/1000.0, 0, Nz)
xx, zz  = np.meshgrid(x, z)

# print(particles_layers)
############################################################################################################################
# Plotting
plot_colorbar = True
h_air = 40.0

start = 0
end = int(xr.open_dataset(f"{model_path}/strain.nc").time.size)
step = 1

make_videos = True
# make_videos = False

make_gifs = True
# make_gifs = False

zip_files = True
# zip_files = False

print("Generating frames...")

color_lower_crust='xkcd:brown'

color_incremental_melt = 'xkcd:bright pink'
color_depleted_mantle='xkcd:bright purple'
#air/crust interface of the properties plotted with plot_property (not 'lithology'):
#False: from surface.nc (topography of Mandyoc); True: from the density (density.nc)
topo_from_density = False
# topo_from_density = True
if(not topo_from_density and not os.path.isfile(f"{model_path}/surface.nc")):
    print("No surface.nc: the air/crust interface is taken from the density (density.nc).")
    topo_from_density = True

linewidth = 0.1
markersize = 4
line_alpha = 1.0
# color_crust='xkcd:grey'

# color_incremental_melt = 'xkcd:bright pink'
# color_depleted_mantle='xkcd:dark grey'

cr = 255.
color_air = (1.,1.,1.) # 5
color_bas = 'xkcd:dark grey blue' #(250./cr,50./cr,50./cr) # 4
color_uc = (228./cr,156./cr,124./cr) # 3
color_lc = (240./cr,209./cr,188./cr) # 2
color_lit = (155./cr,194./cr,155./cr) # 1
color_ast = (207./cr,226./cr,205./cr) # 0


#lithology.nc from tapioca has the Mandyoc layer codes (order of the layers in interfaces.txt, 0 = bottom layer).
#Each code gets the color of its layer (LITHOLOGY_COLORS of functions/mandyocIO.py).
lithology_colors = 'auto' #any number of layers: colors from the names of the layers (infos_*.txt), see get_lithology_colors
# lithology_colors = 'accordion' #models with a basaltic layer (rho = 2900 kg/m3) below the air (e.g. AR models)
# lithology_colors = 'rifting' #rifting models (e.g. Ho35S1)
# lithology_colors = 'double_keel'
# lithology_colors = 'svg_accordion' #SVG accordion scenario of tapioca (10 layers)

colors = get_lithology_colors(model_path, lithology_colors)

#Creating a custom colormap according to the list of colors defined above.
# This colormap will be used to plot the lithology mesh, where each lithology code is represented by a specific color.
# Use vmin=-0.5 and vmax=len(colors)-0.5 so the code k gets colors[k]

cmap = ListedColormap(colors)
np.seterr(divide='ignore')
# time = time[::-1]
steps = np.array(xr.open_dataset(f"{model_path}/strain.nc").step.values)
times = np.array(xr.open_dataset(f"{model_path}/strain.nc").time.values)
print(f'nx, nz, times: {Nx}, {Nz}, {len(times)}, type')
fig_format = 'jpeg'

# Files from tapioca keep the raw values in the air (density < 1200 kg/m3): they are set to zero to not be plotted
def worker(i):
    lithology = xr.open_dataset("lithology.nc", cache=False, engine="h5netcdf")["lithology"].isel(time=i).values[::-1,:]
    density = xr.open_dataset("density.nc", cache=False, engine="h5netcdf")["density"].isel(time=i).values
    strain  = np.where(density < 1200.0, 0.0, xr.open_dataset("strain.nc", cache=False, engine="h5netcdf")["strain"].isel(time=i).values)[::-1,:]
    temperature = xr.open_dataset("temperature.nc", cache=False, engine="h5netcdf")["temperature"].isel(time=i).values
    return lithology, strain, temperature

A = np.zeros((100, 10))

A[:25, :] = 2700
A[25:50, :] = 2800
A[50:75, :] = 3300
A[75:100, :] = 3400

A = A[::-1, :]

xA = np.linspace(-0.5, 0.9, 10)
yA = np.linspace(0, 1.5, 100)

xxA, yyA = np.meshgrid(xA, yA)

chunk_size = 4

isotherms = [500, 600, 700, 800, 900, 1300] #oC

xlims = [0, Lx/1000]
ylims = [-Lz/1000+40, 0+40]
# xlims = [1600, 2200]
# ylims = [-150, 40]
ylims_surface = [-6, 6] #km, topography (prop = 'surface')

def frame_prefix(prop):
    '''Prefix of the frames, video and gif of prop: _output/<model>_<prop>[_MeltFrac]'''
    prefix = f"{model_name}_{prop}"
    if(plot_melt):
        prefix = f"{prefix}_MeltFrac"
    return prefix

#properties plotted with plot_property and the .nc files that it needs
#(density: air and topography; temperature: isotherms; Phi.nc and dPhi.nc: melt)
other_properties = [prop for prop in properties if prop != 'lithology']
variables = ['density', 'temperature']
variables += [prop for prop in other_properties if prop not in ['temperature_anomaly', 'surface']]
if(plot_melt):
    variables += ['melt', 'incremental_melt']

def plot_lithology_frame(axs, i):
    '''Lithology of lithology.nc shaded by the strain, isotherms and melt of the time step i (prop = 'lithology')'''
    air = density.isel(time=i).values < 1200.0 #tapioca keeps the values in the air: set them to zero there

    axs.imshow(lithology.isel(time=i).values[::-1, :], aspect='auto', extent=(0, Lx/1000, -Lz/1000+40, 40), cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, alpha=1.0)

    # axs.imshow(np.log10(strain), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2)
    axs.imshow(np.log10(np.where(air, 0.0, strain.isel(time=i).values)[::-1, :]), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2, aspect='auto')
    if(plot_isotherms):
        axs.contour(x, z+40, temperature.isel(time=i).values, levels=isotherms, colors='r', linewidths=1.0)
    axs.set_ylim(ylims)
    axs.set_xlim(xlims)

    #lithology/strain legend with fixed size [inches], whatever the size of the figure
    plot_lithology_legend(axs, colors=[color_uc, color_lc, color_lit, color_ast], vmin=-0.5, vmax=0.9, alpha=0.2,
                          position=(0.97, 0.25), loc='lower right')

    if(plot_melt):
        incremental_melt_i = np.where(air, 0.0, incremental_melt.isel(time=i).values)[::,:]
        melt_i = np.where(air, 0.0, melt.isel(time=i).values)[::,:]
        #Plotting incremental melt
        color_incremental_melt = 'xkcd:bright pink'
        color_depleted_mantle='xkcd:purple'

        incremental_melt_i[incremental_melt_i == 0] = np.nan # Set zero values to NaN to avoid plotting them
        melt_i[melt_i == 0] = np.nan # Set zero values to NaN to avoid plotting them

        axs.contourf(xx, zz+40, incremental_melt_i, colors=color_incremental_melt, alpha=0.4, zorder=30)
        axs.contourf(xx, zz+40, melt_i, colors=color_depleted_mantle, alpha=0.4, zorder=20)
        #plotting melt legend
        text_fsize = 12
        axs.text(0.10, 1.035, r'Melt Fraction $\left(\frac{\partial \phi}{\partial t}\right)$', color='xkcd:bright pink', fontsize=text_fsize, transform=axs.transAxes, zorder=60)
        axs.text(0.25, 1.035, r'Depleted Mantle ($\phi$)', color='xkcd:bright purple', fontsize=text_fsize, transform=axs.transAxes, zorder=60)

with pymp.Parallel() as p:
    # for i in p.range(start, end-step, step):
    for chunk_start in p.range(start, end, chunk_size):
        if('lithology' in properties):
            lithology = xr.open_dataset("lithology.nc", cache=False, engine="h5netcdf")["lithology"]
            strain = xr.open_dataset("strain.nc", cache=False, engine="h5netcdf")["strain"]
            temperature = xr.open_dataset("temperature.nc", cache=False, engine="h5netcdf")["temperature"]
            density = xr.open_dataset("density.nc", cache=False, engine="h5netcdf")["density"] #used to remove the values in the air
            if(plot_melt):
                incremental_melt = xr.open_dataset("dPhi.nc", cache=False, engine="h5netcdf")["dPhi"]
                melt = xr.open_dataset("Phi.nc", cache=False, engine="h5netcdf")["Phi"]

        if(len(other_properties) > 0):
            dataset = read_tapioca_datasets(model_path, variables, cache=False, engine="h5netcdf")

        for i in range(chunk_start, min(chunk_start + chunk_size, end)):
            steps_model = steps[i] #dataset.steps.values[i]
            current_time = float(times[i])

            if(len(other_properties) > 0):
                dataset_i = dataset.isel(time=i).load()

            for prop in properties:
                fig, axs = plt.subplots(1, 1, figsize=(12, 3), constrained_layout=True)

                axs.text(0.01, 1.035, f'{model_name}', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs.transAxes)
                axs.text(0.5, 1.035, f'Time = {current_time:.2f} Myr', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs.transAxes)

                if(prop == 'lithology'):
                    plot_lithology_frame(axs, i)
                else:
                    plot_property(dataset_i, prop, xlims, ylims_surface if prop == 'surface' else ylims, model_path, fig, axs,
                                  print_time=False,
                                  plot_isotherms=plot_isotherms, isotherms=isotherms,
                                  topo_from_density=topo_from_density,
                                  plot_colorbar=plot_colorbar,
                                  plot_melt=plot_melt,
                                  lithology_colors=lithology_colors)

                # Setting plot details
                fsize = 14
                axs.set_xlabel('Distance [km]', fontsize=fsize)
                axs.set_ylabel('Topography [km]' if prop == 'surface' else 'Depth [km]', fontsize=fsize)
                axs.tick_params(axis='both', labelsize=fsize)

                figname = f"{frame_prefix(prop)}_{str(int(steps_model)).zfill(6)}.{fig_format}"

                proc = psutil.Process(os.getpid())
                # print(f"before - worker ID:{os.getpid()}; Memory Usage: {proc.memory_info().rss/1024**3:.2f} GB")
                fig.savefig(f"_output/{figname}", dpi=300)
                # print(f"after - worker ID:{os.getpid()}; Memory Usage: {proc.memory_info().rss/1024**3:.2f} GB")
                plt.close(fig)

                del fig
                del axs
                gc.collect()

        if('lithology' in properties):
            del lithology
            del strain
            del temperature
            del density
            if(plot_melt):
                del incremental_melt
                del melt
        if(len(other_properties) > 0):
            del dataset
        gc.collect()

        # objs = gc.get_objects()
        # print(f'figures remaining: {sum(isinstance(o, Figure) for o in objs)}')

print("Done!")

##############################################################################################################################################################################
if(make_videos):
    print("Generating videos...")

    fps = 24
    for prop in properties:
        videoname = f'{model_path}/_output/{frame_prefix(prop)}'

        if(plot_particles):
            if(prop == 'viscosity'):
                videoname = f'{videoname}'
            else:
                videoname = f'{videoname}_particles'
                # videoname = f'{videoname}_particles_onlymb'
            
        try:
            comand = f"rm {videoname}.mp4"
            result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
            print(f"\tRemoving previous {prop} video.")
        except:
            print(f"\tNo {prop} video to remove.")

        comand = f"ffmpeg -r {fps} -f image2 -s 1920x1080 -pattern_type glob -i \"{videoname}_[0-9]*.{fig_format}\" -vcodec libx264 -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -r 24 -y -an -crf 25 -pix_fmt yuv420p {videoname}.mp4"
        result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
    print("\tDone!")


##########################################################################################################################################################################

# # Converting videos to gifs
# 
# ss: skip seconds
# 
# t: duration time of the output
# 
# i: inputs format
# 
# vf: filtergraph (video filters)
# 
#     - fps: frames per second
# 
#     - scale: resize accordint to given pixels (e.g. 1080 = 1080p wide)
#     
#     - lanczos: scaling algorithm
#     
#     - palettegen and palette use: filters that generate a custom palette
#     
#     - split: filter that allows everything to be done in one command
# 
# loop: number of loops
# 
#     - 0: infinite
# 
#     - -1: no looping
# 
#     - for numbers n >= 0, create n+1 loops


if(make_gifs):
    print("Converting videos to gifs...")
    for prop in properties:
        gifname = f'{model_path}/_output/{frame_prefix(prop)}'

        if(plot_particles):
            if(prop == 'viscosity'):
                gifname = f'{gifname}'
            else:
                gifname = f'{gifname}_particles'
                # gifname = f'{gifname}_particles_onlymb'
            

        try:
            comand = f"rm {gifname}.gif"
            result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
            print(f"\tRemoving previous {prop} gif.")
        except:
            print(f"\tNo {prop} gif to remove.")
        
        comand = f"ffmpeg -ss 0 -t 15 -i '{gifname}.mp4' -vf \"fps=30,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse\" -loop 0 {gifname}.gif"
        result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True) 
    print("\tDone!")

##########################################################################################################################################################################

if(zip_files):
    #zip plots, videos and gifs
    print('Zipping figures, videos and gifs...')
    outputs_path = f'{model_path}/_output/'
    os.chdir(outputs_path)
    subprocess.run(f"zip {model_name}_imgs.zip *.{fig_format}", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"zip {model_name}_videos.zip *.mp4", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"zip {model_name}_gifs.zip *.gif", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"rm *.{fig_format}", shell=True, check=True, capture_output=True, text=True)
    print('Zipping complete!')
    os.chdir(f'{model_path}')