#!/usr/bin/env python
# coding: utf-8

import os
import gc
import sys
import multiprocessing #needed to run pymp in mac
multiprocessing.set_start_method('fork') #needed to run pymp in mac
import pymp
import subprocess
import numpy as np
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patches as patches

from matplotlib.colors import ListedColormap

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
from functions.mandyocIO import find_nearest, plot_lithology_legend

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

if(plot_isotherms or plot_melt):
    clean_plot=False
else:
    clean_plot = True

datasets = [#Properties from mandyoc. Comment/uncomment to select properties of the dataset
            'density',
            'radiogenic_heat',
            'pressure',
            'strain',
            'strain_rate',### Read ascii outputs and save them as xarray.Datasets,
            # 'surface',
            'temperature',
            'viscosity'
            ]# Read data and convert them to xarray.Dataset

properties = [#Properties from mandyoc. Comment/uncomment to select which ones you would like to plot
            #  'density',
            #  'radiogenic_heat',
             'lithology',
            #  'pressure',
            #  'strain',
            #  'strain_rate',
            #  'temperature',
            #  'temperature_anomaly',
            #  'surface',
            #  'viscosity'
             ]

#########################################
# Get domain and particles informations #
#########################################
#tapioca files: data (time, z, x), x from 0 to Lx and z from -Lz to 0 [m]
with xr.open_dataset(f"{model_path}/strain.nc") as mesh:
    Nx = int(mesh.sizes["x"])
    Nz = int(mesh.sizes["z"])
    Lx = float(mesh.x.values[-1] - mesh.x.values[0])
    Lz = float(mesh.z.values[-1] - mesh.z.values[0])
    steps_mesh = np.array(mesh.step.values)
    times_mesh = np.array(mesh.time.values)

x = np.linspace(0, Lx/1000.0, Nx)
z = np.linspace(-Lz/1000.0, 0, Nz)
xx, zz  = np.meshgrid(x, z)

trackdataset = xr.open_dataset("_track_xzPT_all_steps.nc")
x_track = trackdataset.xtrack.values[::-1]
z_track = trackdataset.ztrack.values[::-1]
P = trackdataset.ptrack.values[::-1]
T = trackdataset.ttrack.values[::-1]
time = trackdataset.time.values[::-1]
steps = trackdataset.step.values[::-1]

n = int(trackdataset.ntracked.values)
nTotal = np.size(x_track)
steps = nTotal//n

print(f"len of:\n x_track: {len(x_track)}\n z_track: {len(z_track)}\n P: {len(P)}\n T: {len(T)}\n time: {len(time)}\n n_tracked: {n}\n steps: {steps}\n")
print(f"n_tracked x len(all_time) = {n}*{len(time)} = {n*len(time)}")
print(f"nTotal: {nTotal}, n: {n}, steps: {steps}")

x_track = np.reshape(x_track,(steps,n))
z_track = np.reshape(z_track,(steps,n))
P = np.reshape(P,(steps,n))/1.0e3 #GPa
T = np.reshape(T,(steps,n))
T_maxs = np.nanmax(T, axis=0) #get the maximum temperature for each particle to categorize by temperature
P_maxs = np.nanmax(P, axis=0) #get the maximum pressure for each particle to categorize by temperature

particles_layers = trackdataset.particles_layers.values[::-1] #code of the tracked layers

double_keel = True
# double_keel = False
print(f"double_keel model") if(double_keel) else print("homogeneous model")

if(double_keel):
    asthenosphere_code = 0
    lower_craton_code = 1
    upper_craton_code = 2
    mantle_lithosphere_code = 3
    lower_crust1_code = 4
    seed_code = 5
    lower_crust2_code = 6
    upper_crust_code = 7
    decolement_code = 8
    sediments_code = 9
    air_code = 10

else:
    asthenosphere_code = 0
    mantle_lithosphere_code = 1
    lower_crust1_code = 2
    seed_code = 3
    lower_crust2_code = 4
    upper_crust_code = 5
    decolement_code = 6
    sediments_code = 7
    air_code = 8

T_initial = T[0]
P_initial = P[0]

if(double_keel==True):
    if(lower_craton_code in particles_layers):
        plot_lower_craton_particles = True
    else:
        plot_lower_craton_particles = False
        cond_lower_craton2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if(upper_craton_code in particles_layers):
        plot_upper_craton_particles = True
    else:
        plot_upper_craton_particles = False
        cond_upper_craton2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

if(mantle_lithosphere_code in particles_layers):
    plot_mantle_lithosphere_particles = True
else:
    plot_mantle_lithosphere_particles = False
    cond_mantle_lithosphere2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1
    
if((lower_crust1_code in particles_layers) | (lower_crust2_code in particles_layers)):
    plot_lower_crust_particles = True
else:
    plot_lower_crust_particles = False
    cond_lower_crust_2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

if(upper_crust_code in particles_layers):
    plot_upper_crust_particles = True
else:
    plot_upper_crust_particles = False
    cond_upper_crust2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

if(decolement_code in particles_layers):
    plot_decolement_particles = True
else:
    plot_decolement_particles = False
    cond_decolement2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

if(sediments_code in particles_layers):
    plot_sediments_particles = True
else:
    plot_sediments_particles = False
    cond_sediments2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

############################################################################################################################
# Plotting
plot_colorbar = True
h_air = 40.0

start = 0
end = int(times_mesh.size)
step = 4

# start = 4
# end = 5
# step = 1

make_videos = True
# make_videos = False

make_gifs = True
# make_gifs = False

zip_files = True
# zip_files = False

print("Generating frames...")
linewidth = 0.1
markersize = 4
line_alpha = 1.0
# color_crust='xkcd:grey'

color_incremental_melt = 'xkcd:bright pink'
color_depleted_mantle='xkcd:bright purple'

def plot_particles_of_a_layer_temperature_coded(axs, i, current_time, time, x_track, z_track, P, T, T_max, P_max, particle, markersize=4, linewidth=0.15, line_alpha=0.6, h_air=40.0e3,
                                                color_wedge='xkcd:orange', #orogenic wedge rocks T<900 and P<2GPa
                                               color_mid_T='xkcd:teal', #subducted crustal rocks 400<T<900 and P>2GPa
                                               color_high_T_wedge='xkcd:raspberry', #orogenic wedge rocks T>900 and P<2GPa
                                               color_high_T_subducted='xkcd:dirty purple', #subducted crustal rocks T>900 and P>2GPa
                                                plot_steps=False):
    '''
    Plot the particles of a specific layer with temperature coding.

    Parameters:
    - axs: The axes to plot on.
    - i: The current time index.
    - current_time: The current time.
    - time: The time array.
    - x_track: The x coordinates of the particles.
    - z_track: The z coordinates of the particles.
    - P: The pressure array.
    - T: The temperature array.
    - T_max: The maximum temperature of the particle.
    - particle: The particle index.
    - markersize: The size of the markers.
    - linewidth: The width of the lines.
    - line_alpha: The alpha transparency of the lines.
    - h_air: The height of the air.
    - color_low_T: The color for low temperature particles.
    - color_mid_T: The color for intermediate temperature particles.
    - color_high_T: The color for high temperature particles.
    - plot_steps: Whether to plot the steps.
    '''

    if not (T_max >= 400.0): # particles that never reach 400 oC (or without temperature) are not plotted
        return

    # if(T_max >= 400.0 and T_max < 900.0):
    #     color = color_mid_T
    # if(T_max >= 900.0):
    #     if(P_max < 2.0):
    #         color = color_high_T_wedge
    #     if(P_max >= 2.0):
    #         color = color_high_T_subducted

    if(T_max >= 400 and T_max < 900.0):
        if(P_max < 2.0):
            color = color_wedge
        if(P_max >= 2.0):
            color = color_mid_T
    if(T_max >= 900.0):
        if(P_max < 2.0):
            color = color_high_T_wedge
        if(P_max >= 2.0):
            color = color_high_T_subducted
    
    # if(T_max >=0 and T_max < 400.0):
    #     color = color_low_T
    # if(T_max >= 400.0 and T_max < 900.0):
    #     color = color_mid_T
    # if(T_max >= 900.0):
    #     color = color_high_T

    axs[0].plot(x_track[i, particle]/1.0e3,
                    z_track[i, particle]/1.0e3+h_air/1.0e3,
                    marker='.',
                    markersize=1.5,#markersize-2,
                    color=color,
                    zorder=60)
    axs[1].plot(T[i, particle], P[i, particle], '.', color=color, markersize=markersize)
    axs[1].plot(T[:i, particle], P[:i, particle], '-', color=color, linewidth=linewidth, alpha=line_alpha, zorder=60)

    if(plot_steps):
        for j in np.arange(0, current_time, 20):
            idx = find_nearest(time, j)
            axs[1].plot(T[idx, particle], P[idx, particle], '.', color='black', markersize=0.5)

def plot_PT_fields(ax, vertices, x_text, y_text, label, color, fsize_text=6, linewidth=0.01, alpha=0.3):
    """
    Plot the P-T fields as a polygon on the given axis.
    
    Parameters:
    - ax: The axis to plot on.
    - vertices: The vertices of the polygon.
    - color: The color of the polygon.
    - label: The label for the polygon.
    """
    polygon = patches.Polygon(
        vertices,
        closed=True,
        facecolor=color,
        alpha=alpha,
        edgecolor=color,
        linewidth=linewidth
        )
    ax.add_patch(polygon)

    ax.text(
            x_text,
            y_text,
            label,   # texto
            fontsize=fsize_text,
            color=color,
            ha='center',           # centraliza horizontalmente
            va='center',           # centraliza verticalmente
        )

#################################################
# Polygon coordenates of ultra high temperature #
#################################################
vertices_UHT = np.array([
    [900, 0/1.0e3],
    [1100, 0/1.0e3],
    [1100, 1400/1.0e3],
    [900, 1150/1.0e3]
])
color_UHT = 'xkcd:red'
label_UHT = 'UHT\nfield'
x_text_UHT = np.mean(vertices_UHT[:, 0])
y_text_UHT = np.mean(vertices_UHT[:, 1])
alpha_UHT = 0.1
fsize_UHT = 6

###################################################
# Polygon coordenates of high pressure granulites #
###################################################
vertices_granulite_HP = np.array([
    [674, 860/1.0e3],
    [1100, 1400/1.0e3],
    [1100, 2400/1.0e3],
    [1000, 2300/1.0e3],
    [615,1360/1.0e3]
])
color_granulite_HP = 'xkcd:brown'
label_granulite_HP = 'HP\ngranulites'
x_text_granulite_HP = np.mean(vertices_granulite_HP[:, 0])
y_text_granulite_HP = np.mean(vertices_granulite_HP[:, 1])
alpha_granulite_HP = 0.3
fsize_granulite_HP = 6

###################################
# Polygon coordenates of eclogite #
###################################
vertices_eclogite = np.array([
    [575, 1265/1.0e3],
    [1000, 2300/1.0e3],
    [1100, 2400/1.0e3],
    [1100, 4500/1.0e3],
    [600,4500/1.0e3],
    [347,2300/1.0e3],
    [468,2159/1.0e3],
    [488,2083/1.0e3],
    [531,1360/1.0e3]
])
color_eclogite = 'xkcd:blue'
label_eclogite = 'Eclogite\nfacies'
x_text_eclogite = np.mean(vertices_eclogite[:, 0])
y_text_eclogite = np.mean(vertices_eclogite[:, 1])
alpha_eclogite = 0.2
fsize_eclogite = 6

#####################################
# Polygon coordenates of granulites #
#####################################
vertices_granulite = np.array([
    [800, 0/1.0e3],
    [900, 0/1.0e3],
    [900, 1180/1.0e3],
    [670, 860/1.0e3],
    [684, 589/1.0e3],
    [710, 318/1.0e3],
    [748, 118/1.0e3]
])
color_granulite = 'xkcd:pink'
label_granulite = 'Granulites'
x_text_granulite = np.mean(vertices_granulite[:, 0]) + 30
y_text_granulite = np.mean(vertices_granulite[:, 1])
alpha_granulite = 0.5
fsize_granulite = 6

# Polygon coordenates of granulites facies anfibolite
vertices_anfibolite = np.array([
    [420, 0/1.0e3],
    [790, 0/1.0e3],
    [748, 118/1.0e3],
    [710, 318/1.0e3],
    [674,861/1.0e3],
    [615, 1360/1.0e3],
    [576, 1268/1.0e3],
    [531, 1360/1.0e3],
    [410, 859/1.0e3]
])
color_anfibolite = 'xkcd:green'
label_anfibolite = 'Anfibolite\nfacies'
x_text_anfibolite = np.mean(vertices_anfibolite[:, 0]) - 40
y_text_anfibolite = np.mean(vertices_anfibolite[:, 1])
alpha_anfibolite = 0.5
fsize_anfibolite = 6

# Polygon coordenates of blueschist
vertices_blueschist = np.array([
    [100, 344/1.0e3],
    [408, 859/1.0e3],
    [531, 1360/1.0e3],
    [488, 2080/1.0e3],
    [468, 2160/1.0e3],
    [348, 2318/1.0e3],
    [100, 668/1.0e3]
])
color_blueschist = 'xkcd:blue'
label_blueschist = 'Blueschist\nfacies'
x_text_blueschist = np.mean(vertices_blueschist[:, 0])
y_text_blueschist = np.mean(vertices_blueschist[:, 1])
alpha_blueschist = 0.4
fsize_blueschist = 6

# Polygon coordenates of greenschiest
vertices_greenschiest = np.array([
    [331, 0/1.0e3],
    [420, 0/1.0e3],
    [408, 859/1.0e3],
    [349, 778/1.0e3]
])
color_greenschiest = 'xkcd:green'
label_greenschiest = 'Greenschist\nfacies'
x_text_greenschiest = np.mean(vertices_greenschiest[:, 0])
y_text_greenschiest = np.mean(vertices_greenschiest[:, 1])
alpha_greenschiest = 0.3
fsize_greenschiest = 4

# Polygon coordenates of sub-greenschiest
vertices_sub_greenschiest = np.array([
    [100, 0/1.0e3],
    [331, 0/1.0e3],
    [349, 778/1.0e3],
    [100, 344/1.0e3]
])
color_sub_greenschiest = 'xkcd:green'
label_sub_greenschiest = 'Sub-Greenschist\nfacies'
x_text_sub_greenschiest = np.mean(vertices_sub_greenschiest[:, 0]) - 0#20
y_text_sub_greenschiest = np.mean(vertices_sub_greenschiest[:, 1])
alpha_sub_greenschiest = 0.1
fsize_sub_greenschiest = 4

color_wedge = 'xkcd:orange' #orogenic wedge rocks 400<T<900 and P<2GPa
color_mid_T = 'xkcd:teal' #subducted crustal rocks 400<T<900 and P>2GPa
color_high_T_wedge = 'xkcd:raspberry' #orogenic wedge rocks T>=900 and P<2GPa
color_high_T_subducted = 'xkcd:dirty purple' #subducted crustal rocks T>900 and P>2GPa


steps = steps_mesh
times = times_mesh

cr = 255.
color_air = (1.,1.,1.)
color_sed = (241./cr,184./cr,68./cr)
color_dec = (137./cr,81./cr,151./cr)
color_uc = (228./cr,156./cr,124./cr)
color_lc = (240./cr,209./cr,188./cr)
color_lit = (155./cr,194./cr,155./cr)
color_mlit_uc = (180. / cr, 194. / cr, 162. / cr)
color_mlit_lc = (155. / cr, 194. / cr, 155. / cr)
color_ast = (207./cr,226./cr,205./cr)


#lithology.nc from tapioca has the Mandyoc layer codes (order of the layers in interfaces.txt, 0 = bottom layer).
#Each code gets the color of its layer.
if(double_keel):
    lithology_colors = {asthenosphere_code: color_ast,
                        lower_craton_code: color_mlit_lc,
                        upper_craton_code: color_mlit_uc,
                        mantle_lithosphere_code: color_lit,
                        lower_crust1_code: color_lc,
                        seed_code: color_lc,
                        lower_crust2_code: color_lc,
                        upper_crust_code: color_uc,
                        decolement_code: color_dec,
                        sediments_code: color_sed,
                        air_code: color_air}
else:
    lithology_colors = {asthenosphere_code: color_ast,
                        mantle_lithosphere_code: color_lit,
                        lower_crust1_code: color_lc,
                        seed_code: color_lc,
                        lower_crust2_code: color_lc,
                        upper_crust_code: color_uc,
                        decolement_code: color_dec,
                        sediments_code: color_sed,
                        air_code: color_air}

colors = [lithology_colors[code] for code in range(air_code+1)]

#Creating a custom colormap according to the list of colors defined above.
# This colormap will be used to plot the lithology mesh, where each lithology code is represented by a specific color.
# Use vmin=-0.5 and vmax=len(colors)-0.5 so the code k gets colors[k]

cmap = ListedColormap(colors)
np.seterr(divide='ignore')

#strain shading limits (log10) used in the plot and in the legend
strain_vmin = -1.5
strain_vmax = 0.0
strain_alpha = 0.2

with pymp.Parallel() as p:
    lithology = xr.open_dataset(f"{model_path}/lithology.nc", cache=False, engine="h5netcdf")["lithology"]
    strain = xr.open_dataset(f"{model_path}/strain.nc", cache=False, engine="h5netcdf")["strain"]
    temperature = xr.open_dataset(f"{model_path}/temperature.nc", cache=False, engine="h5netcdf")["temperature"]
    density = xr.open_dataset(f"{model_path}/density.nc", cache=False, engine="h5netcdf")["density"] #used to remove the strain in the air

    for i in p.range(start, end, step):
        
        # data = dataset.isel(time=i)
        current_time = float(times[i])

        # print(f"Plotting frame {i} of {end} at time {current_time} Myr...")
        # for prop in properties:
        fig, axs = plt.subplots(1, 2, figsize=(12, 3), constrained_layout=True, gridspec_kw={'width_ratios': [1, 0.4]})
        axs[0].text(0.5, 1.035, f'Time = {current_time:.2f} Myr', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs[0].transAxes)

        xlims = [0, Lx / 1.0e3]
        ylims = [-Lz / 1.0e3 + 40, 40]
        # ylims = [-150, 40]
        # ylims = [-400, 40]

        axs[0].imshow(lithology.isel(time=i).values[::-1, :], aspect='auto', extent=(0, Lx/1000, -Lz/1000+40, 40), cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, alpha=1.0)

        # strain  = strain_dataset.strain.values[i][::-1,:]
        # temperature = temperature_dataset.temperature.values[i][::-1,:]

        #tapioca keeps the strain in the air: set it to zero there (density < 1200 kg/m3) to not shade the air
        strain_i = np.where(density.isel(time=i).values < 1200.0, 0.0, strain.isel(time=i).values)

        # axs.imshow(np.log10(strain), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2)
        axs[0].imshow(np.log10(strain_i[::-1, :]), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=strain_vmin, vmax=strain_vmax, alpha=strain_alpha, aspect='auto')
        # axs[0].imshow(strain.isel(time=i).values[::-1, :], extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.7, alpha=0.2, aspect='auto')
        axs[0].contour(x, z+40, temperature.isel(time=i).values, levels=[500, 600, 700, 800, 900, 1300], colors='r', linewidths=1.0)

        if(plot_colorbar):
            #lithology/strain legend with fixed size [inches], whatever the size of the figure
            plot_lithology_legend(axs[0], colors=[color_uc, color_lc, color_lit, color_ast],
                                  vmin=strain_vmin, vmax=strain_vmax, alpha=strain_alpha,
                                  position=(0.05, 0.20), loc='lower left')
              

        for particle, particle_layer, T_max, P_max in zip(range(n), particles_layers, T_maxs, P_maxs):
            #Plot particles in prop subplot


            if((plot_lower_crust_particles == True) & ((particle_layer == lower_crust1_code) | (particle_layer == lower_crust2_code))):
                plot_particles_of_a_layer_temperature_coded(axs, i, current_time, time, x_track, z_track, P, T, T_max, P_max, particle,
                                                        markersize=markersize, linewidth=linewidth, line_alpha=line_alpha,
                                                        color_wedge=color_wedge, #orogenic wedge rocks 400<T<900 and P<2GPa
                                                        color_mid_T=color_mid_T, #subducted crustal rocks 400<T<900 and P>2GPa
                                                        color_high_T_wedge=color_high_T_wedge, #orogenic wedge rocks T>=900 and P<2GPa
                                                        color_high_T_subducted=color_high_T_subducted, #subducted crustal rocks T>900 and P>2GPa
                                                        plot_steps=False)
        
        # Setting plot details

        # Stability fields
        plot_PT_fields(axs[1], vertices_UHT, x_text=x_text_UHT, y_text=y_text_UHT, label=label_UHT,
                        color=color_UHT, fsize_text=fsize_UHT, alpha=alpha_UHT)
        plot_PT_fields(axs[1], vertices_granulite_HP, x_text=x_text_granulite_HP, y_text=y_text_granulite_HP, label=label_granulite_HP,
                        color=color_granulite_HP, fsize_text=fsize_granulite_HP, alpha=alpha_granulite_HP)
        plot_PT_fields(axs[1], vertices_eclogite, x_text=x_text_eclogite, y_text=y_text_eclogite, label=label_eclogite,
                        color=color_eclogite, fsize_text=fsize_eclogite, alpha=alpha_eclogite)
        plot_PT_fields(axs[1], vertices_granulite, x_text=x_text_granulite, y_text=y_text_granulite, label=label_granulite,
                        color=color_granulite, fsize_text=fsize_granulite, alpha=alpha_granulite)
        plot_PT_fields(axs[1], vertices_anfibolite, x_text=x_text_anfibolite, y_text=y_text_anfibolite, label=label_anfibolite,
                        color=color_anfibolite, fsize_text=fsize_anfibolite, alpha=alpha_anfibolite)
        plot_PT_fields(axs[1], vertices_blueschist, x_text=x_text_blueschist, y_text=y_text_blueschist, label=label_blueschist,
                        color=color_blueschist, fsize_text=fsize_blueschist, alpha=alpha_blueschist)
        plot_PT_fields(axs[1], vertices_greenschiest, x_text=x_text_greenschiest, y_text=y_text_greenschiest, label=label_greenschiest,
                        color=color_greenschiest, fsize_text=fsize_greenschiest, alpha=alpha_greenschiest)
        plot_PT_fields(axs[1], vertices_sub_greenschiest, x_text=x_text_sub_greenschiest, y_text=y_text_sub_greenschiest, label=label_sub_greenschiest,
                        color=color_sub_greenschiest, fsize_text=fsize_sub_greenschiest, alpha=alpha_sub_greenschiest)


        fsize = 14
        axs[0].set_xlabel('Distance [km]', fontsize=fsize)
        axs[0].set_ylabel('Depth [km]', fontsize=fsize)
        axs[0].tick_params(axis='both', labelsize=fsize)

        axs[1].set_xlim([100, 1100])
        ylims = np.array([0, 4500])/1.0e3 #GPa
        axs[1].set_ylim(ylims)
        axs[1].set_xlabel(r'Temperature [$^{\circ}$C]', fontsize=fsize)
        axs[1].set_ylabel('Pressure [GPa]', fontsize=fsize)
        # axs[1].yaxis.set_label_position("right")
        # axs[1].tick_params(axis='y', labelright=True, labelleft=False, labelsize=fsize)
        
        axs[1].tick_params(axis='both', labelsize=fsize-2)
        axs[1].grid('-k', alpha=0.7)

        #creating depth axis to PTt plot
        ax1 = axs[1].twinx()
        ax1.set_ylim(ylims*1000/30)
        # ax1.tick_params(axis='y', labelright=False, labelleft=True, labelsize=fsize)
        ax1.set_ylabel('Depth [km]', fontsize=fsize)
        ax1.tick_params(axis='y', labelsize=fsize-2)
        # nticks = len(axs[1].get_yticks())
        # ax1.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(nticks))
        # ax1.yaxis.set_label_position("left")

        axs[1].set_xlabel(r'Temperature [$^{\circ}$C]', fontsize=fsize)
        axs[1].plot([-10,-10], [-10,-10], '-', color=color_wedge, markersize=markersize, label='Wedge LT', zorder=60)
        axs[1].plot([-10,-10], [-10,-10], '-', color=color_mid_T, markersize=markersize, label='Subducted LT', zorder=60)
        axs[1].plot([-10,-10], [-10,-10], '-', color=color_high_T_wedge, markersize=markersize, label='Wedge HT', zorder=60)
        axs[1].plot([-10,-10], [-10,-10], '-', color=color_high_T_subducted, markersize=markersize, label='Subducted HT', zorder=60)

        axs[1].legend(loc='upper left', ncol=1, fontsize=8, handlelength=0, handletextpad=0, labelcolor='linecolor')

        # ax1.legend(loc='center left', ncol=1, fontsize=6)

        if(plot_melt):
            #plotting melt legend
            text_fsize = 12
            axs[0].text(0.01, 0.90, r'Melt Fraction $\left(\frac{\partial \phi}{\partial t}\right)$', color='xkcd:bright pink', fontsize=text_fsize, transform=axs[0].transAxes, zorder=60)
            axs[0].text(0.01, 0.80, r'Depleted Mantle ($\phi$)', color='xkcd:bright purple', fontsize=text_fsize, transform=axs[0].transAxes, zorder=60)

            figname = f"{model_name}_lithology_and_PTt_temperature_coded_MeltFrac_{str(int(steps[i])).zfill(6)}.png"
        else:
            figname = f"{model_name}_lithology_and_PTt_temperature_coded_{str(int(steps[i])).zfill(6)}.png"
        fig.savefig(f"_output/{figname}", dpi=300)
        plt.close('all')

        # del data
        gc.collect()

print("Done!")

##############################################################################################################################################################################
if(make_videos):
    print("Generating videos...")

    fps = 24
    for prop in properties:
        videoname = f'{model_path}/_output/{model_name}_lithology_and_PTt'

        if(plot_melt):
            videoname = f'{videoname}_MeltFrac'

        if(plot_particles):
            if(prop == 'viscosity'):
                videoname = f'{videoname}'
            else:
                videoname = f'{videoname}_particles'
                # videoname = f'{videoname}_particles_onlymb'
            
        try:
            comand = f"rm {videoname}.mp4"
            result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
            print(f"\tRemoving previous lithology video.")
        except:
            print(f"\tNo lithology video to remove.")

        comand = f"ffmpeg -r {fps} -f image2 -s 1920x1080 -pattern_type glob -i \"{videoname}_*.png\" -vcodec libx264 -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -r 24 -y -an -crf 25 -pix_fmt yuv420p {videoname}.mp4"
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
        gifname = f'{model_path}/_output/{model_name}_lithology_and_PTt'

        if(plot_melt):
            gifname = f'{gifname}_MeltFrac'

        if(plot_particles):
            if(prop == 'viscosity'):
                gifname = f'{gifname}'
            else:
                gifname = f'{gifname}_particles'
                # gifname = f'{gifname}_particles_onlymb'
            

        try:
            comand = f"rm {gifname}.gif"
            result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True)
            print(f"\tRemoving previous lithology gif.")
        except:
            print(f"\tNo lithology gif to remove.")
        
        comand = f"ffmpeg -ss 0 -t 15 -i '{gifname}.mp4' -vf \"fps=30,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse\" -loop 0 {gifname}.gif"
        result = subprocess.run(comand, shell=True, check=True, capture_output=True, text=True) 
    print("\tDone!")

##########################################################################################################################################################################

if(zip_files):
    #zip plots, videos and gifs
    print('Zipping figures, videos and gifs...')
    outputs_path = f'{model_path}/_output/'
    os.chdir(outputs_path)
    subprocess.run(f"zip {model_name}_imgs.zip *.png", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"zip {model_name}_videos.zip *.mp4", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"zip {model_name}_gifs.zip *.gif", shell=True, check=True, capture_output=True, text=True)
    subprocess.run(f"rm *.png", shell=True, check=True, capture_output=True, text=True)
    print('Zipping complete!')