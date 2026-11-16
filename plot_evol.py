#!/usr/bin/env python
# coding: utf-8

########################################################################
# BEFORE RUNNING THIS SCRIPT:                                          #
#  - CREATE THE .nc FILES WITH tapioca's meshes2netcdf.jl              #
#  - RUN THIS SCRIPT IN THE SAME FOLDER WHERE THE .nc FILES ARE        #
#                                                                      #
# Evolution of the scenario: the frames of frames_generator.py of      #
# n_snapshots instants (equally spaced steps from the first to the     #
# last one) stacked in a single figure, one per property, saved as     #
# _output/<model>_evol_<prop>[_MeltFrac].png                           #
########################################################################

import os
import sys
import argparse
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

matplotlib.use('agg')

path_to_functions = os.path.dirname(os.path.abspath(__file__)) # post_processing folder, with functions/mandyocIO.py
sys.path.append(os.path.abspath(path_to_functions))

if '' in sys.path:
    sys.path.remove('')
from functions.mandyocIO import plot_property, plot_lithology_legend, get_lithology_colors, read_tapioca_datasets, get_layer_roles, _legend_colors, COLOR_SED

####################################################################################################################################
#Options of the command line (the same as frames_generator.py, plus --n-snapshots): the settings below are used for
#the options that are not given
#e.g. python plot_evol.py --prop lithology temperature --plot-melt --xlims 600 1200 --ylims -150 40
PROPERTIES = ['lithology', 'density', 'radiogenic_heat', 'pressure', 'strain', 'strain_rate', 'temperature',
              'temperature_anomaly', 'surface', 'viscosity']
parser = argparse.ArgumentParser(description="Evolution of a Mandyoc model converted by tapioca: snapshots from the first to the last step in one figure (run in the folder of the scenario).")
parser.add_argument('--prop', nargs='+', choices=PROPERTIES, metavar='PROP',
                    help=f"properties to plot (default: properties of the script). Options: {', '.join(PROPERTIES)}")
parser.add_argument('--plot-melt', action=argparse.BooleanOptionalAction, default=None,
                    help="plot (--plot-melt) or not (--no-plot-melt) the melt (default: plot_melt of the script, 'auto': only if the model has melt)")
parser.add_argument('--xlims', nargs=2, type=float, metavar=('XMIN', 'XMAX'), help="horizontal limits [km] (default: the whole model)")
parser.add_argument('--ylims', nargs=2, type=float, metavar=('ZMIN', 'ZMAX'), help="vertical limits [km], 0 at the initial surface and 40 at the top of the air (default: the whole model)")
parser.add_argument('--n-snapshots', type=int, default=None, metavar='N', help="number of snapshots (default: n_snapshots of the script)")
args = parser.parse_args()

model_path = os.getcwd() # Get local file
model_name = model_path.split('/')[-1]
output_path = '_output'
print(f"Model name: {model_name}\n")
print(f"Model path: {model_path}\n")
print(f"Output path: {output_path}\n")

if not os.path.isdir(output_path):
    os.makedirs(output_path)

n_snapshots = 5 #equally spaced steps from the first to the last one
if(args.n_snapshots is not None):
    n_snapshots = args.n_snapshots

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

properties = [#Properties from mandyoc. Comment/uncomment to select which ones you would like to plot
            #  'density',
            #  'radiogenic_heat',
             'lithology', #default: lithology.nc shaded by the strain
            #  'pressure',
            #  'strain',
            #  'strain_rate',
            #  'temperature',
            #  'temperature_anomaly',
            #  'surface',
            #  'viscosity'
             ]
#'lithology' is plotted by this script (as in frames_generator.py); the other properties are plotted with
#plot_property (functions/mandyocIO.py), reading <property>.nc (heat.nc for radiogenic_heat). One figure per property.
if(args.prop is not None):
    properties = list(dict.fromkeys(args.prop))
print(f"Properties: {properties}")

#########################################
# Get domain and snapshots informations #
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

idxs = np.unique(np.round(np.linspace(0, len(times) - 1, n_snapshots)).astype(int)) #indices of the snapshots
print(f"Snapshots: steps {steps[idxs].tolist()}, times {np.round(times[idxs], 2).tolist()} Myr")

############################################################################################################################
# Plotting
plot_colorbar = True
fig_format = 'png'

#air/crust interface of the properties plotted with plot_property (not 'lithology'):
#False: from surface.nc (topography of Mandyoc); True: from the density (density.nc)
topo_from_density = False
# topo_from_density = True
if(not topo_from_density and not os.path.isfile(f"{model_path}/surface.nc")):
    print("No surface.nc: the air/crust interface is taken from the density (density.nc).")
    topo_from_density = True

#lithology.nc from tapioca has the Mandyoc layer codes (order of the layers in interfaces.txt, 0 = bottom layer).
#Each code gets the color of its layer (LITHOLOGY_COLORS of functions/mandyocIO.py).
lithology_colors = 'auto' #any number of layers: colors from the names of the layers (infos_*.txt), see get_lithology_colors
# lithology_colors = 'accordion'
# lithology_colors = 'rifting'
# lithology_colors = 'double_keel'
# lithology_colors = 'svg_accordion'

colors = get_lithology_colors(model_path, lithology_colors)

#Scenarios with surface processes: each deposited sediment layer gets a color of SEDIMENT_COLORS (True)
#or all of them the same color (False), as in frames_generator.py. No effect without sediments.
sediment_layers = True
# sediment_layers = False
sediment_codes = sorted(get_layer_roles(model_path).get('sediments', []))
if(not sediment_layers):
    colors = [COLOR_SED if code in sediment_codes else color for code, color in enumerate(colors)]
legend_colors = _legend_colors([COLOR_SED if code in sediment_codes else color for code, color in enumerate(colors)])

cmap = ListedColormap(colors)
np.seterr(divide='ignore')

isotherms = [500, 600, 700, 800, 900, 1300] #oC

xlims = [0, Lx/1000]
ylims = [-Lz/1000+40, 0+40]
ylims_surface = [-6, 6] #km, topography (prop = 'surface')
if(args.xlims is not None):
    xlims = list(args.xlims)
if(args.ylims is not None):
    ylims = list(args.ylims)
print(f"xlims: {xlims} km, ylims: {ylims} km")

def fig_prefix(prop):
    '''Prefix of the figure of prop: _output/<model>_evol_<prop>[_MeltFrac]'''
    prefix = f"{model_name}_evol_{prop}"
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
    '''Lithology of lithology.nc shaded by the strain, isotherms and melt of the time step i (as frames_generator.py)'''
    air = density.isel(time=i).values < 1200.0 #tapioca keeps the values in the air: set them to zero there

    axs.imshow(lithology.isel(time=i).values[::-1, :], aspect='auto', extent=(0, Lx/1000, -Lz/1000+40, 40), cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, alpha=1.0)
    axs.imshow(np.log10(np.where(air, 0.0, strain.isel(time=i).values)[::-1, :]), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2, aspect='auto')
    if(plot_isotherms):
        axs.contour(x, z+40, temperature.isel(time=i).values, levels=isotherms, colors='r', linewidths=1.0)
    axs.set_ylim(ylims)
    axs.set_xlim(xlims)

    #lithology/strain legend with fixed size [inches], whatever the size of the figure
    plot_lithology_legend(axs, colors=legend_colors, vmin=-0.5, vmax=0.9, alpha=0.2,
                          position=(0.97, 0.25), loc='lower right')

    if(plot_melt):
        incremental_melt_i = np.where(air, 0.0, incremental_melt.isel(time=i).values)
        melt_i = np.where(air, 0.0, melt.isel(time=i).values)
        color_incremental_melt = 'xkcd:bright pink'
        color_depleted_mantle='xkcd:purple'

        incremental_melt_i[incremental_melt_i == 0] = np.nan # Set zero values to NaN to avoid plotting them
        melt_i[melt_i == 0] = np.nan # Set zero values to NaN to avoid plotting them

        axs.contourf(xx, zz+40, incremental_melt_i, colors=color_incremental_melt, alpha=0.4, zorder=30)
        axs.contourf(xx, zz+40, melt_i, colors=color_depleted_mantle, alpha=0.4, zorder=20)
        text_fsize = 11
        box = dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1.5)
        melt_label = axs.text(0.01, 0.04, r'Melt Fraction ($\partial \phi / \partial t$)', color='xkcd:bright pink',
                              fontsize=text_fsize, ha='left', va='bottom', bbox=box, transform=axs.transAxes, zorder=60)
        axs.annotate(r'Depleted Mantle ($\phi$)', xy=(1, 0.5), xycoords=melt_label, xytext=(8, 0),
                     textcoords='offset points', ha='left', va='center', color='xkcd:bright purple', fontsize=text_fsize,
                     bbox=box, zorder=60)

if('lithology' in properties):
    lithology = xr.open_dataset("lithology.nc")["lithology"]
    strain = xr.open_dataset("strain.nc")["strain"]
    temperature = xr.open_dataset("temperature.nc")["temperature"]
    density = xr.open_dataset("density.nc")["density"] #used to remove the values in the air
    if(plot_melt):
        incremental_melt = xr.open_dataset("dPhi.nc")["dPhi"]
        melt = xr.open_dataset("Phi.nc")["Phi"]

if(len(other_properties) > 0):
    dataset = read_tapioca_datasets(model_path, variables)

print("Generating figures...")
fsize = 14
for prop in properties:
    #one panel per snapshot, each one with the size of a frame of frames_generator.py (12 x 3 in)
    fig, axes = plt.subplots(len(idxs), 1, figsize=(12, 3*len(idxs)), constrained_layout=True)
    axes = np.atleast_1d(axes)

    for axs, i in zip(axes, idxs):
        current_time = float(times[i])
        axs.text(0.01, 1.035, f'{model_name}', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs.transAxes)
        axs.text(0.5, 1.035, f'Time = {current_time:.2f} Myr', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs.transAxes)

        if(prop == 'lithology'):
            plot_lithology_frame(axs, i)
        else:
            plot_property(dataset.isel(time=i).load(), prop, xlims, ylims_surface if prop == 'surface' else ylims, model_path, fig, axs,
                          print_time=False,
                          plot_isotherms=plot_isotherms, isotherms=isotherms,
                          topo_from_density=topo_from_density,
                          plot_colorbar=plot_colorbar,
                          plot_melt=plot_melt,
                          lithology_colors=lithology_colors)

        axs.set_xlabel('Distance [km]', fontsize=fsize)
        axs.set_ylabel('Topography [km]' if prop == 'surface' else 'Depth [km]', fontsize=fsize)
        axs.tick_params(axis='both', labelsize=fsize)

    figname = f"{fig_prefix(prop)}.{fig_format}"
    fig.savefig(f"{output_path}/{figname}", dpi=300)
    plt.close(fig)
    print(f"\tSaved {output_path}/{figname}")

print("Done!")
