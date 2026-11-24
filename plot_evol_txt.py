#!/usr/bin/env python
# coding: utf-8

########################################################################
# BEFORE RUNNING THIS SCRIPT:                                          #
#  - RUN THIS SCRIPT IN THE FOLDER OF THE SCENARIO, WITH THE .txt      #
#    FILES OF MANDYOC (the .nc files of tapioca are not needed)        #
#                                                                      #
# Evolution of the scenario from the .txt files of Mandyoc, to follow  #
# a scenario while it runs: the same figure of plot_evol.py            #
# (n_snapshots instants, equally spaced steps from the first to the    #
# last one, stacked in a single figure, one per property), saved as    #
# _output/<model>_evol_txt_<prop>[_MeltFrac].png                       #
#                                                                      #
# The .txt files can be in the folder of the scenario (Mandyoc running)#
# or in one folder per variable (after 0_organize_outputs.sh):         #
# density_<step>.txt or density/density_<step>.txt, and the lithology  #
# in litho_<step>_<rank>.txt or lithos/litho_<step>_<rank>.txt         #
########################################################################

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

matplotlib.use('agg')

path_to_functions = os.path.dirname(os.path.abspath(__file__)) # post_processing folder, with functions/mandyocIO.py
sys.path.append(os.path.abspath(path_to_functions))

if '' in sys.path:
    sys.path.remove('')
from functions.mandyocIO import plot_property, plot_strain_legend, get_lithology_colors, get_layer_roles, COLOR_SED, get_layers_legend, layers_legend_layout, plot_layers_legend

####################################################################################################################################
#Options of the command line (the same as plot_evol.py): the settings below are used for the options that are not given
#e.g. python plot_evol_txt.py --prop lithology temperature --plot-melt --xlims 600 1200 --ylims -150 40
#'surface' is not an option: plot_property reads it from surface.nc
PROPERTIES = ['lithology', 'density', 'radiogenic_heat', 'pressure', 'strain', 'strain_rate', 'temperature',
              'temperature_anomaly', 'viscosity']
parser = argparse.ArgumentParser(description="Evolution of a Mandyoc model from its .txt files (while it runs): snapshots from the first to the last step in one figure (run in the folder of the scenario).")
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
plot_melt = 'auto' #only if the model has melt (Phi and dPhi .txt files, with melt in the last step), see model_has_melt
# plot_melt = True
# plot_melt = False

properties = [#Properties from mandyoc. Comment/uncomment to select which ones you would like to plot
            #  'density',
            #  'radiogenic_heat',
             'lithology', #default: litho_<step>_<rank>.txt shaded by the strain
            #  'pressure',
            #  'strain',
            #  'strain_rate',
            #  'temperature',
            #  'temperature_anomaly',
            #  'viscosity'
             ]
#'lithology' is plotted by this script (as in plot_evol.py); the other properties are plotted with plot_property
#(functions/mandyocIO.py), reading <property>_<step>.txt (heat_<step>.txt for radiogenic_heat). One figure per property.
if(args.prop is not None):
    properties = list(dict.fromkeys(args.prop))
print(f"Properties: {properties}")

########################################
# Readers of the .txt files of Mandyoc #
########################################

#name of the .txt files of each property (as the .nc files of tapioca, read_tapioca_datasets of functions/mandyocIO.py)
file_names = {'melt': 'Phi', 'incremental_melt': 'dPhi', 'radiogenic_heat': 'heat'}

def read_param(model_path):
    '''Nx, Nz, Lx, Lz [m] of param.txt (as read_param of tapioca's meshes2netcdf.jl)'''
    param = {}
    with open(f"{model_path}/param.txt") as file:
        for line in file:
            line = line.split('#')[0].replace(' ', '').strip() #without the comments and the spaces
            if('=' in line):
                key, value = line.split('=', 1)
                param[key.lower()] = value
    return int(param['nx']), int(param['nz']), float(param['lx']), float(param['lz'])

def txt_path(var, step, rank=None):
    '''Path of the .txt file of var in step (and rank, for the lithology): in the folder of the scenario
    (Mandyoc running) or in the folder of the variable (after 0_organize_outputs.sh, 'lithos' for the lithology)'''
    fname = f"{var}_{step}.txt" if rank is None else f"{var}_{step}_{rank}.txt"
    folder = 'lithos' if var == 'litho' else var
    if(os.path.isfile(f"{model_path}/{folder}/{fname}")):
        return f"{model_path}/{folder}/{fname}"
    return f"{model_path}/{fname}"

def get_steps(model_path):
    '''Steps with time_<step>.txt, sorted'''
    files = glob.glob(f"{model_path}/time_*.txt") + glob.glob(f"{model_path}/time/time_*.txt")
    steps = [os.path.basename(f)[len('time_'):-len('.txt')] for f in files]
    return np.array(sorted({int(step) for step in steps if step.lstrip('-').isdigit()}))

def read_time(step):
    '''Time of step [Myr]: first line of time_<step>.txt, e.g. "[0]  0:   1.3870e+07" (years)'''
    with open(txt_path('time', step)) as file:
        line = file.readline()
    return float(line.split()[-1])/1.0e6

def read_mesh_txt(var, step):
    '''Field var of step in a (Nz, Nx) array, first row at the bottom (as the .nc files of tapioca).
    The file has 2 lines of header ("Vec Object..." and "type: mpi") and a line "Process [n]" before the values of
    each MPI process, with x varying faster (as read_data of tapioca's meshes2netcdf.jl)'''
    values = pd.read_csv(txt_path(var, step), header=None, skiprows=2, comment='P', dtype=np.float64).to_numpy(copy=True)[:, 0] #copy: pandas >= 3 gives read-only arrays
    values[np.abs(values) < 1.0e-50] = 0.0 #as tapioca
    return values.reshape(Nz, Nx)

def number_of_ranks(model_path):
    '''Number of MPI processes of the lithology files (litho_<step>_<rank>.txt) of the first step'''
    return len(glob.glob(txt_path('litho', steps[0], '*')))

def fill_negatives_with_neighbors(litho):
    '''Cells without particles (-1) get the code of a neighbor (up, down, left, right), as
    replace_negatives_with_neighbors! of tapioca's meshes2netcdf.jl'''
    for shift, axis in [(1, 0), (-1, 0), (1, 1), (-1, 1)]:
        neighbor = np.roll(litho, shift, axis=axis)
        #np.roll wraps around: the borders do not get the code of the opposite side
        if(axis == 0):
            neighbor[0 if shift == 1 else -1, :] = -1
        else:
            neighbor[:, 0 if shift == 1 else -1] = -1
        litho = np.where((litho < 0) & (neighbor >= 0), neighbor, litho)
    return litho

def read_litho_txt(step):
    '''Lithology of step (Mandyoc layer codes, order of the layers in interfaces.txt, 0 = bottom layer) in a (Nzl, Nxl)
    array, first row at the bottom (as lithology.nc of tapioca), on a mesh 5 times finer than the other fields.
    Each litho_<step>_<rank>.txt has the columns ix, iz (index of the cell of the finer mesh, iz = 0 at the top) and
    the code of the layer (as convert_litho_to_nc of tapioca's meshes2netcdf.jl)'''
    litho = np.full((Nzl, Nxl), -1, dtype=np.int16) #-1: cell without particles
    for rank in range(nranks):
        fpath = txt_path('litho', step, rank)
        if(not os.path.isfile(fpath)):
            continue
        data = pd.read_csv(fpath, header=None, sep=' ', comment='P', usecols=[0, 1, 2], dtype=np.int32).values
        litho[data[:, 1], data[:, 0]] = data[:, 2]
    litho = fill_negatives_with_neighbors(litho)
    return litho[::-1, :] #first row at the bottom

def step_is_complete(step, variables):
    '''True if step has the .txt files of all variables (and of the lithology of all MPI processes). The last step of
    a scenario that is running can still be being written'''
    for var in variables:
        if(not os.path.isfile(txt_path(file_names.get(var, var), step))):
            return False
    if('lithology' in properties):
        return len(glob.glob(txt_path('litho', step, '*'))) == nranks
    return True

#########################################
# Get domain and snapshots informations #
#########################################

#tapioca convention: data (time, z, x), x from 0 to Lx and z from -Lz to 0 [m]
Nx, Nz, Lx, Lz = read_param(model_path)
Nxl, Nzl = 5*(Nx - 1) + 1, 5*(Nz - 1) + 1 #mesh of the lithology

x = np.linspace(0, Lx/1000.0, Nx)
z = np.linspace(-Lz/1000.0, 0, Nz)
xx, zz  = np.meshgrid(x, z)

steps = get_steps(model_path)
if(len(steps) == 0):
    sys.exit(f"No time_<step>.txt in {model_path} or {model_path}/time")
nranks = number_of_ranks(model_path) if 'lithology' in properties else 0
if('lithology' in properties and nranks == 0):
    sys.exit(f"No litho_{steps[0]}_<rank>.txt in {model_path} or {model_path}/lithos")

def model_has_melt(model_path):
    '''True if the model has the Phi and dPhi .txt files and melt (Phi > 0 out of the air) in the last step'''
    last_step = steps[-1]
    if not (os.path.isfile(txt_path('Phi', last_step)) and os.path.isfile(txt_path('dPhi', last_step))):
        return False
    melt_last = read_mesh_txt('Phi', last_step)
    air = read_mesh_txt('density', last_step) < 1200.0 #Mandyoc keeps the values in the air
    return bool(np.nanmax(np.where(air, 0.0, melt_last)) > 0.0)

if(args.plot_melt is not None):
    plot_melt = args.plot_melt
if(plot_melt == 'auto'):
    plot_melt = model_has_melt(model_path)
    print(f"Melt: {'found, plotting it' if plot_melt else 'not found (no Phi/dPhi .txt files or no melt), not plotting it'}")

#properties plotted with plot_property and the .txt files that it needs
#(density: air and topography; temperature: isotherms; Phi and dPhi: melt)
other_properties = [prop for prop in properties if prop != 'lithology']
variables = ['density', 'temperature']
if('lithology' in properties):
    variables += ['strain'] #shading of the lithology
variables += [prop for prop in other_properties if prop not in ['temperature_anomaly']]
if(plot_melt):
    variables += ['melt', 'incremental_melt']
variables = list(dict.fromkeys(variables))

#only the steps with all the files (the last step can still be being written)
steps = np.array([step for step in steps if step_is_complete(step, variables)])
if(len(steps) == 0):
    sys.exit("No step with all the .txt files of the properties")

idxs = np.unique(np.round(np.linspace(0, len(steps) - 1, n_snapshots)).astype(int)) #indices of the snapshots
steps = steps[idxs] #from here on, only the steps of the snapshots
times = np.array([read_time(step) for step in steps])
idxs = np.arange(len(steps))
print(f"Snapshots: steps {steps.tolist()}, times {np.round(times, 2).tolist()} Myr")

#############################################################################
# Reading the .txt files of the snapshots in a Dataset as the .nc of tapioca #
#############################################################################
#the same structure of read_tapioca_datasets (functions/mandyocIO.py): data (time, z, x) with the coordinates x and z
#[m], time [Myr] and step, so that plot_property and plot_lithology_frame work as with the .nc files
print("Reading the .txt files...")
coords = {'time': times, 'x': np.linspace(0, Lx, Nx), 'z': np.linspace(-Lz, 0, Nz)}
data_vars = {}
for var in variables:
    print(f"\t{file_names.get(var, var)}")
    fields = np.array([read_mesh_txt(file_names.get(var, var), step) for step in steps])
    name = {'melt': 'Phi', 'incremental_melt': 'dPhi'}.get(var, var) #names of the variables in Phi.nc and dPhi.nc
    data_vars[name] = (('time', 'z', 'x'), fields)
dataset = xr.Dataset(data_vars, coords=coords)
dataset = dataset.assign_coords(step=('time', steps))

density = dataset['density'] #used to remove the values in the air
temperature = dataset['temperature']
if('lithology' in properties):
    strain = dataset['strain']
    print("\tlitho")
    lithology = xr.DataArray(np.array([read_litho_txt(step) for step in steps]), dims=('time', 'z', 'x'),
                             coords={'time': times, 'x': np.linspace(0, Lx, Nxl), 'z': np.linspace(-Lz, 0, Nzl)})
    if(plot_melt):
        incremental_melt = dataset['dPhi']
        melt = dataset['Phi']

############################################################################################################################
# Plotting
plot_colorbar = True
fig_format = 'png'

#air/crust interface of the properties plotted with plot_property (not 'lithology'): from the density,
#since plot_property reads the topography of Mandyoc from surface.nc
topo_from_density = True

#litho_<step>_<rank>.txt has the Mandyoc layer codes (order of the layers in interfaces.txt, 0 = bottom layer).
#Each code gets the color of its layer (LITHOLOGY_COLORS of functions/mandyocIO.py).
lithology_colors = 'auto' #any number of layers: colors from the names of the layers (infos_*.txt), see get_lithology_colors
# lithology_colors = 'accordion'
# lithology_colors = 'rifting'
# lithology_colors = 'double_keel'
# lithology_colors = 'svg_accordion'

#'auto' needs the layer codes (number of layers and the layer around the weak seeds): those of the first snapshot,
#since there is no lithology.nc (without the lithology, from infos_*.txt or interfaces.txt)
lithology_first = lithology.isel(time=0) if 'lithology' in properties else None
colors = get_lithology_colors(model_path, lithology_colors, lithology=lithology_first)

#Scenarios with surface processes: each deposited sediment layer gets a color of SEDIMENT_COLORS (True)
#or all of them the same color (False), as in frames_generator.py. No effect without sediments.
sediment_layers = True
# sediment_layers = False
sediment_codes = sorted(get_layer_roles(model_path).get('sediments', []))
if(not sediment_layers):
    colors = [COLOR_SED if code in sediment_codes else color for code, color in enumerate(colors)]

#legend of the layers below the figure of the lithology, as in plot_evol.py (title, scale factor and rheology
#and density of each layer; thickness=True adds the initial thickness); the figure gets the height of the legend
#besides the height of the panels
fig_width, axes_height = 12, 3 #inches, size of each panel (a frame of frames_generator.py)
layers_legend = get_layers_legend(model_path, colors, thickness=False, lithology=lithology_first) if 'lithology' in properties else []
layers_legend_ncol, layers_legend_height = layers_legend_layout(layers_legend, fig_width)

cmap = ListedColormap(colors)
np.seterr(divide='ignore')

isotherms = [550, 700, 800, 1300] #oC

xlims = [0, Lx/1000]
ylims = [-Lz/1000+40, 0+40]
if(args.xlims is not None):
    xlims = list(args.xlims)
if(args.ylims is not None):
    ylims = list(args.ylims)
print(f"xlims: {xlims} km, ylims: {ylims} km")

def fig_prefix(prop):
    '''Prefix of the figure of prop: _output/<model>_evol_txt_<prop>[_MeltFrac]'''
    prefix = f"{model_name}_evol_txt_{prop}"
    if(plot_melt):
        prefix = f"{prefix}_MeltFrac"
    return prefix

def plot_lithology_frame(axs, i):
    '''Lithology of litho_<step>_<rank>.txt shaded by the strain, isotherms and melt of the time step i (as plot_evol.py)'''
    air = density.isel(time=i).values < 1200.0 #Mandyoc keeps the values in the air: set them to zero there

    axs.imshow(lithology.isel(time=i).values[::-1, :], aspect='auto', extent=(0, Lx/1000, -Lz/1000+40, 40), cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, alpha=1.0)
    axs.imshow(np.log10(np.where(air, 0.0, strain.isel(time=i).values)[::-1, :]), extent=(0, Lx/1000, -Lz/1000+40, 40), cmap="Greys", vmin=-0.5, vmax=0.9, alpha=0.2, aspect='auto')
    if(plot_isotherms):
        axs.contour(x, z+40, temperature.isel(time=i).values, levels=isotherms, colors='r', linewidths=1.0)
    axs.set_ylim(ylims)
    axs.set_xlim(xlims)

    #strain legend: a bar with the shading, with fixed size [inches], whatever the size of the figure
    #(the layers are in the legend below the figure, see plot_layers_legend)
    plot_strain_legend(axs, vmin=-0.5, vmax=0.9, alpha=0.2, position=(0.98, 0.2), loc='lower right')

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

print("Generating figures...")
fsize = 14
for prop in properties:
    #one panel per snapshot, each one with the size of a frame of frames_generator.py (12 x 3 in)
    height = axes_height*len(idxs) + (layers_legend_height if prop == 'lithology' else 0)
    fig, axes = plt.subplots(len(idxs), 1, figsize=(fig_width, height), constrained_layout=True)
    axes = np.atleast_1d(axes)

    for axs, i in zip(axes, idxs):
        current_time = float(times[i])
        axs.text(0.01, 1.035, f'{model_name}', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs.transAxes)
        axs.text(0.5, 1.035, f'Time = {current_time:.2f} Myr', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=axs.transAxes)
        axs.text(0.99, 1.035, f'Step {steps[i]}', bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, ha='right', transform=axs.transAxes)

        if(prop == 'lithology'):
            plot_lithology_frame(axs, i)
        else:
            plot_property(dataset.isel(time=i), prop, xlims, ylims, model_path, fig, axs,
                          print_time=False,
                          plot_isotherms=plot_isotherms, isotherms=isotherms,
                          topo_from_density=topo_from_density,
                          plot_colorbar=plot_colorbar,
                          plot_melt=plot_melt,
                          lithology_colors=lithology_colors)

        axs.set_xlabel('Distance [km]', fontsize=fsize)
        axs.set_ylabel('Depth [km]', fontsize=fsize)
        axs.tick_params(axis='both', labelsize=fsize)

    if(prop == 'lithology'):
        plot_layers_legend(fig, layers_legend, layers_legend_ncol)

    figname = f"{fig_prefix(prop)}.{fig_format}"
    fig.savefig(f"{output_path}/{figname}", dpi=300)
    plt.close(fig)
    print(f"\tSaved {output_path}/{figname}")

print("Done!")
