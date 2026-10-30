#!/usr/bin/env python
# coding: utf-8

########################################################################
# BEFORE RUNNING THIS SCRIPT:                                          #
#  - CREATE THE .nc FILES WITH tapioca's meshes2netcdf.jl              #
#    (temperature.nc, density.nc, strain.nc, lithology.nc, ...)        #
#  - RUN THIS SCRIPT IN THE SAME FOLDER WHERE THE .nc FILES ARE        #
#                                                                      #
# Snapshots of the scenario (based on Geometries.ipynb, without the    #
# measurements of the crust): one panel per instant, from the first    #
# to the last step, to see the general evolution of the scenario.     #
# The figure is saved as <model>_snapshots_<prop>.png                  #
########################################################################

import os
import sys
import string
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')

path_to_functions = os.path.dirname(os.path.abspath(__file__)) # post_processing folder, with functions/mandyocIO.py
sys.path.append(os.path.abspath(path_to_functions))

if '' in sys.path:
    sys.path.remove('')
from functions.mandyocIO import read_tapioca_datasets, read_params, read_times_steps, plot_property

####################################################################################################################################
model_path = os.getcwd() # Get local file
model_name = model_path.split('/')[-1]
print(f"Model name: {model_name}\n")
print(f"Model path: {model_path}\n")

##############################################################################################################
#Set mandyoc property

# prop = 'density'
# prop = 'radiogenic_heat'
prop = 'lithology'
# prop = 'pressure'
# prop = 'strain'
# prop = 'strain_rate'
# prop = 'temperature'
# prop = 'temperature_anomaly'
# prop = 'viscosity'

##############################################################################################################
#Colors of the Mandyoc layer codes in the lithology plot: 'auto' (model type from the number of layers in
#lithology.nc: 11 double keel, 7 rifting, 8 accordion), a key of LITHOLOGY_COLORS (functions/mandyocIO.py)
#or a list with one color per code (order of the layers in interfaces.txt, 0 = bottom layer, the air is the last one)

lithology_colors = 'auto'
# lithology_colors = 'double_keel'
# lithology_colors = 'rifting'
# lithology_colors = 'accordion'

##############################################################################################################
#Number of snapshots: equally spaced steps from the first to the last one

n_snapshots = 4

##############################################################################################################
#Reading scenario basic informations

times, steps = read_times_steps(model_path)
Nx, Nz, Lx, Lz = read_params(model_path)

idxs = np.unique(np.round(np.linspace(0, len(times) - 1, n_snapshots)).astype(int)) #indices of the snapshots
print(f"Snapshots: steps {steps[idxs].tolist()}, times {times[idxs].tolist()} Myr")

#Reading the .nc files used by plot_property (density: topography; temperature: isotherms; strain: lithology shading)
variables = ['density', 'temperature', 'strain']
if(prop not in ['lithology', 'temperature_anomaly']):
    variables.append(prop)
dataset = read_tapioca_datasets(model_path, variables)

##############################################################################################################
#Set domain limits of the figure

xlims = [0, Lx / 1.0e3]
ylims = [-Lz / 1.0e3 + 40, 40]

##############################################################################################################
#Creating figure and subplots
plt.close('all')

fig, axs = plt.subplots(len(idxs), 1, constrained_layout=True, figsize=(12, 3*len(idxs)), sharex=True, sharey=True)
axs = np.atleast_1d(axs)

letters = list(string.ascii_lowercase)[0: len(idxs)]

##############################################################################################################
#Plot data

for idx, ax in zip(idxs, axs):
    plot_property(dataset.isel(time=idx), prop, xlims, ylims, model_path, fig, ax,
                  print_time=False,
                  isotherms=[550, 800, 1300], label_isotherms=(ax == axs[0]),
                  plot_colorbar=(ax == axs[-1]),
                  bbox_to_anchor=(0.80, 0.09, 0.18, 0.16), #legend of the lithology: lower right corner at (0.98, 0.09)
                  lithology_colors=lithology_colors)

##############################################################################################################
#Setting plot details

fsize = 14
axs[-1].set_xlabel('x [km]', fontsize=fsize)
for ax, letter, idx in zip(axs, letters, idxs):
    #letter of the subplot and instant outside the plot, in the upper left corner
    label_text = letter + ') ' + f'{times[idx]:.2f}' + ' Myr'
    ax.text(0.0,
            1.02,
            label_text,
            color='k',
            fontsize=fsize,
            horizontalalignment='left',
            verticalalignment='bottom',
            transform=ax.transAxes)

    ax.set_ylabel('Depth [km]', fontsize=fsize)
    ax.tick_params(axis='both', labelsize=fsize)
    ax.set_yticks(np.linspace(ylims[0], 0, 5))

axs[0].set_xlim(xlims)
axs[0].set_ylim(ylims)

##############################################################################################################
#Saving figure inside the scenario directory

figname = f"{model_name}_snapshots_{prop}.png"
fig.savefig(f"{model_path}/{figname}", dpi=300)
print(f"Saved {figname}")
