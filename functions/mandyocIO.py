"""
Functions to plot the outputs of Mandyoc converted to NetCDF by tapioca (meshes2netcdf.jl and particles2netcdf.jl),
used by the scripts and notebooks of post_processing.
"""
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter


def plot_lithology_legend(ax, colors, vmin=-0.5, vmax=0.7, alpha=0.2, position=(0.98, 0.08), loc='lower right',
                          width=0.74, height=0.66, fsize=10):
    '''
    Plot the lithology/strain legend box as an inset of ax. The size of the box is fixed in inches, so it is
    the same whatever the size of the figure (the default is the size of the box in the frames of
    plot_PTT_frames_dk.py). Only the position is relative to ax.

    Parameters
    ----------
    ax: `~matplotlib.axes.Axes`
        Axis where the lithology is plotted.

    colors: list
        Lithology colors from top to bottom (e.g. [upper crust, lower crust, lithospheric mantle, asthenosphere]).

    vmin, vmax: float
        Limits of the log10(strain) shading. Must be the same used to plot the strain.

    alpha: float
        Alpha transparency of the strain shading. Must be the same used to plot the strain.

    position: tuple
        (x, y) of the corner of the box given by loc, in fraction of ax.

    loc: str
        Corner of the box placed in position (e.g. 'lower right', 'lower left', 'upper right').

    width, height: float
        Size of the box [inches].

    fsize: int
        Font size of the label and ticks.

    Returns
    -------
    bv1: `~matplotlib.axes.Axes`
        Axis of the legend box.
    '''

    bv1 = inset_axes(ax,
                     width=width, #inches
                     height=height, #inches
                     loc=loc,
                     bbox_to_anchor=position,
                     bbox_transform=ax.transAxes,
                     borderpad=0
                     )

    #lithology bands (first color on top)
    nbands = len(colors)
    bands = np.repeat(np.arange(nbands), 25)[:, None] * np.ones((1, 10))
    bv1.imshow(bands, extent=[vmin, vmax, 0, 1.5], cmap=ListedColormap(colors), vmin=-0.5, vmax=nbands-0.5, aspect='auto')

    #strain shading
    xA = np.linspace(vmin, vmax, 10)
    bv1.imshow(np.tile(xA, (100, 1)), extent=[vmin, vmax, 0, 1.5], zorder=100, alpha=alpha, cmap='Greys', vmin=vmin, vmax=vmax, aspect='auto')

    bv1.set_xlabel(r"log$(\varepsilon_{II})$", size=fsize)
    bv1.tick_params(axis='x', which='major', labelsize=fsize)
    bv1.set_xticks(np.linspace(vmin, vmax, 3))
    bv1.set_yticks([])
    bv1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    return bv1

def add_colorbar(fig, ax, im, label, position=(0.02, 0.08), loc='lower left', length=1.8, thickness=0.11,
                 orientation='horizontal', fsize=10, fmt=None, ticks=None):
    '''
    Add a colorbar as an inset of ax. The size of the colorbar is fixed in inches, so it is the same whatever the
    size of the figure. Only the position is relative to ax.

    Parameters
    ----------
    fig: `~matplotlib.figure.Figure`
        Figure of ax.

    ax: `~matplotlib.axes.Axes`
        Axis where the data is plotted.

    im: mappable
        Output of imshow, contourf, etc.

    label: str
        Label of the colorbar.

    position: tuple
        (x, y) of the corner of the colorbar given by loc, in fraction of ax.

    loc: str
        Corner of the colorbar placed in position (e.g. 'lower left', 'lower right', 'upper right').

    length, thickness: float
        Size of the colorbar [inches].

    orientation: str
        'horizontal' or 'vertical'.

    fsize: int
        Font size of the label and ticks.

    fmt: str or callable
        Format of the tick labels (e.g. '%.0f').

    ticks: list
        Position of the ticks.

    Returns
    -------
    clb: `~matplotlib.colorbar.Colorbar`
    '''

    if(orientation == 'horizontal'):
        width, height = length, thickness
    else:
        width, height = thickness, length

    cax = inset_axes(ax,
                     width=width, #inches
                     height=height, #inches
                     loc=loc,
                     bbox_to_anchor=position,
                     bbox_transform=ax.transAxes,
                     borderpad=0
                     )

    clb = fig.colorbar(im, cax=cax, orientation=orientation, format=fmt, ticks=ticks)
    clb.set_label(label, fontsize=fsize)
    clb.ax.tick_params(labelsize=fsize)
    clb.minorticks_off()

    return clb

def find_nearest(array, value):
    '''Return the index in array nearest to a given value.
    
    Parameters
    ----------
    
    array: array_like
        1D array used to find the index
        
    value: float
        Value to be seached
    '''
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def _extract_interface(z, Z, Nx, field_datai, value_to_search):
    '''
    Extract interface from Rhoi according to a given density (rho)

    Parameters
    ----------
    z: array_like
        Array representing z direction.

    Z: array_like
        Array representing z direction resampled with higher resolution.

    Nx: int
        Number of points in x direction.

    field_data: array_like (Nz, Nx)
        Mandyoc field from mandyoc (e.g.: density, temperature, etc).

    value_to_search: float
        Value of density to be searched in field_data array.

    Returns
    -------
    mapped_interface: array_like
        Interface extracted from Rhoi field
    '''

    mapped_interface = []

    for j in np.arange(Nx):

        interp_function = interp1d(z, field_datai[:,j]) #return a "function" of interpolation to apply in other array
        
        interface_inverted = interp_function(Z)[::-1] #inverted array: from top to bottom

        idx = np.argmax(interface_inverted > value_to_search) #first occurrence of rho in the inverted array
        idx_corrected = len(interface_inverted) - idx #correcting the index to the original array

        depth = Z[idx_corrected]

        mapped_interface = np.append(mapped_interface, depth)

    return mapped_interface


def _log_fmt(x, pos):
    return "{:.0f}".format(np.log10(x))


def change_dataset(properties, datasets):
    '''
    Create new_datasets based on the properties that will be plotted
    
    Parameters
    ----------
    properties: list of strings
        Properties to plot.

    datasets: list of strings
        List of saved properties.

    Returns
    -------
    new_datasets: list of strings
        New list of properties that will be read.
    '''
    
    new_datasets = []
    for prop in properties:
        if (prop in datasets) and (prop not in new_datasets):
            new_datasets.append(prop)
        if (prop == "lithology") and ("strain" not in new_datasets):
            new_datasets.append("strain")
        if (prop == "temperature_anomaly") and ("temperature" not in new_datasets):
            new_datasets.append("temperature")
        if(prop == 'melt'):
            new_datasets.append('melt')
        if(prop == 'incremental_melt'):
            new_datasets.append('incremental_melt')
        if (prop == "lithology" or prop == 'temperature_anomlay') and ("density" not in new_datasets):
            new_datasets.append("density")
            
    return new_datasets


def plot_tracked_particles(trackdataset, ax, i, plot_other_particles=False, color_other_particles='xkcd:black', size_other_particles=5):
    """
    Plot tracked particles in the subplot ax

    Parameters
    ----------
    trackdataset : xarray.Dataset
        Dataset containing the tracked particles
    ax : matplotlib.axes._subplots.AxesSubplot
        Subplot to plot the tracked particles
    i : int
        Current time step
    """
    x_track = trackdataset.xtrack.values[::-1]
    z_track = trackdataset.ztrack.values[::-1]
    P = trackdataset.ptrack.values[::-1]
    T = trackdataset.ttrack.values[::-1]
    time = trackdataset.time.values[::-1]
    steps = trackdataset.step.values[::-1]
    n = int(trackdataset.ntracked.values)
    nTotal = np.size(x_track)
    steps = nTotal//n

    x_track = np.reshape(x_track,(steps,n))
    z_track = np.reshape(z_track,(steps,n))
    P = np.reshape(P,(steps,n))
    T = np.reshape(T,(steps,n))
    particles_layers = trackdataset.particles_layers.values[::-1] #code of the tracked layers

    asthenosphere_code = 0
    mantle_lithosphere1_code = 1
    seed_code = 2
    mantle_lithosphere2_code = 3
    lower_crust_code = 4
    upper_crust_code = 5

    T_initial = T[0]
    P_initial = P[0]

    if(asthenosphere_code in particles_layers):
        cond_ast = particles_layers == asthenosphere_code
        particles_ast = particles_layers[cond_ast]

        cond_ast2plot = P_initial <= 4000 #only plot asthenosphere particles with initial pressure less than 4000 MPa
        plot_asthenosphere_particles = True
    else:
        plot_asthenosphere_particles = False
        cond_ast2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if((mantle_lithosphere1_code in particles_layers) | (mantle_lithosphere2_code in particles_layers)):
        cond_mlit = (particles_layers == mantle_lithosphere1_code) | (particles_layers == mantle_lithosphere2_code)
        particles_mlit = particles_layers[cond_mlit]

        T_initial_mlit = T_initial[cond_mlit] #initial temperature of lithospheric mantle particles
        T_initial_mlit_sorted = np.sort(T_initial_mlit)

        Ti_mlit_max = np.max(T_initial_mlit_sorted)
        mid_index = len(T_initial_mlit_sorted)//2
        Ti_mlit_mid = T_initial_mlit_sorted[mid_index]
        Ti_mlit_min = np.min(T_initial_mlit_sorted)

        cond_mantle_lithosphere2plot = (T_initial == Ti_mlit_min) | (T_initial == Ti_mlit_mid) | (T_initial == Ti_mlit_max)

        plot_mantle_lithosphere_particles = True

        dict_mlit_markers = {Ti_mlit_max: '*',
                            Ti_mlit_mid: '^',
                            Ti_mlit_min: 'D'}

        dict_mlit_colors = {Ti_mlit_min: 'xkcd:cerulean blue',
                            Ti_mlit_mid: 'xkcd:scarlet',
                            Ti_mlit_max: 'xkcd:dark green'}
    else:
        plot_mantle_lithosphere_particles = False
        cond_mantle_lithosphere2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if(lower_crust_code in particles_layers):
        cond_crust = particles_layers == lower_crust_code
        particles_crust = particles_layers[cond_crust]

        T_initial_lower_crust = T_initial[cond_crust] #initial temperature of crustal particles
        T_initial_lower_crust_sorted = np.sort(T_initial_lower_crust)

        Ti_lower_crust_max = np.max(T_initial_lower_crust_sorted)
        mid_index = len(T_initial_lower_crust_sorted)//2
        Ti_lower_crust_mid = T_initial_lower_crust_sorted[mid_index]
        Ti_lower_crust_min = np.min(T_initial_lower_crust_sorted)

        cond_lower_crust2plot = (T_initial == Ti_lower_crust_min) | (T_initial == Ti_lower_crust_mid) | (T_initial == Ti_lower_crust_max)
        plot_lower_crust_particles = True
    else:
        plot_lower_crust_particles = False
        cond_lower_crust2plot = np.arange(0, n, 1) == np.arange(0,n,1) + 1

    # linewidth = 0.85
    markersize = 8
    color_lower_crust='xkcd:brown'
    h_air = 40.0

    for particle, particle_layer in zip(range(n), particles_layers):
        #Plot particles in prop subplot

        if((plot_lower_crust_particles == True) & (particle_layer == lower_crust_code)): #crustal particles
            # print(particle_layer)
            if(cond_lower_crust2plot[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air, '.', color=color_lower_crust, markersize=markersize-2, zorder=61)


        if((plot_mantle_lithosphere_particles == True) & ((particle_layer == mantle_lithosphere1_code) | (particle_layer == mantle_lithosphere2_code))): #lithospheric mantle particles
            if(cond_mantle_lithosphere2plot[particle]==True):
                # print(f"Particle: {particle}, Layer: {particle_layer}, T_initial: {T_initial[particle]}")
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air,
                            dict_mlit_markers[T_initial[particle]],
                            color=dict_mlit_colors[T_initial[particle]],
                            markersize=markersize-2, zorder=61)
            else:
                if(plot_other_particles == True):
                    ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air, '.', color=color_other_particles, markersize=size_other_particles, zorder=60)

        if((plot_asthenosphere_particles == True) & (particle_layer == asthenosphere_code)):
            if(cond_ast2plot[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air, '.', color='xkcd:violet', markersize=markersize-6, zorder=60)


def plot_tracked_particles_depth_coded(trackdataset, ax, i, hcrust=35.0e3, markersize=4,
                                       plot_lower_crust_particles=False, plot_mantle_lithosphere_particles=True, plot_asthenosphere_particles=True,
                                       color_lower_crust='xkcd:brown',
                                       color_mlit_upper='xkcd:cerulean blue', color_mlit_intermediate='xkcd:scarlet', color_mlit_lower='xkcd:dark green', particles_alpha=0.7):
    """
    Plot tracked particles in the subplot ax

    Parameters
    ----------
    trackdataset : xarray.Dataset
        Dataset containing the tracked particles
    ax : matplotlib.axes._subplots.AxesSubplot
        Subplot to plot the tracked particles
    i : int
        Index of the current step
    hcrust : float
        Crustal thickness in m
    markersize : float
        Size of the markers
    plot_lower_crust_particles : bool
        Whether to plot lower crust particles
    plot_mantle_lithosphere_particles : bool
        Whether to plot mantle lithosphere particles
    plot_asthenosphere_particles : bool
        Whether to plot asthenosphere particles
    color_lower_crust : str
        Color of the lower crust particles
    color_mlit_upper : str
        Color of the upper mantle lithosphere particles
    color_mlit_intermediate : str
        Color of the intermediate mantle lithosphere particles
    color_mlit_lower : str
        Color of the lower mantle lithosphere particles
    particles_alpha : float
        Alpha transparency of the particles
    """
    x_track = trackdataset.xtrack.values[::-1]
    z_track = trackdataset.ztrack.values[::-1] #the values of z are negative, so we need to multiply by -1 to plot them in the correct position
    P = trackdataset.ptrack.values[::-1]
    T = trackdataset.ttrack.values[::-1]
    time = trackdataset.time.values[::-1]
    steps = trackdataset.step.values[::-1]
    n = int(trackdataset.ntracked.values)
    nTotal = np.size(x_track)
    steps = nTotal//n

    x_track = np.reshape(x_track,(steps,n))
    z_track = np.reshape(z_track,(steps,n))
    P = np.reshape(P,(steps,n))
    T = np.reshape(T,(steps,n))
    particles_layers = trackdataset.particles_layers.values[::-1] #code of the tracked layers

    asthenosphere_code = 0
    mantle_lithosphere1_code = 1
    seed_code = 2
    mantle_lithosphere2_code = 3
    lower_crust_code = 4
    upper_crust_code = 5

    T_initial = T[0]
    P_initial = P[0]
    z_initial = z_track[0]
    x_initial = x_track[0]
    
    # color_lower_crust='xkcd:brown'
    # color_mlit_upper = 'xkcd:cerulean blue'
    # color_mlit_intermediate = 'xkcd:scarlet'
    # color_mlit_lower = 'xkcd:violet'

    h_air = 40.0e3 #m
    h_lithosphere = 120.0e3 #m
    #thickness of mantle lithosphere sections
    thickness_upper = 30.0e3
    thickness_intermediate = 30.0e3
    thickness_lower = h_lithosphere - (hcrust + thickness_upper + thickness_intermediate)

    zb_crust_and_air = -1.0 * (hcrust + h_air) #m z depth of the bottom of the crust considering air layer
    zb_upper = zb_crust_and_air - thickness_upper #m. The values of z_initial are negative, so we need to subtract
    zb_intermediate = zb_upper - thickness_intermediate #m. The values of z_initial are negative, so we need to subtract
    zb_lower = -1.0 * (h_lithosphere + h_air) #m z depth of the bottom of the lithosphere considering air layer

    if(asthenosphere_code in particles_layers):
        cond_ast = particles_layers == asthenosphere_code
        particles_ast = particles_layers[cond_ast]

        cond_ast2plot = P_initial <= 4000 #only plot asthenosphere particles with initial pressure less than 4000 MPa

        cond_upper_2plot_ast = (P_initial <= 4000) & (z_initial < zb_crust_and_air) & (z_initial >= zb_upper) #condition for upper particles - the values of z_initial are negative, so the logical condition is inverted
        cond_intermediate_2plot_ast = (P_initial <= 4000) & (z_initial < zb_upper) & (z_initial >= zb_intermediate) #condition for intermediate particles
        cond_lower_2plot_ast = (P_initial <= 4000) & (z_initial < zb_intermediate) & (z_initial >= zb_lower) #condition for lower particles

    else:
        plot_asthenosphere_particles = False
        cond_ast2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if((mantle_lithosphere1_code in particles_layers) | (mantle_lithosphere2_code in particles_layers)):
        cond_mlit = (particles_layers == mantle_lithosphere1_code) | (particles_layers == mantle_lithosphere2_code)
        particles_mlit = particles_layers[cond_mlit]

        # x_initial_mlit = x_initial[cond_mlit] #initial x position of lithospheric mantle particles
        # z_initial_mlit = z_initial[cond_mlit] #initial depth of lithospheric mantle particles

        cond_upper_2plot_mlit = (z_initial < zb_crust_and_air) & (z_initial >= zb_upper) #condition for upper particles - the values of z_initial are negative, so the logical condition is inverted
        cond_intermediate_2plot_mlit = (z_initial < zb_upper) & (z_initial >= zb_intermediate) #condition for intermediate particles
        cond_lower_2plot_mlit = (z_initial < zb_intermediate) & (z_initial >= zb_lower) #condition for lower particles
        
    else:
        plot_mantle_lithosphere_particles = False
        cond_mantle_lithosphere2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if(lower_crust_code in particles_layers):
        cond_crust = particles_layers == lower_crust_code
        particles_crust = particles_layers[cond_crust]

        T_initial_lower_crust = T_initial[cond_crust] #initial temperature of crustal particles
        T_initial_lower_crust_sorted = np.sort(T_initial_lower_crust)

        Ti_lower_crust_max = np.max(T_initial_lower_crust_sorted)
        mid_index = len(T_initial_lower_crust_sorted)//2
        Ti_lower_crust_mid = T_initial_lower_crust_sorted[mid_index]
        Ti_lower_crust_min = np.min(T_initial_lower_crust_sorted)

        cond_lower_crust2plot = (T_initial == Ti_lower_crust_min) | (T_initial == Ti_lower_crust_mid) | (T_initial == Ti_lower_crust_max)

    else:
        plot_lower_crust_particles = False
        cond_lower_crust2plot = np.arange(0, n, 1) == np.arange(0,n,1) + 1


    for particle, particle_layer in zip(range(n), particles_layers):
        #Plot particles in prop subplot

        if((plot_lower_crust_particles == True) & (particle_layer == lower_crust_code)): #crustal particles
            if(cond_lower_crust2plot[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_lower_crust, markersize=markersize, zorder=61, alpha=particles_alpha)

        if((plot_mantle_lithosphere_particles == True) & ((particle_layer == mantle_lithosphere1_code) | (particle_layer == mantle_lithosphere2_code))): #lithospheric mantle particles
            # if(cond_mantle_lithosphere2plot[particle]==True):
                # print(f"Particle: {particle}, Layer: {particle_layer}, T_initial: {T_initial[particle]}")
            if(cond_upper_2plot_mlit[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_mlit_upper, markersize=markersize, zorder=60, alpha=particles_alpha)
            elif(cond_intermediate_2plot_mlit[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_mlit_intermediate, markersize=markersize, zorder=60, alpha=particles_alpha)
            elif(cond_lower_2plot_mlit[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_mlit_lower, markersize=markersize, zorder=60, alpha=particles_alpha)
            # else:
            #     if(plot_other_particles == True):
            #         ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air, '.', color=color_other_particles, markersize=size_other_particles, zorder=60)

        if((plot_asthenosphere_particles == True) & (particle_layer == asthenosphere_code)):
            if(cond_upper_2plot_ast[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_mlit_upper, markersize=markersize, zorder=60, alpha=particles_alpha)
            elif(cond_intermediate_2plot_ast[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_mlit_intermediate, markersize=markersize, zorder=60, alpha=particles_alpha)
            elif(cond_lower_2plot_ast[particle] == True):
                ax.plot(x_track[i, particle]/1.0e3, z_track[i, particle]/1.0e3+h_air/1.0e3, '.', color=color_mlit_lower, markersize=markersize, zorder=60, alpha=particles_alpha)


def plot_ptt_paths_depth_coded_frame(trackdataset, ax, i, current_time, hcrust=35.0e3, alpha_lines=1.0, markersize=6, linewidth=0.35,
                                     plot_lower_crust_particles=False, plot_mantle_lithosphere_particles=True, plot_asthenosphere_particles=True,
                                     color_lower_crust='xkcd:brown',
                                     color_mlit_upper='xkcd:cerulean blue', color_mlit_intermediate='xkcd:scarlet', color_mlit_lower='xkcd:dark green',
                                     plot_steps=False,
                                     alpha=0.7):
    """
    Plot PTt path of tracked particles in the subplot ax

    Parameters
    ----------
    trackdataset : xarray.Dataset
        Dataset containing the tracked particles
    ax : matplotlib.axes._subplots.AxesSubplot
        Subplot to plot the PTt path
    i : int
        Current step index to plot the PTt path
    current_time : float
        Current time to plot the PTt path
    plot_lower_crust_particles : bool
        Whether to plot lower crust particles
    plot_mantle_lithosphere_particles : bool
        Whether to plot mantle lithosphere particles
    plot_asthenosphere_particles : bool
        Whether to plot asthenosphere particles
    color_lower_crust : str
        Color for lower crust particles
    color_mlit_upper : str
        Color for upper mantle lithosphere particles
    color_mlit_intermediate : str
        Color for intermediate mantle lithosphere particles
    color_mlit_lower : str
        Color for lower mantle lithosphere particles
    plot_steps : bool
        Whether to plot steps as black dots along the PTt path
    linewidth : float
        Width of the lines to plot
    markersize : float
        Size of the markers to plot
    alpha : float
        Transparency of the lines to plot
    """
    x_track = trackdataset.xtrack.values[::-1]
    z_track = trackdataset.ztrack.values[::-1]
    P = trackdataset.ptrack.values[::-1]
    T = trackdataset.ttrack.values[::-1]
    time = trackdataset.time.values[::-1]
    steps = trackdataset.step.values[::-1]
    n = int(trackdataset.ntracked.values)
    nTotal = np.size(x_track)
    steps = nTotal//n

    x_track = np.reshape(x_track,(steps,n))
    z_track = np.reshape(z_track,(steps,n))
    P = np.reshape(P,(steps,n))
    T = np.reshape(T,(steps,n))
    particles_layers = trackdataset.particles_layers.values[::-1] #code of the tracked layers
    
    z_initial = z_track[0]
    x_initial = x_track[0]
    
    # color_lower_crust='xkcd:brown'
    # color_mlit_upper = 'xkcd:cerulean blue'
    # color_mlit_intermediate = 'xkcd:scarlet'
    # color_mlit_lower = 'xkcd:violet'

    h_air = 40.0e3 #m
    h_lithosphere = 120.0e3 #m
    #thickness of mantle lithosphere sections
    thickness_upper = 30.0e3
    thickness_intermediate = 30.0e3
    thickness_lower = h_lithosphere - (hcrust + thickness_upper + thickness_intermediate)

    zb_crust_and_air = -1.0 * (hcrust + h_air) #m z depth of the bottom of the crust considering air layer
    zb_upper = zb_crust_and_air - thickness_upper #m. The values of z_initial are negative, so we need to subtract
    zb_intermediate = zb_upper - thickness_intermediate #m. The values of z_initial are negative, so we need to subtract
    zb_lower = -1.0 * (h_lithosphere + h_air) #m z depth of the bottom of the lithosphere considering air layer

    asthenosphere_code = 0
    mantle_lithosphere1_code = 1
    seed_code = 2
    mantle_lithosphere2_code = 3
    lower_crust_code = 4
    upper_crust_code = 5
    T_initial = T[0]
    P_initial = P[0]

    if(asthenosphere_code in particles_layers):
        cond_ast = particles_layers == asthenosphere_code
        particles_ast = particles_layers[cond_ast]

        cond_ast2plot = P_initial <= 4000 #only plot asthenosphere particles with initial pressure less than 4000 MPa

        cond_upper_2plot_ast = (P_initial <= 4000) & (z_initial < zb_crust_and_air) & (z_initial >= zb_upper) #condition for upper particles - the values of z_initial are negative, so the logical condition is inverted
        cond_intermediate_2plot_ast = (P_initial <= 4000) & (z_initial < zb_upper) & (z_initial >= zb_intermediate) #condition for intermediate particles
        cond_lower_2plot_ast = (P_initial <= 4000) & (z_initial < zb_intermediate) & (z_initial >= zb_lower) #condition for lower particles

    else:
        plot_asthenosphere_particles = False
        cond_ast2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if((mantle_lithosphere1_code in particles_layers) | (mantle_lithosphere2_code in particles_layers)):
        cond_mlit = (particles_layers == mantle_lithosphere1_code) | (particles_layers == mantle_lithosphere2_code)
        particles_mlit = particles_layers[cond_mlit]

        # x_initial_mlit = x_initial[cond_mlit] #initial x position of lithospheric mantle particles
        # z_initial_mlit = z_initial[cond_mlit] #initial depth of lithospheric mantle particles

        cond_upper_2plot_mlit = (z_initial < zb_crust_and_air) & (z_initial >= zb_upper) #condition for upper particles - the values of z_initial are negative, so the logical condition is inverted
        cond_intermediate_2plot_mlit = (z_initial < zb_upper) & (z_initial >= zb_intermediate) #condition for intermediate particles
        cond_lower_2plot_mlit = (z_initial < zb_intermediate) & (z_initial >= zb_lower) #condition for lower particles

    else:
        plot_mantle_lithosphere_particles = False
        cond_mantle_lithosphere2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if(lower_crust_code in particles_layers):
        cond_crust = particles_layers == lower_crust_code
        particles_crust = particles_layers[cond_crust]

        T_initial_lower_crust = T_initial[cond_crust] #initial temperature of crustal particles
        T_initial_lower_crust_sorted = np.sort(T_initial_lower_crust)

        Ti_lower_crust_max = np.max(T_initial_lower_crust_sorted)
        mid_index = len(T_initial_lower_crust_sorted)//2
        Ti_lower_crust_mid = T_initial_lower_crust_sorted[mid_index]
        Ti_lower_crust_min = np.min(T_initial_lower_crust_sorted)

        cond_lower_crust2plot = (T_initial == Ti_lower_crust_min) | (T_initial == Ti_lower_crust_mid) | (T_initial == Ti_lower_crust_max)

    else:
        plot_lower_crust_particles = False
        cond_lower_crust2plot = np.arange(0, n, 1) == np.arange(0,n,1) + 1


    for particle, particle_layer in zip(range(n), particles_layers):
        #Plot particles in prop subplot

        if((plot_lower_crust_particles == True) & (particle_layer == lower_crust_code)):
            if(cond_lower_crust2plot[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_lower_crust, linewidth=linewidth, alpha=alpha, zorder=60) #PTt path
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_lower_crust, markersize=markersize) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=2, zorder=60)

        if((plot_mantle_lithosphere_particles == True) & ((particle_layer == mantle_lithosphere1_code) | (particle_layer == mantle_lithosphere2_code))): #lithospheric mantle particles
            if(cond_upper_2plot_mlit[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_mlit_upper, linewidth=linewidth, alpha=alpha, zorder=60) #PTt path
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_mlit_upper, markersize=markersize, zorder=61) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=2, zorder=60)
            if(cond_intermediate_2plot_mlit[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_mlit_intermediate, linewidth=linewidth, alpha=alpha, zorder=60)
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_mlit_intermediate, markersize=markersize, zorder=61) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle]/1.0e3, '.', color='xkcd:black', markersize=2, zorder=60)
            if(cond_lower_2plot_mlit[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_mlit_lower, linewidth=linewidth, alpha=alpha, zorder=60)
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_mlit_lower, markersize=markersize, zorder=61) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle]/1.0e3, '.', color='xkcd:black', markersize=2, zorder=60)

        if((plot_asthenosphere_particles == True) & (particle_layer == asthenosphere_code)):
            if(cond_upper_2plot_ast[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_mlit_upper, linewidth=linewidth-1, alpha=alpha, zorder=60)
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_mlit_upper, markersize=4, zorder=61) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=2, zorder=60)
            if(cond_intermediate_2plot_ast[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_mlit_intermediate, linewidth=linewidth-1, alpha=alpha, zorder=60)
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_mlit_intermediate, markersize=4, zorder=61) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=2, zorder=60)
            if(cond_lower_2plot_ast[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle]/1.0e3, '-', color=color_mlit_lower, linewidth=linewidth-1, alpha=alpha, zorder=60)
                ax.plot(T[i, particle], P[i, particle]/1.0e3, '.', color=color_mlit_lower, markersize=4, zorder=61) #current step
                #plotting points at each 5 Myr
                if(plot_steps==True):
                    for j in np.arange(0, current_time, 5):
                        idx = find_nearest(time, j)
                        ax.plot(T[idx, particle], P[idx, particle]/1.0e3, '.', color='xkcd:black', markersize=2, zorder=60)


def plot_ptt_paths_three_particles(trackdataset, ax, i, current_time, plot_other_particles=True, color_other_particles='xkcd:black', size_other_particles=0.7):
    """
    Plot PTt path of tracked particles in the subplot ax

    Parameters
    ----------
    trackdataset : xarray.Dataset
        Dataset containing the tracked particles
    ax : matplotlib.axes._subplots.AxesSubplot
        Subplot to plot the PTt path
    instants : list
        List of instants to plot the PTt path
    """
    x_track = trackdataset.xtrack.values[::-1]
    z_track = trackdataset.ztrack.values[::-1]
    P = trackdataset.ptrack.values[::-1]
    T = trackdataset.ttrack.values[::-1]
    time = trackdataset.time.values[::-1]
    steps = trackdataset.step.values[::-1]
    n = int(trackdataset.ntracked.values)
    nTotal = np.size(x_track)
    steps = nTotal//n

    x_track = np.reshape(x_track,(steps,n))
    z_track = np.reshape(z_track,(steps,n))
    P = np.reshape(P,(steps,n))
    T = np.reshape(T,(steps,n))
    particles_layers = trackdataset.particles_layers.values[::-1] #code of the tracked layers

    asthenosphere_code = 0
    mantle_lithosphere1_code = 1
    seed_code = 2
    mantle_lithosphere2_code = 3
    lower_crust_code = 4
    upper_crust_code = 5
    T_initial = T[0]
    P_initial = P[0]

    if(asthenosphere_code in particles_layers):
        cond_ast = particles_layers == asthenosphere_code
        particles_ast = particles_layers[cond_ast]

        cond_ast2plot = P_initial <= 4000 #only plot asthenosphere particles with initial pressure less than 4000 MPa
        plot_asthenosphere_particles = True
    else:
        plot_asthenosphere_particles = False
        cond_ast2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if((mantle_lithosphere1_code in particles_layers) | (mantle_lithosphere2_code in particles_layers)):
        cond_mlit = (particles_layers == mantle_lithosphere1_code) | (particles_layers == mantle_lithosphere2_code)
        particles_mlit = particles_layers[cond_mlit]

        T_initial_mlit = T_initial[cond_mlit] #initial temperature of lithospheric mantle particles
        T_initial_mlit_sorted = np.sort(T_initial_mlit)

        Ti_mlit_max = np.max(T_initial_mlit_sorted)
        mid_index = len(T_initial_mlit_sorted)//2
        Ti_mlit_mid = T_initial_mlit_sorted[mid_index]
        Ti_mlit_min = np.min(T_initial_mlit_sorted)

        cond_mantle_lithosphere2plot = (T_initial == Ti_mlit_min) | (T_initial == Ti_mlit_mid) | (T_initial == Ti_mlit_max)

        plot_mantle_lithosphere_particles = True

        dict_mantle_lithosphere_markers = {Ti_mlit_max: '*',
                            Ti_mlit_mid: '^',
                            Ti_mlit_min: 'D'}

        dict_mantle_lithosphere_colors = {Ti_mlit_min: 'xkcd:cerulean blue',
                            Ti_mlit_mid: 'xkcd:scarlet',
                            Ti_mlit_max: 'xkcd:dark green'}
    else:
        plot_mantle_lithosphere_particles = False
        cond_mantle_lithosphere2plot = np.arange(0, n, 1) == np.arange(0, n, 1) + 1

    if(lower_crust_code in particles_layers):
        cond_crust = particles_layers == lower_crust_code
        particles_crust = particles_layers[cond_crust]

        T_initial_crust = T_initial[cond_crust] #initial temperature of crustal particles
        T_initial_crust_sorted = np.sort(T_initial_crust)

        Ti_crust_max = np.max(T_initial_crust_sorted)
        mid_index = len(T_initial_crust_sorted)//2
        Ti_crust_mid = T_initial_crust_sorted[mid_index]
        Ti_crust_min = np.min(T_initial_crust_sorted)

        cond_crust2plot = (T_initial == Ti_crust_min) | (T_initial == Ti_crust_mid) | (T_initial == Ti_crust_max)
        plot_lower_crust_particles = True
    else:
        plot_lower_crust_particles = False
        cond_crust2plot = np.arange(0, n, 1) == np.arange(0,n,1) + 1

    linewidth = 0.85
    markersize = 8
    color_lower_crust='xkcd:brown'
    ast_color = 'xkcd:violet'

    for particle, particle_layer in zip(range(n), particles_layers):
        #Plot particles in prop subplot

        if((plot_lower_crust_particles == True) & (particle_layer == lower_crust_code)):
            if(cond_crust2plot[particle] == True):
                ax.plot(T[:i, particle], P[:i, particle], '-', color=color_lower_crust, linewidth=linewidth, alpha=1.0, zorder=60) #PTt path
                #plotting points at each 5 Myr
                for j in np.arange(0, current_time, 5):
                    idx = find_nearest(time, j)
                    ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=2, zorder=60)

        if((plot_mantle_lithosphere_particles == True) & ((particle_layer == mantle_lithosphere1_code) | (particle_layer == mantle_lithosphere2_code))): #lithospheric mantle particles

                    if(cond_mantle_lithosphere2plot[particle]==True):
                        # print(f"Particle: {particle}, Layer: {particle_layer}, T_initial: {T_initial[particle]}"
                        ax.plot(T[i, particle], P[i, particle],
                                    dict_mantle_lithosphere_markers[T_initial[particle]],
                                    color=dict_mantle_lithosphere_colors[T_initial[particle]], markersize=10, zorder=61) #current PTt point

                        ax.plot(T[:i, particle], P[:i, particle], '-', color=dict_mantle_lithosphere_colors[T_initial[particle]], linewidth=linewidth, alpha=0.8, zorder=61) #PTt path
                        #plotting points at each 5 Myr
                        
                        for j in np.arange(0, current_time, 5):
                            idx = find_nearest(time, j)
                            ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=2, zorder=60)
                    else:
                        if(plot_other_particles == True): #plotting the other lithospheric mantle particles instead only 3
                            ax.plot(T[i, particle], P[i, particle],
                                        '.', color='xkcd:black',
                                        markersize=int(markersize/2), zorder=60) #current PTt point
                            ax.plot(T[:i, particle], P[:i, particle],
                                        '-', color='xkcd:black',
                                        linewidth=0.1, alpha = 1.0, zorder=60) #PTt path
                                #plotting points at each 5 Myr
                            for j in np.arange(0, current_time, 5):
                                idx = find_nearest(time, j)
                                ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=0.7, zorder=60)

        if((plot_asthenosphere_particles == True) & (particle_layer == asthenosphere_code)):
                    if(cond_ast2plot[particle] == True):
                        ax.plot(T[i, particle], P[i, particle], '.', color=ast_color, markersize=int(markersize/2), zorder=60)
                        ax.plot(T[:i, particle], P[:i, particle], '-', color=ast_color, linewidth=linewidth-1, alpha=1.0, zorder=60) #PTt path
                        #plotting points at each 5 Myr
                        for j in np.arange(0, current_time, 5):
                            idx = find_nearest(time, j)
                            ax.plot(T[idx, particle], P[idx, particle], '.', color='xkcd:black', markersize=0.5, zorder=59)

def plot_property(dataset, prop, xlims, ylims, model_path,
                fig,
                ax,
                print_time=True,
                correction_factor=0.0,
                plot_isotherms=True, isotherms=[400, 600, 800, 1000, 1300],
                topo_from_density=True,
                plot_particles=False,
                particle_size=0.2,
                particle_marker="o",
                plot_colorbar=True,
                bbox_to_anchor=(0.85,#horizontal position respective to parent_bbox or "loc" position
                                0.3,# vertical position
                                0.12,# width
                                0.35),
                ncores=None,
                step_plot=1,
                plot_melt=False,
                color_incremental_melt = 'xkcd:bright pink',
                color_depleted_mantle='xkcd:bright purple',):
    '''
    Plot data from mandyoc according to a given property and domain limits.

    Parameters
    ----------
    dataset: :class:`xarray.Dataset`
        Dataset containing mandyoc data for a single time step, read from the .nc files created by tapioca
        (strain.nc, temperature.nc, density.nc, Phi.nc, ...): data (z, x), x from 0 to Lx and z from -Lz to 0 [m].
    
    prop: str
        Property from mandyoc.

    xlims: list
        List with the limits of x axis
        
    ylims: list
        List with the limits of y axis

    model_path: str
        Path to model (with surface.nc and particles_trajectories.nc, used for prop='surface' with
        topo_from_density=False and for plot_particles=True)

    output_path: str
        Path to save outputs
        
    save_frames: bool
        True to save frame by frames
        False to do not save the frames
    '''
    
    props_label = {'density':              r'$\mathrm{[kg/m^3]}$',
                   'radiogenic_heat':       'log(W/kg)',
                   'lithology':            r'log$(\epsilon_{II})$',
                   'pressure':              'P [GPa]',
                   'strain':               r'Accumulated strain [$\varepsilon$]',
                   'strain_rate':          r'log($\dot{\varepsilon}$)',
#                    'strain_rate':          r'$\dot{\varepsilon}$',
                   'temperature':          r'$^{\circ}\mathrm{[C]}$',
                   'temperature_anomaly':  r'Temperature anomaly $^{\circ}\mathrm{[C]}$',
                   'topography':            'Topography [km]',
                   'viscosity':             'log(Pa.s)'
                   }
    
    props_cmap = {'density': 'viridis',
                  'radiogenic_heat': 'inferno',
                  'lithology': 'viridis',
                  'pressure': 'viridis',
#                   'strain': 'viridis', #Default. Comment this line and uncomment one of the options bellow
#                   'strain': 'cividis',
#                   'strain': 'Greys',
                  'strain': 'inferno',
#                   'strain': 'magma',
                  'strain_rate': 'viridis',
                  'temperature': 'viridis',
                  'temperature_anomaly': 'RdBu_r',
                  'topography': '',
                  'viscosity': 'viridis'
                   }

    #limits of colorbars
    vals_minmax = {'density':             [0.0, 3378.],
                   'radiogenic_heat':     [1.0E-13, 1.0E-9],
                   'lithology':           [None, None],
                   'pressure':            [-1.0E-3, 1.0],
                   'strain':              [None, None],
                   'strain_rate':         [1.0E-18, 1.0E-14],
                #    'strain_rate':         [1.0E-20, 1.0E-15],
#                    'strain_rate':         [np.log10(1.0E-19), np.log10(1.0E-14)],
                   'temperature':         [0, 1600],
                   'temperature_anomaly': [-150, 150],
                   'surface':             [-6, 6],
                   'viscosity':           [1.0E18, 1.0E25],
#                    'viscosity':           [np.log10(1.0E18), np.log10(1.0E25)]
                  }

    Nx = int(dataset.sizes["x"])
    Nz = int(dataset.sizes["z"])
    Lx = float(dataset.x.values[-1] - dataset.x.values[0])
    Lz = float(dataset.z.values[-1] - dataset.z.values[0])
    instant = np.round(float(dataset.time), 2)

    #tapioca keeps the values in the air (density < 1200 kg/m3): they are set to zero where they are plotted
    air = dataset.density.values < 1200.0 if "density" in dataset else np.zeros((Nz, Nx), dtype=bool)

    def tapioca_surface():
        '''surface.nc of the time step of dataset [km], in the x of the mesh (surface.nc has more points)'''
        surface = xr.open_dataset(f"{model_path}/surface.nc")
        surface = surface.surface.isel(time=int(np.where(surface.step.values == int(dataset.step))[0][0]))
        return np.interp(xi*1.0e3, surface.x.values, surface.values)
    
    xi = np.linspace(0, Lx/1000, Nx)
    zi = np.linspace(-Lz/1000+40, 0+40, Nz) #km, +40 to compensate the air layer above sea level
    xx, zz = np.meshgrid(xi, zi)
    
    #creating Canvas
    label_size=12
    plt.rc('xtick', labelsize = label_size)
    plt.rc('ytick', labelsize = label_size)
    
    #plot Time in Myr
    # ax.text(0.68, 1.035, ' {:01} Myr'.format(instant), bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=ax.transAxes)
    if(print_time == True):
        ax.text(0.78, 1.035, ' {:01} Myr'.format(instant), bbox=dict(facecolor='white', edgecolor='white', alpha=0.0), fontsize = 14, zorder=52, transform=ax.transAxes)
    
    val_minmax = vals_minmax[prop]
    
    if(plot_isotherms == True and prop != 'surface'): #Plot isotherms
        Temperi = dataset.temperature.values
        
        isot_colors = []
        for isotherm in isotherms:
            isot_colors.append('red')
            
        cs = ax.contour(xx, zz, Temperi, 100, levels=isotherms, colors=isot_colors)
        
        # if(instant == instants[0]):
        #     fmt = {}
        #     for level, isot in zip(cs.levels, isotherms):
        #         fmt[level] = str(level) + r'$^{\circ}$C'

        #     ax.clabel(cs, cs.levels, fmt=fmt, inline=True, use_clabeltext=True)

    if(plot_melt == True and prop != 'surface'):

        #TRYING TO PLOT MELT
        melt = np.where(air, 0.0, dataset.Phi.values) #depleted mantle

        levels_contourf = np.arange(0.02, 0.5, 0.02)
        for i in levels_contourf:
            alpha = 0.9/levels_contourf[-1]*i + 0.1

            cs0 = ax.contourf(xx, zz, melt, levels=[i, i+0.02], colors='xkcd:black', alpha=alpha, zorder=30)


        levels_melt = [0.1, 0.2] #np.arange(0.2, 0.7, 0.4)#[0.2, 0.6]
        cs = ax.contour(xx,
                    zz,
                    melt,
                    levels = levels_melt,
                    colors='xkcd:blue',
                    # cmap = 'inferno',
                    zorder=30,
                    )
        
        
        incremental_melt = np.where(air, 0.0, dataset.dPhi.values)
        scale=1.0e4
        
        levels_incremental_melt = [5.0e-10, 1.0e-6]

        ax.contourf(xx,
                    zz,
                    incremental_melt/scale,#, 'dashdot', 'dashed'],
                    levels = levels_incremental_melt,
                    colors='xkcd:bright pink',#['xkcd:bright pink', 'xkcd:pink'],
                    # colors=colors,
                    alpha=0.4,
                    zorder=30)
        
        ax.contour(xx,
                    zz,
                    incremental_melt/scale,
                    linestyles=['solid'],#, 'dashdot', 'dashed'],
                    levels = [5.0E-10,1.0E-8],
                    colors='xkcd:bright pink',#['xkcd:bright pink', 'xkcd:pink'],
                    # colors=colors,
                    alpha=1.0,
                    zorder=30)
        
    #dealing with special data
    if(prop == 'lithology'):
        with np.errstate(divide='ignore'): #log10(0) in the air
            data = np.log10(np.where(air, 0.0, dataset['strain'].values)) #shaded log of the strain, as in the frame scripts
        
    elif(prop == 'temperature_anomaly'):
        #removing horizontal mean temperature
        B = dataset['temperature'].values #shape: (Nz, Nx)
        C = np.mean(B, axis=1) #shape: Nz
        D = (B - C[:, None]) #(Nz, Nx)
        data = D
        
    elif(prop == 'surface'):
        # print('Dealing with data')
        # topo_from_density = True
        # topo_from_density = False
        
        if(topo_from_density == True):
            Rhoi = dataset.density.values
            interfaces=[2900, 3365]
            ##Extract layer topography
            z = np.linspace(Lz/1000.0, 0, Nz)
            Z = np.linspace(Lz/1000.0, 0, 8001) #zi
            x = np.linspace(Lx/1000.0, 0, Nx)

            topo_interface = _extract_interface(z, Z, Nx, Rhoi, 300.) #200 kg/m3 = air/crust interface
            
            condx = (xi >= 100) & (xi <= 600)
            z_mean = np.mean(topo_interface[condx])
            
            topo_interface -= np.abs(z_mean)
            topo_interface = -1.0*topo_interface

            data = topo_interface
        else:
            surface = tapioca_surface()
            condx = (xi >= 100) & (xi <= 600)
            z_mean = np.mean(surface[condx])/1.0e3 + 40.0

            data = surface/1.0e3 + 40.0 + np.abs(z_mean) #km + air layer correction
            
            
    elif(prop == 'pressure'):
        data = dataset[prop].values/1.0E9 #GPa
        
    else:
        data = dataset[prop].values
        
    if(prop == 'strain_rate' or prop == 'radiogenic_heat' or prop == 'strain_rate' or prop == 'viscosity'): #properties that need a lognorm colorbar
        im = ax.imshow(data,
                       cmap = props_cmap[prop],
                       origin='lower',
                       extent = (0, Lx / 1.0E3, -Lz / 1.0E3 + 40, 0 + 40),
                       norm = LogNorm(vmin=val_minmax[0], vmax=val_minmax[1]),
#                        vmin=val_minmax[0], vmax=val_minmax[1],
                       aspect = 'auto')
        
        #creating colorbar
        if(plot_colorbar):
            #fixed size [inches], whatever the size of the figure
            add_colorbar(fig, ax, im, props_label[prop], position=(0.7, 0.3), loc='lower left', fmt=_log_fmt)
    
    elif (prop == 'density' or prop == 'pressure' or prop == 'temperature' or prop == 'temperature_anomaly'): #properties that need a regular colorbar
        im = ax.imshow(data,
                       cmap = props_cmap[prop],
                       origin='lower',
                       extent = (0, Lx / 1.0E3, -Lz / 1.0E3 + 40, 0 + 40),
                       vmin=val_minmax[0], vmax=val_minmax[1],
                       aspect = 'auto')
        
        if(plot_colorbar):
            #precision of colorbar ticks
            if(prop == 'pressure'): 
                fmt = '%.2f'
            else:
                fmt = '%.0f'

            #fixed size [inches], whatever the size of the figure
            add_colorbar(fig, ax, im, props_label[prop], position=(0.05, 0.3), loc='lower left', fmt=fmt)
        
    elif(prop == 'strain'):
        im = ax.imshow(data,
                       cmap = props_cmap[prop],
                       origin='lower',
                       extent = (0, Lx / 1.0E3, -Lz / 1.0E3 + 40, 0 + 40),
                       vmin = float(data.min()),
                       vmax = float(data.max()),
                       aspect = 'auto')
        
    elif(prop == 'surface'):
        ax.plot(xi, data, alpha = 1, linewidth = 2.0, color = "blueviolet")
        
    elif(prop == 'lithology'): #shaded lithology plot
        cr = 255.
        color_uc = (228. / cr, 156. / cr, 124. / cr)
        color_lc = (240. / cr, 209. / cr, 188. / cr)
        color_lit = (155. / cr, 194. / cr, 155. / cr)
        color_ast = (207. / cr, 226. / cr, 205. / cr)

        # cr = 255.
        # color_sed = (241./cr,184./cr,68./cr)
        # color_dec = (137./cr,81./cr,151./cr)
        # color_uc = (228./cr,156./cr,124./cr)
        # color_lc = (240./cr,209./cr,188./cr)
        # color_lit = (155./cr,194./cr,155./cr)
        # color_mlit_uc = (180. / cr, 194. / cr, 162. / cr)
        # color_mlit_lc = (155. / cr, 194. / cr, 155. / cr)
        # color_ast = (207./cr,226./cr,205./cr)
        
        Rhoi = dataset.density.values
        # interfaces=[2900, 3365]
        
        # ##Extract layer topography
        # z = np.linspace(Lz/1000.0, 0, Nz)
        # Z = np.linspace(Lz/1000.0, 0, 8001) #zi
        # x = np.linspace(Lx/1000.0, 0, Nx)
            
        # topo_interface = _extract_interface(z, Z, Nx, Rhoi, 300.) #200 kg/m3 = air/crust interface
        # condx = (xi >= 100) & (xi <= 400)
        # z_mean = np.mean(topo_interface[condx])
        # topo_interface -= np.abs(z_mean)
        # topo_interface = -1.0*topo_interface
        
        #Density field
        ax.contourf(xx,
                    zz+correction_factor,
                    Rhoi,
                    levels = [200., 2750, 2900, 3365, 3900],
                    colors = [color_uc, color_lc, color_lit, color_ast],
                    # levels=[200.,2350,2450,2750,2900,3325,3355,3365,3378],
                    # colors=[color_sed,color_dec,color_uc,color_lc,color_lit,color_mlit_uc,color_mlit_lc,color_ast],
                    )
        #Strain shaded areas
        im=ax.imshow(data,
                     cmap = 'Greys',
                     origin = 'lower',
                     extent = (0, Lx / 1.0E3, -Lz / 1.0E3 + 40+correction_factor, 0 + 40+correction_factor),
                     # extent = (xlims[0], xlims[1], ylims[0], ylims[1]),
                     zorder = 50,
                     alpha = 0.2, vmin=-0.5,
                     vmax = 0.7,
                     aspect = 'auto')
        #legend box: fixed size [inches], whatever the size of the figure.
        #Its lower right corner is at the lower right corner of bbox_to_anchor (x, y, width, height in fraction of ax)
        if(plot_colorbar == True):
            plot_lithology_legend(ax, colors=[color_uc, color_lc, color_lit, color_ast], vmin=-0.5, vmax=0.7, alpha=0.2,
                                  position=(bbox_to_anchor[0] + bbox_to_anchor[2], bbox_to_anchor[1]), loc='lower right')
            
    
    if(plot_particles == True):
        if(prop != 'surface'):
            #particles of the time step of dataset from particles_trajectories.nc (tapioca); -999 where the particle does not exist
            particles = xr.open_dataset(f"{model_path}/particles_trajectories.nc")
            particles = particles.isel(time=int(np.where(particles.step.values == int(dataset.step))[0][0]))
            exists = particles["x"].values != -999.0
            data_x = particles["x"].values[exists]
            data_z = particles["z"].values[exists]
            data_lithology = particles["lithology"].values[exists]
            # ax.scatter(data_x/1000, data_z/1000, 2, c='xkcd:black', marker='.', zorder=30)

            cond_litho = data_lithology > 1
            cond_mb = data_lithology == 1
            cond_ast = data_lithology == 0

            # if(prop=='lithology'):
            #     color_litho = 'xkcd:bright pink'
            #     color_ast = 'xkcd:black'
            # else:
            #     color_litho = 'xkcd:bright green'
            #     color_ast = 'xkcd:black'

            color_litho = 'xkcd:black'
            # color_mb = 'xkcd:neon green'
            color_mb = 'xkcd:black'
            color_ast = 'xkcd:bright pink'
            color_mb = 'xkcd:black'
            # color_mb = 'xkcd:neon green'

            ax.plot(data_x[cond_litho][::step_plot]/1000, data_z[cond_litho][::step_plot]/1000+40, particle_marker, color=color_litho, markersize=particle_size, alpha=1.0, zorder=30)
            ax.plot(data_x[cond_mb][::step_plot]/1000, data_z[cond_mb][::step_plot]/1000+40, particle_marker, color=color_mb, markersize=particle_size, alpha=1.0, zorder=30)
            ax.plot(data_x[cond_ast][::step_plot]/1000, data_z[cond_ast][::step_plot]/1000+40, particle_marker, color=color_ast, markersize=particle_size, alpha=1.0, zorder=30)
            # ax.plot(data_x[cond_ast][::step_plot*4]/1000, data_z[cond_ast][::step_plot*4]/1000+40, particle_marker, color=color_ast, markersize=particle_size-0.95, alpha=1.0, zorder=30)
    
    #Fill above topography with white color
    if(prop != 'surface'):
        # topo_from_density = True
        # topo_from_density = False
        if(topo_from_density==True):
            Rhoi = dataset.density.values
            # interfaces=[2900, 3365]
            # ##Extract layer topography
            z = np.linspace(Lz/1000.0, 0, Nz)
            Z = np.linspace(Lz/1000.0, 0, 8001) #zi
            x = np.linspace(Lx/1000.0, 0, Nx)

            topo_interface = _extract_interface(z, Z, Nx, Rhoi, 200.) #200 kg/m3 = air/crust interface
            # condx = (xi >= 100) & (xi <= 400)
            # z_mean = np.mean(topo_interface[condx])
            # topo_interface -= np.abs(z_mean)
            topo_interface = -1.0*topo_interface + 40.0 - correction_factor
        else:
            topo_interface = tapioca_surface()/1.0e3 + 40.0

        xaux = xx[0]
        condaux = (xaux>=xlims[0]) & (xaux<=xlims[1])
        xaux = xaux[condaux]
        xaux[0] += 2
        xaux[-1] -= 2
        ax.fill_between(xaux, topo_interface[condaux], ylims[-1]-0.8, color='white', alpha=1.0, zorder=51)

        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        # ax.set_xlabel("Distance (km)", fontsize = label_size)
        # ax.set_ylabel("Depth (km)", fontsize = label_size)
        
    else:
        ax.grid('-k', alpha=0.7)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
