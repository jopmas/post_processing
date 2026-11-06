########################################################################
# BEFORE RUNNING THIS SCRIPT:                                          #
#  - CREATE THE .nc FILES WITH tapioca's meshes2netcdf.jl              #
#    (density.nc, temperature.nc, surface.nc)                          #
#  - RUN THIS SCRIPT IN THE SAME FOLDER WHERE THE .nc FILES ARE        #
#                                                                      #
# Forward models the vertical component of gravity (g_z) of the        #
# effective density in each time step with harmonica.prism_gravity and #
# writes _gz_evolution.nc, with g_z (time, x) in mGal.                 #
#                                                                      #
# - Reference: median, in each depth, of the far-field columns         #
#   (x_far_begin to x_far_end) in the same time step, without water.   #
# - Surface: the top of the rock follows surface.nc inside the cells.  #
# - Water (post-processing): fills the columns below the sea level,    #
#   which sink by Airy isostasy, so the mass of each column does not   #
#   change.                                                            #
# - Observation points: obs_height above the water or the rock, as in  #
#   surface free-air anomalies.                                        #
########################################################################

###########################################################################
# Decisions behind the corrections                                        #
#                                                                         #
# The numbers below were measured in Ho35S1_Clc40_Cmlit025_Ts0_1090 at    #
# about 40 Myr (except where noted). The same patterns were found in      #
# Ho35S1_dz05_Clc10_Cmlit025_Ts0_1080. They explain the choices: they are #
# not valid for other scenarios.                                          #
#                                                                         #
# Effective density                                                       #
#  density.nc has only the density of the layers, without the thermal     #
#  effect: rho_eff = rho*(1 - alpha*T). Only the contrasts matter, so the #
#  temperature of reference (0 oC) has little effect.                     #
#                                                                         #
# From 2D to 3D                                                           #
#  The model is 2D, so the prisms are long in the northing direction.     #
#  half_y of 1e7 m and 1e9 m differ by only 0.02 mGal, and half_y is      #
#  3.2e7 m with half_y_factor = 20.                                       #
#                                                                         #
# Reference density: far-field columns in the same time step              #
#  - Not the first time step: the far field changes during the model. The #
#    Moho rises about 1.6 km and the lithosphere cools, which creates a   #
#    plateau of about +84 mGal away from the rift.                        #
#  - Not one column at the side of the model: the temperature is fixed    #
#    there, while the interior cools. With x = 0, the plateau is also     #
#    about +100 mGal.                                                     #
#  - Not the median of the whole model: it is close to the far field      #
#    until 40 Myr, but differs by about 14 mGal at 60 Myr, when the rifted #
#    region is larger.                                                    #
#  - Not PREM (Dziewonski and Anderson, 1981): a 1D reference only shifts #
#    the profile. With the model extended in x, the difference to the     #
#    far-field reference is about -2700 mGal, with a standard deviation   #
#    of 1.8 mGal, and matches the integral of the density difference. The #
#    crust of PREM is also a global average, not a continental column.    #
#    Without the extension in x, a 1D reference also bends the profile    #
#    near the sides.                                                      #
#  Consequence: the far field is zero by construction, so the absolute    #
#  level and its slow changes cannot be compared between time steps.      #
#                                                                         #
# Surface inside the cells                                                #
#  The nodes around the surface mix rock and air, and a mask on the nodes #
#  leaves noise of up to about 1000 kg/m3 in one row, which changes g_z   #
#  by up to about 25 mGal. So the rock is extended up to the surface, and #
#  each cell gets the fraction of rock, water and air from surface.nc.    #
#                                                                         #
# Water with Airy isostasy                                                #
#  The model has air, not water, above the basin, so its surface is in    #
#  equilibrium without the water load. Water added without any adjustment #
#  is extra mass: the mean excess of the oceanic columns over the far     #
#  field goes from -44 to +166 mGal. Each column sinks so that its mass   #
#  does not change (the excess stays at -44 mGal). rho_m is the density   #
#  of the base of each column, the rock that leaves the model. The mean   #
#  depth of the ocean goes from 5.0 to 7.3 km.                            #
#  This is local isostasy, done after the simulation: there is no flexure #
#  and it does not change the evolution of the rift. The lithology of the #
#  model is not moved in the figures.                                     #
#  Sea level: the mean far-field surface in the first time step, computed #
#  for each scenario instead of a fixed value. Columns slightly below it  #
#  also get water.                                                        #
#                                                                         #
# Observation points                                                      #
#  At the water or rock surface, because the surface free-air anomalies   #
#  of real data (e.g. WGM2012) are given there. A flat level does not     #
#  compare with them, because it is farther from the sources where the    #
#  topography is low. The surface is rough and the mesh is coarse, so the #
#  peaks have steps. The profile is not filtered to the resolution of the #
#  data here.                                                             #
#                                                                         #
# Not included                                                            #
#  Sediments, surface processes and flexure. The peaks at the margins are #
#  still high (about 300 mGal at 20-30 Myr). In a test after the          #
#  simulation (t = 20, 40 and 60 Myr), most of the difference comes from  #
#  the scarp: without the rock above the sea level, the peaks are about   #
#  90-130 mGal, or about 150 mGal with isostatic rebound. A wedge of      #
#  sediments up to 8 km thick changes the peaks by about 20 mGal or less. #
#  Peaks of about 150 mGal and minimums of about -140 mGal remain.        #
###########################################################################

import numpy as np
import os
import sys
import xarray as xr
import harmonica as hm
import multiprocessing #needed to run pymp in macOS
multiprocessing.set_start_method('fork') #needed to run pymp in macOS arm architectures
import pymp

model_path = os.getcwd() # Get local file
model_name = model_path.split('/')[-1]

###########################################################################
# Settings                                                                #
###########################################################################

alpha = 3.28e-5 #thermal expansion coefficient [1/oC]

#Far-field columns: reference density and sea level
x_far_begin = 200.0e3 #m
x_far_end = 400.0e3 #m

water = True #fills the columns below the sea level with water, with Airy isostasy
rho_water = 1000.0 #kg/m3
rho_air = 1.0 #kg/m3

obs_height = 100.0 #m, above the water or the rock

half_y_factor = 20.0 #half length of the prisms in the northing direction, in units of Lx: extends the 2D model

stride = 1 #uses one of each stride time steps

###########################################################################
# Functions                                                               #
###########################################################################

def get_rock(density_i, rho_eff, surface_i):
    '''
    Effective density of the rock, extended upwards from the last node inside the rock of each column: the values
    around the surface mix rock and air.
    '''
    rock = rho_eff.copy()
    for j in range(Nx):
        inside = np.where((z <= surface_i[j] - dz) & (density_i[:, j] >= 1200.0))[0]
        rock[inside[-1]+1:, j] = rock[inside[-1], j]
    return rock

def get_columns(rock, surface_i, z_water=None):
    '''Density of the cells: rock below the surface inside the cells, water up to z_water (if any) and air.'''
    fraction_rock = np.clip((surface_i[None, :] - z_bottom[:, None])/dz, 0.0, 1.0)
    if(z_water is None):
        fraction_water = np.zeros_like(fraction_rock)
    else:
        fraction_water = np.clip((np.minimum(z_top[:, None], z_water) - np.maximum(z_bottom[:, None], surface_i[None, :]))/dz, 0.0, 1.0)
    return fraction_rock*rock + fraction_water*rho_water + (1.0 - fraction_rock - fraction_water)*rho_air

def sink_columns(rock, surface_i, z_water):
    '''
    Airy isostasy of the water load: each column below z_water sinks by rho_water*depth/(rho_m - rho_water), where
    rho_m is the density of the rock that leaves the column at the base of the model, so its mass does not change.
    '''
    depth = np.maximum(z_water - surface_i, 0.0)
    shift = rho_water*depth/(rock[0] - rho_water)
    sunk = np.empty_like(rock)
    for j in range(Nx):
        sunk[:, j] = np.interp(z + shift[j], z, rock[:, j])
    return sunk, surface_i - shift

###########################################################################
# Reading the mesh and the surface                                        #
###########################################################################

#tapioca files: data (time, z, x), x from 0 to Lx and z from -Lz to 0 [m]
with xr.open_dataset(f"{model_path}/density.nc", engine="h5netcdf") as mesh:
    x = mesh.x.values.astype(float)
    z = mesh.z.values.astype(float)
    time_mesh = mesh.time.values.astype(float)
    step_mesh = mesh.step.values

#surface.nc has its own x, with more points than the mesh
with xr.open_dataset(f"{model_path}/surface.nc", engine="h5netcdf") as surface_dataset:
    surface_x = surface_dataset.x.values.astype(float)
    surface_step = surface_dataset.step.values
    surface_all = surface_dataset.surface.values.astype(float)

if(not np.array_equal(surface_step, step_mesh)):
    raise ValueError("The steps of surface.nc and density.nc are different.")

Nx = x.size
Nz = z.size
Lx = x[-1] - x[0]
Lz = z[-1] - z[0]
dx = Lx/(Nx-1)
dz = Lz/(Nz-1)
z_bottom = z - dz/2
z_top = z + dz/2

print(
    "nx:", Nx, "\n",
    "nz:", Nz, "\n",
    "Lx:", Lx, "\n",
    "Lz:", Lz
)

idx_times = np.arange(0, time_mesh.size, stride)
time = time_mesh[idx_times]
step = step_mesh[idx_times]
n_times = idx_times.size

#Surface on the mesh x
surface = np.array([np.interp(x, surface_x, surface_all[i]) for i in idx_times])

###########################################################################
# Far-field columns and sea level                                         #
###########################################################################

far_columns = (x >= x_far_begin) & (x <= x_far_end)
if(not far_columns.any()):
    raise ValueError("No columns between x_far_begin and x_far_end.")

#Sea level: the mean surface of the far-field columns in the first time step
cond_sea = (surface_x >= x_far_begin) & (surface_x <= x_far_end)
sea_level = float(surface_all[0, cond_sea].mean())

print(f"Far-field columns between {x_far_begin/1000} and {x_far_end/1000} km")
print(f"Sea level: {sea_level} m. Water: {water}")

###########################################################################
# Prisms                                                                  #
###########################################################################

half_y = Lx*half_y_factor

X, Z = np.meshgrid(x, z)

#west, east, south, north, bottom, top
prisms = np.column_stack([
                        X.ravel() - dx/2,
                        X.ravel() + dx/2,
                        np.full(X.size, -half_y),
                        np.full(X.size, half_y),
                        Z.ravel() - dz/2,
                        Z.ravel() + dz/2,
                        ])

print(f"Number of prisms: {prisms.shape[0]}, observation points: {Nx}, time steps: {n_times}")

###########################################################################
# Forward modelling                                                       #
###########################################################################

#Compiling the numba function before the fork, so each process does not need to compile it again
hm.prism_gravity((0.0, 0.0, 1.0), [-1.0, 1.0, -1.0, 1.0, -2.0, -1.0], 1.0, field="g_z", parallel=False)

#We need shared pymp arrays to store data.
gz = pymp.shared.array((n_times, Nx), dtype=np.float64)
surface_sunk = pymp.shared.array((n_times, Nx), dtype=np.float64) #surface with the water load
z_obs = pymp.shared.array((n_times, Nx), dtype=np.float64)
gz[:] = np.nan

sys.stdout.flush() #otherwise, when the output goes to a file, each forked process writes the buffered text again

with pymp.Parallel() as p:
    #Opening the files inside the parallel region: each process has its own handles
    density = xr.open_dataset("density.nc", cache=False, engine="h5netcdf")["density"]
    temperature = xr.open_dataset("temperature.nc", cache=False, engine="h5netcdf")["temperature"]

    for it in p.range(0, n_times):
        i = idx_times[it]
        density_i = density.isel(time=i).values.astype(float)
        temperature_i = temperature.isel(time=i).values.astype(float)

        rho_eff = density_i*(1 - alpha*temperature_i) #effective density
        rock = get_rock(density_i, rho_eff, surface[it])

        columns_dry = get_columns(rock, surface[it])
        reference = np.median(columns_dry[:, far_columns], axis=1, keepdims=True) #reference density in each depth

        if(water):
            rock_sunk, surface_sunk[it] = sink_columns(rock, surface[it], sea_level)
            columns = get_columns(rock_sunk, surface_sunk[it], sea_level)
            z_obs[it] = np.maximum(surface_sunk[it], sea_level) + obs_height
        else:
            surface_sunk[it] = surface[it]
            columns = columns_dry
            z_obs[it] = surface[it] + obs_height

        drho = columns - reference #density contrast

        #parallel=False: the time steps are already split between the pymp processes
        coordinates = (x, np.zeros_like(x), z_obs[it])
        gz[it] = hm.prism_gravity(coordinates, prisms, drho.ravel(), field="g_z", parallel=False)

        p.print(f"Thread {p.thread_num}: time {time[it]:.2f} Myr, step {step[it]} ({it+1}/{n_times})")

if(np.isnan(gz).any()):
    raise ValueError("g_z has NaN values: a time step failed.")

###########################################################################
# Writing the dataset                                                     #
###########################################################################

ds = xr.Dataset(
    {
        "gz": (["time", "x"], np.array(gz), {"long_name": "downward component of gravity", "units": "mGal"}),
        "step": (["time"], step),
        "surface": (["time", "x"], surface, {"long_name": "surface of the model on the mesh x", "units": "m"}),
        "surface_sunk": (["time", "x"], np.array(surface_sunk), {"long_name": "surface with the water load", "units": "m"}),
        "z_obs": (["time", "x"], np.array(z_obs), {"long_name": "height of the observation points", "units": "m"}),
    },
    coords={
        "time": ("time", time, {"units": "Myr"}),
        "x": ("x", x, {"units": "m"}),
    },
    attrs={
        "model_name": model_name,
        "sea_level": sea_level,
        "x_far_begin": x_far_begin,
        "x_far_end": x_far_end,
        "water": int(water),
        "rho_water": rho_water,
        "obs_height": obs_height,
        "alpha": alpha,
        "half_y": half_y,
    },
)

ds.to_netcdf(f"{model_path}/_gz_evolution.nc", engine="h5netcdf")

print(f"gz min/max: {np.min(gz):.2f}/{np.max(gz):.2f} mGal")
print("_gz_evolution.nc written")
