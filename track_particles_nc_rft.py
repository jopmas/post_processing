########################################################################
# BEFORE RUNNING THIS SCRIPT:                                          #
#  - CREATE THE .nc FILES WITH tapioca's meshes2netcdf.jl              #
#    (temperature.nc, pressure.nc, ...)                                #
#  - CREATE particles_trajectories.nc WITH tapioca's                   #
#    particles2netcdf.jl                                               #
#  - RUN THIS SCRIPT IN THE SAME FOLDER WHERE THE .nc FILES ARE        #
#                                                                      #
# Rifting models. Same as track_particles.py, but reads the particles  #
# from particles_trajectories.nc instead of the steps/step*.txt files  #
#                                                                      #
# particles_trajectories.nc only has the particles that exist in its   #
# reference time step: create it with the final step of the tracking   #
# (step of instant_to_take below) as reference to select all           #
# particles there.                                                     #
########################################################################

import numpy as np
import os
import xarray as xr
import subprocess
import multiprocessing #needed to run pymp in macOS
multiprocessing.set_start_method('fork') #needed to run pymp in macOS arm architectures
import pymp

model_path = os.getcwd() # Get local file
model_name = model_path.split('/')[-1]

#reading datasets created by tapioca
Tdataset = xr.open_dataset(f"{model_path}/temperature.nc") #temperature (time, z, x)
Pdataset = xr.open_dataset(f"{model_path}/pressure.nc") #pressure (time, z, x)
Sdataset = xr.open_dataset(f"{model_path}/particles_trajectories.nc") #particles x, z, lithology (time, id)

#getting the dimensions (tapioca: x from 0 to Lx and z from -Lz to 0 [m])
Nx = int(Tdataset.sizes["x"])
Nz = int(Tdataset.sizes["z"])
Lx = float(Tdataset.x.values[-1] - Tdataset.x.values[0])
Lz = float(Tdataset.z.values[-1] - Tdataset.z.values[0])

print(
    "nx:", Nx, "\n",
    "nz:", Nz, "\n",
    "Lx:", Lx, "\n",
    "Lz:", Lz
)

xi = np.linspace(0,Lx,Nx)
zi = np.linspace(-Lz,0,Nz)
xx,zz = np.meshgrid(xi,zi)

dx = Lx/(Nx-1)
dz = Lz/(Nz-1)

if('He' in model_name):
    take_asthenosphere = True
else:
    take_asthenosphere = False

#select the extension rate:
if('W1' in model_name or 'S1' in model_name):
    v=1 #cm/yr #For groups 1,2 and 3
if('W2' in model_name or 'S2' in model_name):
    v=2 #cm/yr #For groups 4,5 and 6
if('W3' in model_name or 'S3' in model_name):
    v=3 #cm/yr #For groups 7,8 and 9

if(v==1):
    instant_to_take = 40 #Myr
if(v==2):
    instant_to_take = 40#20 #Myr
if(v==3):
    instant_to_take = 40#8.55 #Myr

if(Tdataset.time.values[-1] < instant_to_take):
    print(f"WARNING: the model only reaches {Tdataset.time.values[-1]:.2f} Myr (< {instant_to_take} Myr). The last step will be used.")

idx = (np.abs(Tdataset.time.values - instant_to_take)).argmin()
step_final = Tdataset.step.values[idx]

print(f"Track particles until {instant_to_take} Myr, idx: {idx}, time idx: {Tdataset.time.values[idx]} Myr")

print("Final step: ", step_final)

#index of each step of Tdataset inside particles_trajectories.nc
steps_nc = Sdataset.step.values
ids_nc = Sdataset["id"].values #particles ids, sorted
missing_steps = np.setdiff1d(Tdataset.step.values[:idx+1], steps_nc)
if(missing_steps.size > 0):
    raise ValueError(f"Steps not found in particles_trajectories.nc: {missing_steps}")

def steps_nc_index(step):
    return int(np.where(steps_nc == step)[0][0])

def particles_positions(k):
    '''x and z [m] of all particles in the time index k of particles_trajectories.nc, NaN where the particle does not exist (-999 in the file)'''
    x_k = Sdataset["x"].isel(time=k).values.astype(float)
    z_k = Sdataset["z"].isel(time=k).values.astype(float)
    missing = x_k == -999.0
    x_k[missing] = np.nan
    z_k[missing] = np.nan
    return x_k, z_k

reference_step = int(Sdataset.attrs["reference_timestep"])
if(reference_step != step_final):
    print(f"WARNING: particles_trajectories.nc was created with reference step {reference_step}, not with the final step {step_final}.")
    print(f"         Particles that do not exist in step {reference_step} cannot be selected.")

#Reading the data of final step
x_final, z_final = particles_positions(steps_nc_index(step_final))
exists = np.isfinite(x_final)
x = x_final[exists]
z = z_final[exists]
id_vec = ids_nc[exists]
layer_vec = Sdataset["lithology"].isel(time=steps_nc_index(step_final)).values[exists] #layer code of each particle (constant in time)


###########################################################################
#Selecting the particles to track                                         #
# All particles above search_thickness and between x_begin and x_end [km] #
# and does not belong to air                                              #
###########################################################################

h_air = 40.0e3

#Determining the thickness to search for lithospheric mantle particles
search_thickness = 12.0e3

x_begin = 600.0e3
x_end = 1600.0e3

asthenosphere_code = 0
mantle_lithosphere1_code = 1
seed_code = 2
mantle_lithosphere2_code = 3
lower_crust_code = 4
upper_crust_code = 5
air_code = 6

print(f"Selecting particles between {x_begin/1000} and {x_end/1000} km, above {search_thickness/1000} km")

if(take_asthenosphere):
    cond = (z>-(h_air + search_thickness)) & (x>=x_begin) & (x<=x_end) & (layer_vec<upper_crust_code) & (layer_vec>=asthenosphere_code) #to take particles from asthenosphere and lithosphere.
else:
    cond = (z>-(h_air + search_thickness)) & (x>=x_begin) & (x<=x_end) & (layer_vec<upper_crust_code) & (layer_vec>asthenosphere_code) #to take particles from lithosphere

part_selec = id_vec[cond]
layers_selec = layer_vec[cond] #get the layer numbers of the selected particles
n_tracked = np.size(part_selec)

print(f"x size: {np.size(x)}")
print(f"Number of tracked particles (n): {np.size(part_selec)}, where:")

n_tracked_ast = layers_selec[layers_selec==asthenosphere_code].size
n_tracked_mantle_lithosphere1 = layers_selec[layers_selec==mantle_lithosphere1_code].size
n_tracked_mantle_lithosphere2 = layers_selec[layers_selec==mantle_lithosphere2_code].size
n_tracked_lower_crust = layers_selec[layers_selec==lower_crust_code].size
n_tracked_upper_crust = layers_selec[layers_selec==upper_crust_code].size

print(f"Asthenosphere: {n_tracked_ast}")
print(f"Mantle Lithosphere: {n_tracked_mantle_lithosphere1+n_tracked_mantle_lithosphere2}")
print(f"Upper Crust: {n_tracked_upper_crust}")
print(f"Lower Crust: {n_tracked_lower_crust}\n")

if(n_tracked == 0):
    raise ValueError("No particles selected. Check search_thickness, x_begin, x_end and instant_to_take.")

#Reading the whole history of the selected particles at once
start = 0
end = idx
print(f'Start idx: {start}, End idx: {end}, time of end: {Tdataset.time.values[end]}')

#Reading one step at a time and indexing with numpy is much faster than .sel(id=...) on the file
part_idx = np.searchsorted(ids_nc, part_selec)
xtrack = np.full((end-start+1, n_tracked), np.nan) #(time, particle)
ztrack = np.full((end-start+1, n_tracked), np.nan)
for it, s in enumerate(Tdataset.step.values[start:end+1]):
    x_k, z_k = particles_positions(steps_nc_index(s))
    xtrack[it] = x_k[part_idx]
    ztrack[it] = z_k[part_idx]

n_missing = np.count_nonzero(~np.isfinite(xtrack))
if(n_missing > 0):
    print(f"WARNING: {n_missing} (time, particle) positions do not exist in particles_trajectories.nc. P and T will be NaN there.")

#We need shared pymp arrays to store data.
ptrack = pymp.shared.array(xtrack.shape, dtype=np.float64)
ttrack = pymp.shared.array(xtrack.shape, dtype=np.float64)
ptrack[:] = np.nan
ttrack[:] = np.nan

#Interpolating P and T at the particle positions
with pymp.Parallel() as p:
    for i in p.range(start, end+1):
        it = i - start
        temperature = Tdataset.temperature[i].values
        pressure = Pdataset.pressure[i].values

        vecx_track = xtrack[it]
        vecz_track = ztrack[it]
        ok = np.isfinite(vecx_track)

        jm = np.clip(((vecx_track[ok]-0)//dx).astype(int), 0, Nx-2)
        im = np.clip(((vecz_track[ok]+Lz)//dz).astype(int), 0, Nz-2)

        xw = (vecx_track[ok] - xx[im,jm])/dx
        zw = (vecz_track[ok] - zz[im,jm])/dz

        present_temperature =   temperature[im,jm]*    (1.-xw)*(1.-zw)+\
                                temperature[im,jm+1]*  (xw   )*(1.-zw)+\
                                temperature[im+1,jm]*  (1.-xw)*(   zw)+\
                                temperature[im+1,jm+1]*(   xw)*(   zw)

        present_pressure =  pressure[im,jm]*    (1.-xw)*(1.-zw)+\
                            pressure[im,jm+1]*  (xw   )*(1.-zw)+\
                            pressure[im+1,jm]*  (1.-xw)*(   zw)+\
                            pressure[im+1,jm+1]*(   xw)*(   zw)

        ptrack[it, ok] = present_pressure / 1.0e6
        ttrack[it, ok] = present_temperature

        np.savetxt(f"track_xzPT_step_{str(Tdataset.step.values[i]).zfill(6)}.txt", np.c_[vecx_track, vecz_track, ptrack[it], ttrack[it]], fmt="%.5f")

all_time = Tdataset.time.values[start:end+1]
all_step = Tdataset.step.values[start:end+1]

#Flattening (time, particle) as in track_particles.py: time ascending, then reversed
all_vecx_track = xtrack.ravel()
all_vecz_track = ztrack.ravel()
all_present_pressure = ptrack.ravel()
all_present_temperature = ttrack.ravel()

print(f"len of:\n all_vecx_track {len(all_vecx_track)}\n all_vecz_track {len(all_vecz_track)}\n all_present_pressure {len(all_present_pressure)}\n all_present_temperature {len(all_present_temperature)}\n all_time {len(all_time)}\n all_step {len(all_step)}")
# Creatiing the xarray dataset with the tracked particles
ds = xr.Dataset(
    {
        "xtrack": (["index"], all_vecx_track[::-1]),
        "ztrack": (["index"], all_vecz_track[::-1]),
        "ptrack": (["index"], all_present_pressure[::-1]),
        "ttrack": (["index"], all_present_temperature[::-1]),
        "step": all_step[::-1],
        "time": all_time[::-1],
        "ntracked": int(n_tracked),
        "particles_layers": layers_selec[::-1]
    },
    coords={
        "index": np.arange(len(all_vecx_track))
    }
)

ds.to_netcdf("_track_xzPT_all_steps.nc")


print(f"n_tracked x len(all_time) = {n_tracked}*{len(all_time)} = {n_tracked*len(all_time)}")

# Compressing the files
print("Zipping files")
#-m moves the files into the zip: each track*.txt is deleted only after it is in the zip
subprocess.run(f"zip -m {model_name}.zip track*.txt", shell=True, check=True, capture_output=True, text=True)
print("Files zipped and track*.txt removed")
