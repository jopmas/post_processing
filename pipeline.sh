#!/usr/bin/env bash

########################################################################
# Post-processing pipeline of a Mandyoc scenario from its zip file:    #
#                                                                      #
#  1. extracts the zip in the folder of the zip                        #
#  2. organizes the outputs in sub-folders (0_organize_outputs.sh of   #
#     tapioca; nothing to do if they are already organized)            #
#  3. creates the .nc files with the Julia scripts of tapioca          #
#     (meshes2netcdf.jl, and particles2netcdf.jl with the last step as #
#     reference step if there are particles)                           #
#  4. removes the extracted folders of the outputs (the zip is kept)   #
#  5. snapshots of the scenario (plot_snapshots.py)                    #
#  6. frames and video of the lithology (frames_generator.py)          #
#                                                                      #
# Usage, with the tapioca environment active:                          #
#   bash ~/opt/post_processing/pipeline.sh /path/to/scenario.zip [threads]
#                                                                      #
# threads: Julia threads (default: auto). The folder of tapioca is     #
# ~/opt/tapioca, or the variable TAPIOCA_DIR.                          #
########################################################################

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash pipeline.sh /path/to/scenario.zip [threads]"
    exit 1
fi

ZIP_FILE="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
THREADS="${2:-auto}"
TAPIOCA_DIR="${TAPIOCA_DIR:-$HOME/opt/tapioca}"
POST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # this repository
SCENARIO_DIR="$(dirname "$ZIP_FILE")"

# Folders of the Mandyoc outputs created by the extraction and by 0_organize_outputs.sh (removed in step 4)
OUTPUT_FOLDERS=(time density heat pressure surface strain_rate strain temperature viscosity velocity
                thermal_diffusivity X_depletion Phi dPhi lithos steps)

# Files needed by the snapshots and the video (checked before removing the extracted folders)
NEEDED_FILES=(temperature.nc density.nc strain.nc lithology.nc)

if [ ! -f "$ZIP_FILE" ]; then
    echo "File not found: $ZIP_FILE"
    exit 1
fi

cd "$SCENARIO_DIR"
echo "Scenario: $SCENARIO_DIR"

echo "========== 1. Extracting $(basename "$ZIP_FILE")"
# -n: never overwrite existing files. unzip returns 1 for warnings (e.g. skipped files)
unzip -q -n "$ZIP_FILE" || [ $? -eq 1 ]

echo "========== 2. Organizing the outputs"
bash "$TAPIOCA_DIR/src/postrunning/0_organize_outputs.sh"

echo "========== 3. Creating the .nc files"
julia -t "$THREADS" "$TAPIOCA_DIR/src/postrunning/meshes2netcdf.jl" "$SCENARIO_DIR"

if [ -d steps ]; then
    # last step of the model: the tracking scripts select the particles there
    last_step=$(find time -maxdepth 1 -name "time_*.txt" | sed -E 's|.*time_([0-9]+)\.txt|\1|' | grep -E '^[0-9]+$' | sort -n | tail -1)
    echo "Particles with reference step $last_step"
    julia -t "$THREADS" "$TAPIOCA_DIR/src/postrunning/particles2netcdf.jl" "$SCENARIO_DIR" "$last_step"
else
    echo "No steps/ folder: particles_trajectories.nc is not created"
fi

for file in "${NEEDED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "$file was not created (check export_lithology = True in param.txt for lithology.nc)."
        echo "Stopping before removing the extracted folders."
        exit 1
    fi
done

echo "========== 4. Removing the extracted folders"
for folder in "${OUTPUT_FOLDERS[@]}"; do
    if [ -d "$folder" ]; then
        rm -rf "./$folder"
        echo "Removed $folder/"
    fi
done

echo "========== 5. Snapshots"
python "$POST_DIR/plot_snapshots.py"

echo "========== 6. Frames and video"
python "$POST_DIR/frames_generator.py"

echo "========== Done: $SCENARIO_DIR"
