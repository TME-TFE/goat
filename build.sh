#!/bin/bash
# Simple install script for goat. Should be ran in the main 
# directory (to be checked if we can automate this...)

# Construct build directories for the different configurations
mkdir build-debug
mkdir build-release

# Add top directory to the path variable, if not there already
if [[ ! -v GOATTOP ]]; then
    echo "GOATTOP variable is not defined, not adding to path. Do source setup.csh first"
elif [[ ":$PATH:" != *":$GOATTOP:"* ]]; then
    export PATH="$GOATTOP:$PATH"
fi

# Execute cmake
cmake -S . -B build-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --target goat --verbose
cmake --build build-debug --target shapeopt --verbose
cmake --install build-debug --prefix ./executables/Debug  --verbose
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target goat --verbose
cmake --build build-release --target shapeopt --verbose
cmake --install build-release --prefix ./executables/Release

# Construct symlinks to executables
ln -sfn ./executables/Debug/goat goat_debug.exe
ln -sfn ./executables/Release/goat goat.exe
ln -sfn ./executables/Debug/shapeopt shapeopt_debug.exe
ln -sfn ./executables/Release/shapeopt shapeopt.exe

