#!/bin/bash
# Simple install script for goat. Should be ran in the main 
# directory (to be checked if we can automate this...)

# Construct build directories for the different configurations
mkdir build-debug
mkdir build-release

# Execute cmake
cmake -S . -B build-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug --target goat --verbose
cmake --install build-debug --prefix ./executables/Debug  --verbose
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release 
cmake --install build-release --prefix ./executables/Release

# Construct symlinks to executables
ln -sfn ./executables/Debug/goat goat_debug.exe
ln -sfn ./executables/Release/goat goat.exe
