# Description
#------------
# Main CMakeLists file for building goat. Here, only some main general
# options are specified, e.g. compiler options etc (depending on the build),
# all dependencies etc are specified through cmakelists.txt files per 
# directory. External libraries are added as external projects in order
# to expose the required .mod files and enable automatic recompilation 
# if their sources should change. 

# Note: it is assumed that this cmake file is ran in single configuration
# mode. Configurations may differ of course, but multi-config builds
# are not supported (mainly due to lack of support for .mod files of
# fortran modules, which would result in rebuilds all the time)

# Preliminaries
#--------------
# Check that CMAKE_BUILD_TYPE has been set
if(NOT CMAKE_BUILD_TYPE)
    message([FATAL_ERROR] "CMAKE_BUILD_TYPE not set")
endif()

# CMake compilation instructions (test) 
cmake_minimum_required(VERSION 3.23)

# Project declaration
project(goat)

# Languages
enable_language(Fortran C CXX)

# Compiler standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Shortcut for build type
set(cfg ${CMAKE_BUILD_TYPE})

# Aesthetics
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    add_compile_options(-fdiagnostics-color=always)
endif()

# Include cmake helper functions
#-------------------------------
include(cmake/SetupFortranTarget.cmake)

# Include required libraries
#---------------------------
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)


find_library(UMFPACK_LIB umfpack REQUIRED)
find_library(AMD_LIB amd REQUIRED)
find_library(CHOLMOD_LIB cholmod REQUIRED)
find_library(COLAMD_LIB colamd REQUIRED)
find_library(SUITESPARSECONFIG_LIB suitesparseconfig REQUIRED)
find_library(CXSPARSE_LIB cxsparse REQUIRED)

find_path(SUITESPARSE_INCLUDE_DIR
    NAMES cs.h
    PATH_SUFFIXES suitesparse
    REQUIRED
)


# Program dependencies
#---------------------
# tcsh shell 
find_program(TCSH_EXECUTABLE tcsh REQUIRED
    DOC "tcsh required for legacy submodule environment setup")

# Compiler options
#-----------------
# Add compiler option "library" (to be used by all subsequent libraries)
#add_library(psiduqc_flags INTERFACE)

# General compilation options for gfortran
#target_compile_options(psiduqc_flags INTERFACE 
#    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wall -Wextra -fimplicit-none -fdiagnostics-color=always>)

# Gfortran options for different configs
#target_compile_options(psiduqc_flags INTERFACE
#    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:
#        $<$<CONFIG:Debug>:-O0 -g -fcheck=all -fbacktrace>
#        $<$<CONFIG:Release>:-O3 -g>
#    >)

# Include source
#---------------
# Add subdirectory
add_subdirectory(src)
add_subdirectory(Runfiles)
