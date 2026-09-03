# Interpolant.F90
add_library(Interpolant STATIC
    Interpolant.F90)

setup_fortran_target(Interpolant)

target_link_libraries(Interpolant PUBLIC 
    StructuredInterpolant2D)

# Interpolant1D.F90
add_library(Interpolant1D STATIC
    Interpolant1D.F90)

setup_fortran_target(Interpolant1D)

target_link_libraries(Interpolant1D PUBLIC 
    mod_precision
    mod_constants
    mod_errorhandler
    Interpolant2D_auxiliaries)

# Interpolant2D_auxiliaries
add_library(Interpolant2D_auxiliaries STATIC
    Interpolant2D_auxiliaries.F90)

setup_fortran_target(Interpolant2D_auxiliaries)

target_link_libraries(Interpolant2D_auxiliaries PUBLIC 
    mod_precision
    mod_errorhandler
    mod_sort
    mod_triangulation)

# Interpolant2D
add_library(Interpolant2D STATIC
    Interpolant2D.F90)

setup_fortran_target(Interpolant2D)

target_link_libraries(Interpolant2D PUBLIC 
    mod_precision
    mod_sparseinterface
    mod_plotter
    mod_structured2Dgridding)

# StructuredInterpolant2D
add_library(StructuredInterpolant2D STATIC
    StructuredInterpolant2D.F90)

setup_fortran_target(StructuredInterpolant2D)

target_link_libraries(StructuredInterpolant2D PUBLIC 
    Interpolant2D_auxiliaries
    Interpolant2D
    Interpolant1D
    mod_constants
    mod_sparseinterface
    mod_linearsolverinterface)

# UnstructuredInterpolation2D
add_library(UnstructuredInterpolation2D STATIC
    UnstructuredInterpolation2D.F90)

setup_fortran_target(UnstructuredInterpolation2D)

target_link_libraries(UnstructuredInterpolation2D PUBLIC 
    mod_triangulation  
    mod_utility
    Interpolant2D
    Interpolant2D_auxiliaries)

# Library of all routines in Interpolation
add_library(Interpolation INTERFACE)
target_link_libraries(Interpolation INTERFACE 
    Interpolant
    Interpolant1D
    Interpolant2D_auxiliaries
    Interpolant2D
    StructuredInterpolant2D
    UnstructuredInterpolation2D)