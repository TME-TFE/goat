# ggmod_gridgeneration2D
add_library(ggmod_gridgeneration2D STATIC 
    ggmod_gridgeneration2D.F90)
setup_fortran_target(ggmod_gridgeneration2D)
target_link_libraries(ggmod_gridgeneration2D PUBLIC 
    General
    Constants
    Auxiliary
    Numerics
    goatmod_types
    goatmod_userinput
    ggmod_vertexdistribution2D
    ggmod_topology2D)
target_link_libraries(ggmod_gridgeneration2D PRIVATE 
    goat_openmp)

# ggmod_gridgenerator
add_library(ggmod_gridgenerator STATIC 
    ggmod_gridgenerator.F90)
setup_fortran_target(ggmod_gridgenerator)
target_link_libraries(ggmod_gridgenerator PUBLIC 
    goatmod_types 
    goatmod_userinput
    ggmod_gridgeneration2D
    ggmod_topology2D
    mod_contour2D
    mod_streamlinetracing2D
    mod_errorhandler)

# ggmod_topology2D
add_library(ggmod_topology2D STATIC 
    ggmod_topology2D.F90)
setup_fortran_target(ggmod_topology2D)
target_link_libraries(ggmod_topology2D PUBLIC 
    General
    Constants
    Auxiliary
    Numerics
    goatmod_types
    goatmod_userinput)
target_link_libraries(ggmod_topology2D PRIVATE 
    goat_openmp)

# ggmod_vertexdistribution2D
add_library(ggmod_vertexdistribution2D STATIC 
    ggmod_vertexdistribution2D.F90)
setup_fortran_target(ggmod_vertexdistribution2D)
target_link_libraries(ggmod_vertexdistribution2D PUBLIC 
    mod_precision
    mod_errorhandler
    goatmod_types
    Interpolant1D
    StructuredInterpolant2D
    mod_lagrangefunctions
    mod_constants
    DistributionFunction)
