# ComputeGradient2DStructured
add_library(ComputeGradient2DStructured STATIC 
    ComputeGradient2DStructured.F90)
setup_fortran_target(ComputeGradient2DStructured)
target_link_libraries(ComputeGradient2DStructured PUBLIC 
    mod_precision)

# DistributionFunction
add_library(DistributionFunction STATIC 
    DistributionFunction.F90)
setup_fortran_target(DistributionFunction)
target_link_libraries(DistributionFunction PUBLIC 
    mod_precision 
    Interpolant
    mod_linearsolverinterface
    mod_polygon
    PolygonLevelsetFunction2D
    mod_structured2Dgridding)

# mod_lagrangefunctions
add_library(mod_lagrangefunctions STATIC 
    mod_lagrangefunction.F90) # mind the s! should fix this...
setup_fortran_target(mod_lagrangefunctions)
target_link_libraries(mod_lagrangefunctions PUBLIC 
    mod_precision
    mod_errorhandler)

# PolygonLevelsetFunction
add_library(PolygonLevelsetFunction2D STATIC 
    PolygonLevelsetFunction2D.F90)
setup_fortran_target(PolygonLevelsetFunction2D)
target_link_libraries(PolygonLevelsetFunction2D PUBLIC 
    mod_precision
    mod_constants
    mod_plotter
    Interpolant
    mod_polygon
    mod_sparseinterface
    mod_structured2Dgridding)

# Library containing all libraries in Numerics
add_library(Numerics INTERFACE)
target_link_libraries(Numerics INTERFACE 
    ComputeGradient2DStructured
    DistributionFunction
    mod_lagrangefunctions
    PolygonLevelsetFunction2D)