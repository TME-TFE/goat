# optmod_constraints
add_library(optmod_constraints STATIC 
    optmod_constraints.F90)
setup_fortran_target(optmod_constraints)
target_link_libraries(optmod_constraints PUBLIC 
    mod_precision)

# optmod_costfunction
add_library(optmod_costfunction STATIC 
    optmod_costfunction.F90)
setup_fortran_target(optmod_costfunction)
target_link_libraries(optmod_costfunction PUBLIC 
    mod_precision)

# optmod_designvariables
add_library(optmod_designvariables STATIC 
    optmod_designvariables.F90)
setup_fortran_target(optmod_designvariables)
target_link_libraries(optmod_designvariables PUBLIC 
    mod_precision)

# optmod_hessapprox
add_library(optmod_hessapprox STATIC 
    optmod_hessapprox.F90)
setup_fortran_target(optmod_hessapprox)
target_link_libraries(optmod_hessapprox PUBLIC 
    mod_precision
    mod_constants
    mod_sparseinterface
    mod_linearsolverinterface)

# optmod_monitor
add_library(optmod_monitor STATIC 
    optmod_monitor.F90)
setup_fortran_target(optmod_monitor)
target_link_libraries(optmod_monitor PUBLIC 
    mod_precision)

# optmod_numerics
add_library(optmod_numerics STATIC 
    optmod_numerics.F90)
setup_fortran_target(optmod_numerics)
target_link_libraries(optmod_numerics PUBLIC 
    mod_precision
    mod_readwrite
    mod_inputfileparser)

# optmod_optimizationengine
add_library(optmod_optimizationengine STATIC 
    optmod_optimizationengine.F90)
setup_fortran_target(optmod_optimizationengine)
target_link_libraries(optmod_optimizationengine PUBLIC 
    mod_precision
    mod_constants
    mod_utility
    mod_sparseinterface
    mod_diagnostics
    mod_linearsolverinterface
    mod_plotter
    mod_errorhandler
    optmod_designvariables
    optmod_costfunction
    optmod_constraints
    optmod_numerics
    optmod_monitor)
