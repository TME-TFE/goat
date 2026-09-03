# gdmod_constraints
add_library(gdmod_constraints STATIC
    gdmod_constraints.F90)
setup_fortran_target(gdmod_constraints)
target_link_libraries(gdmod_constraints PUBLIC 
    mod_precision
    mod_utility
    mod_sparseinterface
    optmod_constraints
    gdmod_types
    gdmod_designvariables
    gdmod_utility_optimization
    gdmod_plots
    gdmod_userinput
    mod_definitions
    PolygonLevelsetFunction2D)

# gdmod_costfunction
add_library(gdmod_costfunction STATIC
    gdmod_costfunction.F90)
setup_fortran_target(gdmod_costfunction)
target_link_libraries(gdmod_costfunction PUBLIC 
    mod_precision
    mod_utility
    mod_sparseinterface
    gdmod_types 
    gdmod_designvariables
    gdmod_utility_optimization
    optmod_costfunction
    gdmod_plots
    DistributionFunction)

# gdmod_designvariables
add_library(gdmod_designvariables STATIC
    gdmod_designvariables.F90)
setup_fortran_target(gdmod_designvariables)
target_link_libraries(gdmod_designvariables PUBLIC 
    gdmod_types
    gdmod_userinput 
    optmod_designvariables)

# gdmod_optimizationengine
add_library(gdmod_optimizationengine STATIC
    gdmod_optimizationengine.F90)
setup_fortran_target(gdmod_optimizationengine)
target_link_libraries(gdmod_optimizationengine PUBLIC 
    mod_precision
    optmod_optimizationengine
    gdmod_designvariables
    gdmod_costfunction
    gdmod_constraints
    gdmod_types
    gdmod_userinput
    gdmod_plots
    mod_plotter
    mod_specialchars)

# gdmod_plots
add_library(gdmod_plots STATIC
    gdmod_plots.F90)
setup_fortran_target(gdmod_plots)
target_link_libraries(gdmod_plots PUBLIC 
    mod_plotter
    gdmod_types
    gdmod_userinput )

# gdmod_types
add_library(gdmod_types STATIC
    gdmod_types.F90)
setup_fortran_target(gdmod_types)
target_link_libraries(gdmod_types PUBLIC 
    mod_precision
    goatmod_types)

# gdmod_userinput
add_library(gdmod_userinput STATIC
    gdmod_userinput.F90)
setup_fortran_target(gdmod_userinput)
target_link_libraries(gdmod_userinput PUBLIC 
    gdmod_types
    goatmod_userinput
    mod_inputfileparser)

# gdmod_utility_optimization
add_library(gdmod_utility_optimization STATIC
    gdmod_utility_optimization.F90)
setup_fortran_target(gdmod_utility_optimization)
target_link_libraries(gdmod_utility_optimization PUBLIC 
    gdmod_types
    gdmod_userinput 
    mod_definitions)