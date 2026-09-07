# somod_constraints
add_library(somod_constraints STATIC
    somod_constraints.F90)
setup_fortran_target(somod_constraints)
target_link_libraries(somod_constraints PUBLIC 
    optmod_constraints
    optmod_optimizationengine
    somod_designvariables
    somod_userinput
    gdmod_constraints
    gdmod_optimizationengine
    PolygonLevelsetFunction2D)

# somod_costfunction
add_library(somod_costfunction STATIC
    somod_costfunction.F90)
setup_fortran_target(somod_costfunction)
target_link_libraries(somod_costfunction PUBLIC 
    optmod_costfunction
    gdmod_optimizationengine
    somod_userinput
    somod_designvariables
    PolygonLevelsetFunction2D
    mod_linearsolverinterface
    optmod_hessapprox)

# somod_designvariables
add_library(somod_designvariables STATIC
    somod_designvariables.F90)
setup_fortran_target(somod_designvariables)
target_link_libraries(somod_designvariables PUBLIC 
    optmod_designvariables
    gdmod_optimizationengine)

# somod_optimizationengine
add_library(somod_optimizationengine STATIC
    somod_optimizationengine.F90)
setup_fortran_target(somod_optimizationengine)
target_link_libraries(somod_optimizationengine PUBLIC 
    mod_precision
    optmod_optimizationengine
    somod_designvariables
    somod_costfunction
    somod_constraints
    gdmod_optimizationengine
    gdmod_designvariables
    gdmod_costfunction
    gdmod_constraints
    gdmod_types
    gdmod_plots
    mod_specialchars
    mod_plotter
    GDinitialize)

# somod_userinput
add_library(somod_userinput STATIC
    somod_userinput.F90)
setup_fortran_target(somod_userinput)
target_link_libraries(somod_userinput PUBLIC 
    mod_precision
    mod_readwrite
    mod_inputfileparser
    mod_constants)

# sosmod_costfunction - currently not used since requires SOLPS
#add_library(sosmod_costfunction STATIC
#    sosmod_costfunction.F90)
#setup_fortran_target(sosmod_costfunction)
#target_link_libraries(sosmod_costfunction PUBLIC 
#    libsolps_b25
#    somod_costfunction 
#    gdmod_optimizationengine
#    mod_precision)