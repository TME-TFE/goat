# A library per driver since each driver may depend on different 
# libraries. One library then to collect all driver libraries

# GADriver
add_library(GADriver STATIC GADriver.F90)
setup_fortran_target(GADriver)
target_link_libraries(GADriver PUBLIC 
    goatmod_types 
    goatmod_userinput
    gamod_driver)

# GAGDDriver
add_library(GAGDDriver STATIC GAGDDriver.F90)
setup_fortran_target(GAGDDriver)
target_link_libraries(GAGDDriver PUBLIC 
    goatmod_types 
    goatmod_userinput
    gamod_driver)

# GDDriver
add_library(GDDriver STATIC GDDriver.F90)
setup_fortran_target(GDDriver)
target_link_libraries(GDDriver PUBLIC 
    goatmod_types 
    goatmod_userinput
    ggmod_gridgeneration2D)

# GDinitialize
add_library(GDinitialize STATIC GDinitialize.F90)
setup_fortran_target(GDinitialize)
target_link_libraries(GDinitialize PUBLIC 
    goatmod_types
    goatmod_userinput
    gdmod_optimizationengine)

# GDtestdriver
add_library(GDtestdriver STATIC GDtestdriver.F90)
setup_fortran_target(GDtestdriver)
target_link_libraries(GDtestdriver PUBLIC 
    goatmod_types 
    goatmod_userinput
    gdmod_optimizationengine
    mod_constants
    GDinitialize)

# GGDriver
add_library(GGDriver STATIC GGDriver.F90)
setup_fortran_target(GGDriver)
target_link_libraries(GGDriver PUBLIC 
    goatmod_types 
    goatmod_userinput
    ggmod_topology2D
    ggmod_gridgeneration2D
    mod_contour2D
    mod_streamlinetracing2D
    mod_structured2Dgridding)

# GGGDDriver
add_library(GGGDDriver STATIC GGGDDriver.F90)
setup_fortran_target(GGGDDriver)
target_link_libraries(GGGDDriver PUBLIC 
    goatmod_types 
    goatmod_userinput
    ggmod_topology2D
    ggmod_gridgeneration2D
    mod_contour2D
    mod_streamlinetracing2D
    mod_structured2Dgridding)

# RunGridOptimization
add_library(RunGridOptimization STATIC RunGridOptimization.F90)
setup_fortran_target(RunGridOptimization)
target_link_libraries(RunGridOptimization PUBLIC 
    gdmod_types
    gdmod_userinput 
    gdmod_plots
    mod_plotter
    gdmod_optimizationengine)


# Library with all drivers
add_library(goat_drivers INTERFACE)
target_link_libraries(goat_drivers INTERFACE
    GADriver
    GAGDDriver
    GDDriver
    GDinitialize
    GDtestdriver
    GGDriver
    GGGDDriver
    RunGridOptimization)
