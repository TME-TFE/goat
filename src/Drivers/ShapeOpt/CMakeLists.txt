# ShapeOptDriver
add_library(ShapeOptDriver STATIC ShapeOptDriver.F90)
setup_fortran_target(ShapeOptDriver)
target_link_libraries(ShapeOptDriver PUBLIC 
    somod_userinput
    somod_optimizationengine
    goatmod_types
    ggmod_gridgenerator)

