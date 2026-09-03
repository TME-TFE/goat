# goatmod_types
add_library(goatmod_types STATIC 
    goatmod_types.F90)
setup_fortran_target(goatmod_types)
target_link_libraries(goatmod_types PUBLIC 
    goatmod_userinput
    mod_precision
    mod_definitions
    mod_constants
    mod_polygon
    mod_inputfileparser
    mod_plotter
    mod_dynamicarrays
    mod_readwrite
    mod_sort
    Interpolant
    PolygonLevelsetFunction2D
    mod_std_formatspecs
    PolygonShapeFunction)

# goatmod_userinput
add_library(goatmod_userinput STATIC 
    goatmod_userinput.F90)
setup_fortran_target(goatmod_userinput)
target_link_libraries(goatmod_userinput PUBLIC 
    mod_precision
    mod_readwrite
    mod_inputfileparser
    mod_global_environment
    mod_constants)