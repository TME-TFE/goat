add_library(mod_contour2D STATIC
    mod_contour2D.F90)

setup_fortran_target(mod_contour2D)

target_link_libraries(mod_contour2D PUBLIC 
    mod_precision
    mod_errorhandler
    mod_dynamicarrays
    mod_structured2Dgridding
    mod_plotter
    mod_sort
    mod_constants
    Interpolant1D
    Interpolant2D_auxiliaries
    mod_specialchars
    mod_definitions
    mod_inputfileparser)

# Library of all routines in Contour
add_library(Contour INTERFACE)
target_link_libraries(Contour INTERFACE 
    mod_contour2D)