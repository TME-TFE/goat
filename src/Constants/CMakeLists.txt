# mod_constants
add_library(mod_constants STATIC
    mod_constants.F90)

setup_fortran_target(mod_constants)

target_link_libraries(mod_constants PUBLIC 
    mod_precision)

# mod_definitions
add_library(mod_definitions STATIC
    mod_definitions.F90)

setup_fortran_target(mod_definitions)

target_link_libraries(mod_definitions PUBLIC 
    mod_precision)

# mod_global_environment
add_library(mod_global_environment STATIC
    mod_global_environment.F90)

setup_fortran_target(mod_global_environment)

# mod_precision
add_library(mod_precision STATIC
    mod_precision.F90)

setup_fortran_target(mod_precision)

# mod_specialchars
add_library(mod_specialchars STATIC
    mod_specialchars.F90)

setup_fortran_target(mod_specialchars)

# mod_std_formatspecs
add_library(mod_std_formatspecs STATIC
    mod_std_formatspecs.F90)

setup_fortran_target(mod_std_formatspecs)

# Library with all libraries in Constants
add_library(Constants INTERFACE)
target_link_libraries(Constants INTERFACE 
    mod_constants
    mod_definitions 
    mod_global_environment
    mod_precision
    mod_specialchars
    mod_std_formatspecs)