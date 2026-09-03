# Function to setup fortran targets properly. Currently, it handles the
# following things:
#   1)  .mod file inclusion: separate module directory per target and per 
#       config to handle multi-config and multi-target builds
#   2)  include directories also automated
#   3)  TO DO: compiler flag propagation
# Function to handle .mod files consistently. It places the .mod files
# of each target into a separate folder that is configuration dependent.
# As such, it should work for multi-target, multi-config builds. The
# generator expressions should expand correctly here. 
function(setup_fortran_target target)
    # Set module directory
    set(mod_base "${CMAKE_BINARY_DIR}/mod/${target}")

    # Construct the folder - to prevent warnings of nonexistent directories
    file(MAKE_DIRECTORY ${mod_base})

    # Set module directory as property
    set_property(TARGET ${target} PROPERTY
        Fortran_MODULE_DIRECTORY "${mod_base}"
    ) 

    # Determine include directories of this target
    target_include_directories(${target} PUBLIC
        ${mod_base}
    )
endfunction()
