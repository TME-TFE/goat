# Goat program
#-------------
# Add main program as target
add_executable(goat
    Goat.F90)

setup_fortran_target(goat)

# Link
target_link_libraries(goat PRIVATE 
    goatmod_userinput 
    mod_global_environment 
    mod_plotter 
    goat_drivers)

# Link to additional external libraries
target_link_libraries(goat PRIVATE 
    ${CXSPARSE_LIB}
    ${UMFPACK_LIB}
    ${CHOLMOD_LIB}
    ${AMD_LIB}
    ${COLAMD_LIB}
    ${SUITESPARSECONFIG_LIB}

    BLAS::BLAS
    LAPACK::LAPACK

    goat_openmp)

# Installation
install(TARGETS goat
    RUNTIME DESTINATION .)

# Shapeopt program
#-----------------
# Add shape optimization as target
add_executable(shapeopt
    ShapeOptimization.F90)

setup_fortran_target(shapeopt)

# Link
target_link_libraries(shapeopt PRIVATE 
    somod_userinput
    mod_global_environment
    mod_plotter
    ShapeOptDriver)

# Link to additional external libraries
target_link_libraries(shapeopt PRIVATE 
    ${CXSPARSE_LIB}
    ${UMFPACK_LIB}
    ${CHOLMOD_LIB}
    ${AMD_LIB}
    ${COLAMD_LIB}
    ${SUITESPARSECONFIG_LIB}

    BLAS::BLAS
    LAPACK::LAPACK

    goat_openmp)

# Installation
install(TARGETS shapeopt
    RUNTIME DESTINATION .)