
# mod_diagnostics
add_library(mod_diagnostics STATIC mod_diagnostics.F90)
setup_fortran_target(mod_diagnostics)
target_link_libraries(mod_diagnostics PUBLIC 
    mod_precision
    mod_sparseinterface
    mod_plotter
    mod_specialchars)

# mod_dynamicarrays
add_library(mod_dynamicarrays STATIC mod_dynamicarrays.F90)
setup_fortran_target(mod_dynamicarrays)
target_link_libraries(mod_dynamicarrays PUBLIC 
    mod_precision
    mod_errorhandler
    mod_sparseinterface)

# mod_errorhandler
add_library(mod_errorhandler STATIC mod_errorhandler.F90)
setup_fortran_target(mod_errorhandler)
target_link_libraries(mod_errorhandler PUBLIC 
    mod_precision)

# mod_inputfileparser
add_library(mod_inputfileparser STATIC mod_inputfileparser.F90)
setup_fortran_target(mod_inputfileparser)
target_link_libraries(mod_inputfileparser PUBLIC 
    mod_precision
    mod_specialchars 
    mod_readwrite 
    mod_errorhandler)

# mod_linearsolverinterface
add_library(mod_linearsolverinterface STATIC mod_linearsolverinterface.F90)
setup_fortran_target(mod_linearsolverinterface)
target_link_libraries(mod_linearsolverinterface PUBLIC 
    mod_precision
    mod_sparseinterface
    mod_errorhandler
    mod_inputfileparser)
target_link_libraries(mod_linearsolverinterface PRIVATE # Need to link the linear solvers too here
    ${UMFPACK_LIB}
    ${CHOLMOD_LIB}
    ${AMD_LIB}
    ${COLAMD_LIB}
    ${SUITESPARSECONFIG_LIB}

    BLAS::BLAS
    LAPACK::LAPACK

    OpenMP::OpenMP_Fortran) 

# mod_plotter
add_library(mod_plotter STATIC mod_plotter.F90)
setup_fortran_target(mod_plotter)
target_link_libraries(mod_plotter PUBLIC 
    mod_precision
    mod_specialchars)

# mod_polygon
add_library(mod_polygon STATIC mod_polygon.F90)
setup_fortran_target(mod_polygon)
target_link_libraries(mod_polygon PUBLIC 
    mod_precision
    mod_errorhandler
    mod_plotter
    mod_sort
    mod_constants)

# mod_readwrite
add_library(mod_readwrite STATIC mod_readwrite.F90)
setup_fortran_target(mod_readwrite)
target_link_libraries(mod_readwrite PUBLIC 
    mod_precision
    mod_errorhandler
    mod_specialchars )

# mod_search
add_library(mod_search STATIC mod_search.F90)
setup_fortran_target(mod_search)
target_link_libraries(mod_search PUBLIC 
    mod_precision
    mod_errorhandler)

# mod_sort
add_library(mod_sort STATIC mod_sort.F90)
setup_fortran_target(mod_sort)
target_link_libraries(mod_sort PUBLIC 
    mod_precision
    mod_errorhandler)

# mod_sparseinterface
add_library(mod_sparseinterface STATIC mod_sparseinterface.F90)
setup_fortran_target(mod_sparseinterface)
target_link_libraries(mod_sparseinterface PUBLIC 
    mod_precision
    mod_errorhandler
    clayer)

# mod_utility
add_library(mod_utility STATIC mod_utility.F90)
setup_fortran_target(mod_utility)
target_link_libraries(mod_utility PUBLIC 
    mod_precision)

# Library with all libraries in General
add_library(General INTERFACE)
target_link_libraries(General INTERFACE 
    mod_diagnostics
    mod_dynamicarrays
    mod_errorhandler
    mod_inputfileparser
    mod_linearsolverinterface
    mod_plotter
    mod_polygon
    mod_readwrite
    mod_search
    mod_sort
    mod_sparseinterface
    mod_utility)

