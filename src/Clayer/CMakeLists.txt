# C routines
add_library(clayer_c STATIC
    ConstructMyCSparse.c ConvertMyCSparseToCS.c SpMM.c SpMMBasic.c SpMMV.c)

# Add the headers
target_include_directories(clayer_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Include)

# Link to SuiteSparse
target_include_directories(clayer_c PRIVATE ${SUITESPARSE_INCLUDE_DIR}) # to make sure it finds cs.h
target_link_libraries(clayer_c PRIVATE ${CXSPARSE_LIB}) #to make sure the suitesparse routines are found

# Fortran interface
add_library(clayer STATIC Clayer.F90 CSparseF.F90 CUtilities.F90)

setup_fortran_target(clayer)

target_link_libraries(clayer PUBLIC clayer_c)