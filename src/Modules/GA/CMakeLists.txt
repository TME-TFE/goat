# gamod_driver
add_library(gamod_driver STATIC
    gamod_driver.F90)
setup_fortran_target(gamod_driver)
target_link_libraries(gamod_driver PUBLIC 
    mod_triangulation
    UnstructuredInterpolation2D
    goatmod_types
    goatmod_userinput
    gamod_utils)

# gamod_utils - here we lump all utilities of grid adaptation modules
add_library(gamod_utils STATIC
    gamod_math.F90 gamod_types.F90 gamod_utility.F90)
setup_fortran_target(gamod_utils)
target_link_libraries(gamod_utils PUBLIC 
    mod_precision
    mod_sort   
    mod_readwrite
    mod_polygon
    mod_dynamicarrays
    mod_gradient 
    mod_structured2Dgridding 
    mod_triangulation  
    goatmod_types 
    DistributionFunction
    UnstructuredInterpolation2D
    goatmod_userinput)