
# Add libraries
#--------------
# mod_gradient
add_library(mod_gradient STATIC
    mod_gradient.F90)

setup_fortran_target(mod_gradient)

target_link_libraries(mod_gradient PUBLIC 
    mod_precision mod_errorhandler mod_linearsolverinterface)

# mod_streamlinetracer
add_library(mod_streamlinetracing2D STATIC
    mod_streamlinetracing2D.F90)

setup_fortran_target(mod_streamlinetracing2D)

target_link_libraries(mod_streamlinetracing2D PUBLIC 
    mod_precision mod_errorhandler mod_linearsolverinterface mod_dynamicarrays mod_structured2Dgridding)
target_link_libraries(mod_streamlinetracing2D PRIVATE 
    goat_openmp)

# mod_structured2Dgridding
add_library(mod_structured2Dgridding STATIC
    mod_structured2Dgridding.F90)

setup_fortran_target(mod_structured2Dgridding)

target_link_libraries(mod_structured2Dgridding PUBLIC 
    mod_precision mod_errorhandler mod_sort)

# mod_triangulation
add_library(mod_triangulation STATIC
    mod_triangulation.F90)

setup_fortran_target(mod_triangulation)

target_link_libraries(mod_triangulation PUBLIC 
    mod_precision
    mod_errorhandler
    mod_readwrite
    mod_gradient
    mod_sort
    mod_utility
    mod_plotter 
    mod_specialchars)

# PolygonShapeFunction
add_library(PolygonShapeFunction STATIC
    PolygonShapeFunction.F90)

setup_fortran_target(PolygonShapeFunction)

target_link_libraries(PolygonShapeFunction PUBLIC 
    mod_precision
    mod_plotter
    mod_errorhandler
    mod_structured2Dgridding)

# Library with all libraries of the Auxiliary folder (subfolders included)
add_library(Auxiliary INTERFACE)
target_link_libraries(Auxiliary INTERFACE 
    mod_gradient
    mod_streamlinetracing2D
    mod_structured2Dgridding
    mod_triangulation
    Contour
    Graphs
    Interpolation)

# Add subdirectories
#-------------------
add_subdirectory(Contour)
add_subdirectory(Graphs)
add_subdirectory(Interpolation)
