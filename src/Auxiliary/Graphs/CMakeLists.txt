add_library(mod_graph STATIC
    mod_graph.F90)

setup_fortran_target(mod_graph)

target_link_libraries(mod_graph PUBLIC 
    mod_precision
    mod_errorhandler
    mod_sort
    mod_dynamicarrays)

# Library of all routines in Graphs
add_library(Graphs INTERFACE)
target_link_libraries(Graphs INTERFACE 
    mod_graph)