list(APPEND src_paraview_headers ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp)
list(APPEND src_paraview_sources ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.cpp)

list(APPEND src_paraview_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_paraview_sources ${src_paraview_sources} PARENT_SCOPE)
set(src_paraview_headers ${src_paraview_headers} PARENT_SCOPE)
set(src_paraview_includes ${src_paraview_includes} PARENT_SCOPE)
