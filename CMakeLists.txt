# based on https://github.com/BrunoLevy/GraphiteThree/blob/main/CMakeLists.txt

cmake_minimum_required(VERSION 3.0)

project(automatic_polycube)

set(
  GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ext/geogram/"
  CACHE PATH "full path to the Geogram (or Vorpaline) installation"
)

include(${CMAKE_SOURCE_DIR}/cmake/automatic_polycube.cmake)

file(GLOB SRCFILES src/*.cpp )
include_directories(include)

add_subdirectory(ext/geogram)
include_directories(ext/DisjointSet)
add_subdirectory(ext/fmt)
include_directories(ext/fmt/include) # why needed ?
link_libraries(fmt::fmt) # why needed ?

add_executable(labeling_viewer ${SRCFILES} app/labeling_viewer)
target_link_libraries(labeling_viewer geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt)

add_executable(hex_mesh_viewer ${SRCFILES} app/hex_mesh_viewer)
target_link_libraries(hex_mesh_viewer geogram geogram_gfx ${GLFW_LIBRARIES}  fmt::fmt)
