# based on https://github.com/BrunoLevy/GraphiteThree/blob/main/CMakeLists.txt

cmake_minimum_required(VERSION 3.11)

project(automatic_polycube)

set (CMAKE_CXX_STANDARD 20)

set(
  GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ext/geogram/"
  CACHE PATH "full path to the Geogram (or Vorpaline) installation"
)

include(${CMAKE_SOURCE_DIR}/cmake/automatic_polycube.cmake)

file(GLOB SRCFILES src/*.cpp )
include_directories(include)

set_source_files_properties(app/graphcut_labeling.cpp PROPERTIES COMPILE_OPTIONS "-Wdouble-promotion")

#################
# DEPENDENCIES
#################

# Geogram
add_subdirectory(ext/geogram)

# DisjointSet
include_directories(ext/DisjointSet)

# {fmt}
add_subdirectory(ext/fmt)
include_directories(ext/fmt/include) # why needed ?

# GraphCutOptimization
add_subdirectory(ext/GraphCutOptimization)
include_directories(ext/GraphCutOptimization/src) # why needed ?
set_source_files_properties(ext/GraphCutOptimization/src/GCoptimization.h PROPERTIES COMPILE_OPTIONS "-Wnon-virtual-dtor;-Wunused-parameter")

# nlohmann/json
add_subdirectory(ext/json)

#################
# VARIABLES
#################

add_compile_definitions(WINDOW_SIZE="3840x2160")

#################
# EXECUTABLES
#################

add_executable(automatic_polycube ${SRCFILES} app/automatic_polycube)
target_link_libraries(automatic_polycube geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(hex_mesh_viewer ${SRCFILES} app/hex_mesh_viewer)
target_link_libraries(hex_mesh_viewer geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(labeling_viewer ${SRCFILES} app/labeling_viewer)
target_link_libraries(labeling_viewer geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(labeling_painter ${SRCFILES} app/labeling_painter)
target_link_libraries(labeling_painter geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(extract_surface ${SRCFILES} app/extract_surface)
target_link_libraries(extract_surface geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(volume_labeling ${SRCFILES} app/volume_labeling)
target_link_libraries(volume_labeling geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(naive_labeling ${SRCFILES} app/naive_labeling)
target_link_libraries(naive_labeling geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(graphcut_labeling ${SRCFILES} app/graphcut_labeling)
target_link_libraries(graphcut_labeling geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(ovm_io ${SRCFILES} app/ovm_io)
target_link_libraries(ovm_io geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(areas_to_tilt ${SRCFILES} app/areas_to_tilt)
target_link_libraries(areas_to_tilt geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)