# based on https://github.com/BrunoLevy/GraphiteThree/blob/main/CMakeLists.txt

cmake_minimum_required(VERSION 3.11)

project(automatic_polycube)

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

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

# nlohmann/json
add_subdirectory(ext/json)

# optional test subdirectory (using googletest)
option(AUTOMATIC_POLYCUBE_BUILD_TESTS "Build automatic_polycube tests" OFF)
if(AUTOMATIC_POLYCUBE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

#################
# VARIABLES
#################

add_compile_definitions(WINDOW_SIZE="3840x2160")

#################
# EXECUTABLES
#################

add_executable(automatic_polycube ${SRCFILES} app/automatic_polycube.cpp)
target_link_libraries(automatic_polycube geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(hex_mesh_viewer ${SRCFILES} app/hex_mesh_viewer.cpp)
target_link_libraries(hex_mesh_viewer geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(labeling_viewer ${SRCFILES} app/labeling_viewer.cpp)
target_link_libraries(labeling_viewer geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(labeling_painter ${SRCFILES} app/labeling_painter.cpp)
target_link_libraries(labeling_painter geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(extract_surface ${SRCFILES} app/extract_surface.cpp)
target_link_libraries(extract_surface geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(volume_labeling ${SRCFILES} app/volume_labeling.cpp)
target_link_libraries(volume_labeling geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(naive_labeling ${SRCFILES} app/naive_labeling.cpp)
target_link_libraries(naive_labeling geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(graphcut_labeling ${SRCFILES} app/graphcut_labeling.cpp)
target_link_libraries(graphcut_labeling geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(ovm_io ${SRCFILES} app/ovm_io.cpp)
target_link_libraries(ovm_io geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(areas_to_tilt ${SRCFILES} app/areas_to_tilt.cpp)
target_link_libraries(areas_to_tilt geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(mesh_stats ${SRCFILES} app/mesh_stats.cpp)
target_link_libraries(mesh_stats geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(labeling_stats ${SRCFILES} app/labeling_stats.cpp)
target_link_libraries(labeling_stats geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(halfedge_movements ${SRCFILES} app/halfedge_movements.cpp)
target_link_libraries(halfedge_movements geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)

add_executable(flip_normals ${SRCFILES} app/flip_normals.cpp)
target_link_libraries(flip_normals geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json)