cmake_minimum_required(VERSION 3.11)

project(automatic_polycube)

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# autofill VORPALINE_PLATFORM
include(${CMAKE_SOURCE_DIR}/ext/geogram/cmake/geo_detect_platform.cmake)

set (GEOGRAM_LIB_ONLY ON)
set (GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ext/geogram/" CACHE PATH "full path to the Geogram (or Vorpaline) installation")

# based on https://github.com/BrunoLevy/GraphiteThree/blob/main/CMakeLists.txt
include(${CMAKE_SOURCE_DIR}/cmake/automatic_polycube.cmake)

set(
  SRCFILES
  src/containers_Geogram.cpp
  src/geometry.cpp
  src/geometry_hexahedra.cpp
  src/geometry_distortion.cpp
  src/labeling.cpp
  src/labeling_io.cpp
  src/labeling_graphcuts.cpp
  src/labeling_generators.cpp
  src/labeling_graph.cpp
  src/labeling_operators_on_invalidity.cpp
  src/labeling_operators_on_distortion.cpp
  src/io_glTF.cpp
  src/io_dump.cpp
  src/gui_labeling.cpp
  src/gui_base.cpp
)

#################
# DEPENDENCIES
#################

# Geogram
add_subdirectory(ext/geogram)

# DisjointSet
include_directories(ext/DisjointSet)

# {fmt}
add_subdirectory(ext/fmt)

# GraphCutOptimization
add_subdirectory(ext/GraphCutOptimization)

# nlohmann/json
add_subdirectory(ext/json)

# dbg-macro
include_directories(ext/dbg-macro)

# TinyGLTF
set(TINYGLTF_BUILD_LOADER_EXAMPLE OFF)
set(TINYGLTF_BUILD_GL_EXAMPLES OFF)
set(TINYGLTF_BUILD_VALIDATOR_EXAMPLE OFF)
set(TINYGLTF_BUILD_BUILDER_EXAMPLE OFF)
set(TINYGLTF_HEADER_ONLY OFF)
set(TINYGLTF_INSTALL OFF)
add_subdirectory(ext/tinygltf)
include_directories(ext/tinygltf)

# Eigen
# Thanks https://github.com/cpm-cmake/CPM.cmake/wiki/More-Snippets#eigen
include(cmake/CPM.cmake)
CPMAddPackage(
  NAME Eigen
  VERSION 3.2.8
  URL https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
  # Eigen's CMakelists are not intended for library use
  DOWNLOAD_ONLY YES 
)
if(Eigen_ADDED)
  add_library(Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

set(
  EXTERNAL_LIBS
  geogram
  geogram_gfx
  ${GLFW_LIBRARIES}
  fmt::fmt
  gco
  nlohmann_json
  tinygltf
  Eigen
)

#################
# COMPILE OPTIONS
#################

add_compile_definitions(WINDOW_SIZE="3840x2160") # bigger window size for Geogram-based applications
# colored output with Ninja
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-fcolor-diagnostics)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fdiagnostics-color)
endif()


#################
# LIBRARY
#################

add_library(automatic_polycube_lib STATIC ${SRCFILES})
target_include_directories(automatic_polycube_lib PRIVATE include)
target_link_libraries(automatic_polycube_lib PRIVATE ${EXTERNAL_LIBS})

#################
# EXECUTABLES
#################

add_executable(automatic_polycube app/automatic_polycube.cpp)
target_include_directories(automatic_polycube PRIVATE include)
target_link_libraries(automatic_polycube PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(hex_mesh_viewer app/hex_mesh_viewer.cpp)
target_include_directories(hex_mesh_viewer PRIVATE include)
target_link_libraries(hex_mesh_viewer PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(labeling_viewer app/labeling_viewer.cpp)
target_include_directories(labeling_viewer PRIVATE include)
target_link_libraries(labeling_viewer PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(labeling_painter app/labeling_painter.cpp)
target_include_directories(labeling_painter PRIVATE include)
target_link_libraries(labeling_painter PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(extract_surface app/extract_surface.cpp)
target_include_directories(extract_surface PRIVATE include)
target_link_libraries(extract_surface PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(volume_labeling app/volume_labeling.cpp)
target_include_directories(volume_labeling PRIVATE include)
target_link_libraries(volume_labeling PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(naive_labeling app/naive_labeling.cpp)
target_include_directories(naive_labeling PRIVATE include)
target_link_libraries(naive_labeling PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(graphcut_labeling app/graphcut_labeling.cpp)
target_include_directories(graphcut_labeling PRIVATE include)
target_link_libraries(graphcut_labeling PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(areas_to_tilt app/areas_to_tilt.cpp)
target_include_directories(areas_to_tilt PRIVATE include)
target_link_libraries(areas_to_tilt PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(mesh_stats app/mesh_stats.cpp)
target_include_directories(mesh_stats PRIVATE include)
target_link_libraries(mesh_stats PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(labeling_stats app/labeling_stats.cpp)
target_include_directories(labeling_stats PRIVATE include)
target_link_libraries(labeling_stats PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(halfedge_movements app/halfedge_movements.cpp)
target_include_directories(halfedge_movements PRIVATE include)
target_link_libraries(halfedge_movements PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(flip_normals app/flip_normals.cpp)
target_include_directories(flip_normals PRIVATE include)
target_link_libraries(flip_normals PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(reorient_mesh app/reorient_mesh.cpp)
target_include_directories(reorient_mesh PRIVATE include)
target_link_libraries(reorient_mesh PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(to_glTF app/to_glTF.cpp)
target_include_directories(to_glTF PRIVATE include)
target_link_libraries(to_glTF PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(polycube_distortion app/polycube_distortion.cpp)
target_include_directories(polycube_distortion PRIVATE include)
target_link_libraries(polycube_distortion PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

add_executable(convert_labeled_obj app/convert_labeled_obj.cpp)
target_include_directories(convert_labeled_obj PRIVATE include)
target_link_libraries(convert_labeled_obj PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS})

# optional test subdirectory (using googletest)
option(AUTOMATIC_POLYCUBE_BUILD_TESTS "Build automatic_polycube tests" OFF)
if(AUTOMATIC_POLYCUBE_BUILD_TESTS)

  # http://google.github.io/googletest/quickstart-cmake.html
  # https://github.com/matheusgomes28/cmake-google-tests

  enable_testing()

  option(BUILD_GMOCK OFF)
  option(INSTALL_GTEST OFF)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/a6d7fa8c0cebaa8b493a1cecdfa05018f6adc8df.zip
  )
  FetchContent_MakeAvailable(googletest)
  include_directories("${googletest_SOURCE_DIR}/googletest/include")
  add_library(GTest::GTest INTERFACE IMPORTED)
  target_link_libraries(GTest::GTest INTERFACE gtest_main)
  set_target_properties(GTest::GTest PROPERTIES COMPILE_OPTIONS "-w") # disable warnings. in particular -Wdouble-promotion

  add_executable(tests tests/main.cpp)
  target_include_directories(tests PRIVATE include)
  target_link_libraries(tests PRIVATE automatic_polycube_lib ${EXTERNAL_LIBS} GTest::GTest)

  include(GoogleTest)
  gtest_discover_tests(tests)

endif()
