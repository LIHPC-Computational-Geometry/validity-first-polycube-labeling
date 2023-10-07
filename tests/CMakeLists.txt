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

add_executable(tests main.cpp)
target_link_libraries(tests PRIVATE geogram geogram_gfx ${GLFW_LIBRARIES} fmt::fmt gco nlohmann_json GTest::GTest)

include(GoogleTest)
gtest_discover_tests(tests)