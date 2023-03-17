  # based on https://github.com/BrunoLevy/GraphiteThree/blob/main/cmake/graphite.cmake
  
  ##############################################################################
  
  set(AUTOMATIC_POLYCUBE_SOURCE_DIR ${CMAKE_SOURCE_DIR})
  include(${AUTOMATIC_POLYCUBE_SOURCE_DIR}/cmake/automatic_polycube_config.cmake)
  
  # CMakeOptions is included after (so that it is 
  # possible to override user-editable variables in
  # it instead of using CMakeGUI)
  
  include(${AUTOMATIC_POLYCUBE_SOURCE_DIR}/ext/geogram/cmake/geogram.cmake)
  
  ##############################################################################
  
  link_directories(${AUTOMATIC_POLYCUBE_SOURCE_DIR}/${RELATIVE_LIB_DIR})
  
  ##############################################################################
  