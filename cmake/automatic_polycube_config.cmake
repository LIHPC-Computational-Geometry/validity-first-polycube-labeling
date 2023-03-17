# based on https://github.com/BrunoLevy/GraphiteThree/blob/main/cmake/graphite_config.cmake

# User-configurable variables:
# - Path to Geogram

# Detecting Python
macro(automatic_polycube_find_Python)
   # Use version of Python compiled with:
   # ./configure --with-pydebug --without-pymalloc --prefix /opt/debugpython/
   #  --enable-shared
   #set(WHERE_ARE_PYTHON_INCLUDES /opt/debugpython/include/python3.7d/)
   #set(WHERE_IS_PYTHON_LIB /opt/debugpython/lib/libpython3.7d.so) 
   find_package(PythonLibs 3 QUIET)
   if(NOT PYTHONLIBS_FOUND)
      message(
         STATUS
	"CMake did not find Python library, 
         using default fallbacks (edit WHERE_IS... in CMakeGUI if need be)."	
      )
      set(PYTHON_INCLUDE_DIRS ${WHERE_ARE_PYTHON_INCLUDES})
      set(PYTHON_LIBRARIES ${WHERE_IS_PYTHON_LIB})
   endif()
   if(
      NOT "${PYTHON_INCLUDE_DIRS}" STREQUAL "" AND
      NOT "${PYTHON_LIBRARIES}" STREQUAL ""
   )
      set(AUTOMATIC_POLYCUBE_FOUND_PYTHON TRUE)
   endif()
endmacro()

# Path to geogram
#################

if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/ext/geogram/)
   set(
      GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ext/geogram/"
      CACHE PATH "full path to the Geogram (or Vorpaline) installation"
   )
   set(USE_BUILTIN_GEOGRAM TRUE)
else()
   message(
      SEND_ERROR
      "CMake did not find Geogram in ${CMAKE_SOURCE_DIR}/ext/geogram/"	
      )
endif()