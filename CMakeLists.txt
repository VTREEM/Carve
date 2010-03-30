cmake_minimum_required(VERSION 2.6)

project(carve)

set(CARVE_VERSION_MAJOR 1)
set(CARVE_VERSION_MINOR 1)
set(CARVE_VERSION_PATCH 0)

set(CARVE_VERSION ${CARVE_VERSION_MAJOR}.${CARVE_VERSION_MINOR}.${CARVE_VERSION_PATCH})



if(WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-DNOMINMAX)
endif(WIN32)

option(CARVE_WITH_GUI          "Compile gui code"                    ON)
option(CARVE_SYSTEM_BOOST      "Compile with system installed boost" ON)
option(CARVE_BOOST_COLLECTIONS "Compile with boost collections"      ON)
option(CARVE_DEBUG             "Compile in debug code"               OFF)



set(HAVE_TR1_UNORDERED_COLLECTIONS FALSE)
set(HAVE_STD_UNORDERED_COLLECTIONS FALSE)
set(HAVE_LIBSTDCPP_UNORDERED_COLLECTIONS FALSE)

set(HAVE_BOOST_UNORDERED_COLLECTIONS FALSE)


if(CARVE_SYSTEM_BOOST)
  find_package(BOOST 1.40)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
  else(Boost_FOUND)
    set(CARVE_SYSTEM_BOOST OFF)
  endif(Boost_FOUND)
endif(CARVE_SYSTEM_BOOST)
  
if(CARVE_BOOST_COLLECTIONS)
  set(HAVE_BOOST_UNORDERED_COLLECTIONS TRUE)

else(CARVE_BOOST_COLLECTIONS)
  # attempt to work out which unordered collection implementation we can use
  try_compile(_HAVE_STD_UNORDERED_COLLECTIONS
              ${CMAKE_BINARY_DIR} 
              "${carve_SOURCE_DIR}/cmake/test_std_unordered.cpp"
              OUTPUT_VARIABLE OUTPUT)
  try_compile(_HAVE_TR1_UNORDERED_COLLECTIONS
              ${CMAKE_BINARY_DIR} 
              "${carve_SOURCE_DIR}/cmake/test_tr1_unordered.cpp"
              OUTPUT_VARIABLE OUTPUT)
  try_compile(_HAVE_LIBSTDCPP_UNORDERED_COLLECTIONS
              ${CMAKE_BINARY_DIR} 
              "${carve_SOURCE_DIR}/cmake/test_libstdcpp_unordered.cpp"
              OUTPUT_VARIABLE OUTPUT)

  if(_HAVE_STD_UNORDERED_COLLECTIONS)
    set(HAVE_STD_UNORDERED_COLLECTIONS TRUE)
    message(STATUS "Using std::unordered_map")
  elseif(_HAVE_TR1_UNORDERED_COLLECTIONS)
    set(HAVE_TR1_UNORDERED_COLLECTIONS TRUE)
    message(STATUS "Using tr1::unordered_map")
  elseif(_HAVE_LIBSTDCPP_UNORDERED_COLLECTIONS)
    set(HAVE_LIBSTDCPP_UNORDERED_COLLECTIONS TRUE)
    message(STATUS "Using __gnu_cxx::unordered_map ")
  endif(_HAVE_STD_UNORDERED_COLLECTIONS)

endif(CARVE_BOOST_COLLECTIONS)

configure_file (
  "${carve_SOURCE_DIR}/include/carve/cmake-config.h.in"
  "${carve_BINARY_DIR}/include/carve/config.h"
)
include_directories(${carve_BINARY_DIR}/include)

if(CARVE_WITH_GUI)
  find_package(OpenGL)
  find_package(GLUT)

  if(NOT OPENGL_FOUND)
    message(WARNING "Unable to locate OpenGL")
    set(CARVE_WITH_GUI OFF)

  elseif(NOT GLUT_FOUND)
    message(WARNING "Unable to locate GLUT")
    set(CARVE_WITH_GUI OFF)

  else(OENGL_FOUND AND GLUT_FOUND)
    message(STATUS "Found OpenGL and GLUT")
    include_directories(${OPENGL_INCLUDE_DIR})
    include_directories(${GLUT_INCLUDE_DIR})
    add_subdirectory(external/GLEW) 
    add_subdirectory(external/GLUI) 

  endif(NOT OPENGL_FOUND)

endif(CARVE_WITH_GUI)

add_subdirectory(external/GLOOP) 

add_subdirectory(lib) 
add_subdirectory(include) 
add_subdirectory(common) 
add_subdirectory(src) 
add_subdirectory(examples) 
add_subdirectory(tests) 
