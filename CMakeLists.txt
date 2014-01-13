cmake_minimum_required(VERSION 2.6)

project(carve)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(carve_VERSION_MAJOR 2)
set(carve_VERSION_MINOR 0)
set(carve_VERSION_PATCH 0a)

set(CARVE_VERSION ${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}.${carve_VERSION_PATCH})

option(BUILD_COVERAGE                    "Compile with gcov"                                 OFF)
option(BUILD_SHARED_LIBS                 "Compile libcarve as shared"                        ON)
option(CARVE_WITH_GUI                    "Compile gui code"                                  ON)
option(CARVE_SYSTEM_BOOST                "Compile with system installed boost"               ON)
option(CARVE_BOOST_COLLECTIONS           "Compile with boost collections"                    ON)
option(CARVE_DEBUG                       "Compile in debug code"                             OFF)
option(CARVE_DEBUG_WRITE_PLY_DATA        "Write geometry output during debug"                OFF)
option(CARVE_USE_EXACT_PREDICATES        "Use Shewchuk's exact predicates, where possible"   OFF)
option(CARVE_INTERSECT_GLU_TRIANGULATOR  "Include support for GLU triangulator in intersect" OFF)
option(CARVE_GTEST_TESTS                 "Complie gtest, and dependent tests"                ON)

if (MSVC)
  # For MSVC, CMake sets certain flags to defaults we want to override.
  # This replacement code is taken from sample in the CMake Wiki at
  # http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace.
  foreach (flag_var
           CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
           CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    # In hermetic build environments, tests may not have access to MS runtime
    # DLLs, so this replaces /MD (CRT libraries in DLLs) with /MT (static CRT
    # libraries).
    string(REPLACE "/MD" "-MT" ${flag_var} "${${flag_var}}")
    # We prefer more strict warning checking for building Google Test.
    # Replaces /W3 with /W4 in defaults.
    string(REPLACE "/W3" "-W4" ${flag_var} "${${flag_var}}")
  endforeach()
endif()

SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (BUILD_COVERAGE)
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_BUILD_TYPE DEBUG)
  IF (CMAKE_COMPILER_IS_GNUCC)
    SET(CMAKE_C_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
  ENDIF()
endif(BUILD_COVERAGE)

if(WIN32)
  set(BUILD_SHARED_LIBS OFF) # until everything is exported
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-DNOMINMAX)
endif(WIN32)

find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(HAVE_TR1_UNORDERED_COLLECTIONS FALSE)
set(HAVE_STD_UNORDERED_COLLECTIONS FALSE)
set(HAVE_LIBSTDCPP_UNORDERED_COLLECTIONS FALSE)

set(HAVE_BOOST_UNORDERED_COLLECTIONS FALSE)

if(CARVE_SYSTEM_BOOST)
  find_package(Boost 1.40)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "Using system boost")
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

if(CARVE_WITH_GUI)
  find_package(OpenGL)
  find_package(GLUT)
  find_package(GLEW)

  if(NOT OPENGL_FOUND)
    message(WARNING "Unable to locate OpenGL")
    set(CARVE_WITH_GUI OFF)

  elseif(NOT GLUT_FOUND)
    message(WARNING "Unable to locate GLUT")
    set(CARVE_WITH_GUI OFF)

  else(OPENGL_FOUND AND GLUT_FOUND)
    message(STATUS "Found OpenGL and GLUT")
    include_directories(${OPENGL_INCLUDE_DIR})
    include_directories(${GLUT_INCLUDE_DIR})

    if(GLEW_FOUND)
    else(GLEW_FOUND)
      if(WIN32)
        add_definitions(-DGLEW_STATIC)
      endif(WIN32)
      add_subdirectory(external/GLEW) 
    endif(GLEW_FOUND)

    if(WIN32)
      add_definitions(-DGLUI_NO_LIB_PRAGMA)
      add_definitions(-DGLUI_USE_STATIC_LIB)
    endif(WIN32)
    add_subdirectory(external/GLUI) 

  endif(NOT OPENGL_FOUND)

endif(CARVE_WITH_GUI)

add_subdirectory(external/GLOOP) 

if (CARVE_GTEST_TESTS)
  add_subdirectory(external/gtest-1.5.0)
endif(CARVE_GTEST_TESTS)

configure_file (
  "${carve_SOURCE_DIR}/include/carve/cmake-config.h.in"
  "${carve_BINARY_DIR}/include/carve/config.h"
)
add_definitions(-DCMAKE_BUILD)

include_directories(${carve_BINARY_DIR}/include)

add_subdirectory(lib) 
add_subdirectory(include) 
add_subdirectory(common) 
add_subdirectory(src) 
add_subdirectory(examples) 
add_subdirectory(tests) 


include(CTest)

if(BUILD_TESTING)
  ENABLE_TESTING()
endif(BUILD_TESTING)
