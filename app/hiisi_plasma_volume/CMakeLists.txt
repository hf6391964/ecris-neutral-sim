# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( main )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package
include_directories( BEFORE "${PROJECT_SOURCE_DIR}/../../include"
                     "${PROJECT_SOURCE_DIR}/../../../gsl-2.2.1")


# Creating entries for target: main
# ############################
add_executable( main main.cpp )
target_link_libraries(main CGAL ${CGAL_LIBRARIES} pthread cblas "${PROJECT_SOURCE_DIR}/../../libneutrals.a" "${PROJECT_SOURCE_DIR}/../../../gsl-2.2.1/.libs/libgsl.a")

