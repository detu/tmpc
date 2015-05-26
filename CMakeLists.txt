cmake_minimum_required (VERSION 2.8)
project (RTMC)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
#SET(WABS_ROOT ${PROJECT_SOURCE_DIR}/.. CACHE PATH "WABS root directory")

set(RTMC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
#set(EXTERN_DIR ${WABS_ROOT}/extern CACHE STRING "3rd-party software directory")
#set(EXTERN_INCLUDE_DIR ${EXTERN_DIR}/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#SET(MATLAB_ROOT ${EXTERN_DIR}/MATLAB CACHE PATH "Path to MATLAB")
#FIND_PACKAGE(MATLAB REQUIRED)

find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

find_package( qpDUNES REQUIRED )
include_directories(${qpDUNES_INCLUDE_DIRS})

include_directories(${RTMC_INCLUDE_DIR})

add_subdirectory(src)
