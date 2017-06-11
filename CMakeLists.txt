cmake_minimum_required (VERSION 2.8)
project (tmpc)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

find_package(Eigen3 3.2.10 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(qpOASES)
find_package(HPMPC)

add_subdirectory(src)
add_subdirectory(test)
