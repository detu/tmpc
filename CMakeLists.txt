cmake_minimum_required (VERSION 2.8)
project (tmpc)

# Enable modern C++
set(CMAKE_CXX_STANDARD 14)

# Don't use compiler's C++ extensions
set(CMAKE_CXX_EXTENSIONS False)

# Include CMakeToolsHelpers for easier debugging from VS Code
include(CMakeToolsHelpers OPTIONAL)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# Cached CMake variables.
set(WITH_BLASFEO OFF CACHE BOOL "Link to BLASFEO")
set(WITH_HPMPC OFF CACHE BOOL "Link to HPMPC")
set(WITH_HPIPM OFF CACHE BOOL "Link to HPIPM")

message(STATUS "WITH_HPMPC=${WITH_HPMPC}")
message(STATUS "WITH_HPIPM=${WITH_HPIPM}")
message(STATUS "WITH_BLASFEO=${WITH_BLASFEO}")

if (WITH_HPIPM AND NOT WITH_BLASFEO)
    message(SEND_ERROR "If you choose WITH_HPIPM, you have to choose WITH_BLASFEO as well.")
endif ()

find_package(Eigen3 3.2.10)

if (EIGEN3_FOUND)
    message(STATUS "Eigen3 Found in ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen3 NOT FOUND")
endif()

if (BLAZE_DIR)
    include_directories(${BLAZE_DIR})
    message(STATUS "BLAZE_DIR=${BLAZE_DIR}")
endif()

find_package(qpOASES)

if (WITH_HPMPC)
    find_package(HPMPC REQUIRED)
    add_definitions("-DWITH_HPMPC")
endif ()

if (WITH_HPIPM)
    find_package(HPIPM REQUIRED)
    add_definitions("-DWITH_HPIPM")
endif ()

if (WITH_BLASFEO)
    find_package(BLASFEO REQUIRED)
    add_definitions("-DWITH_BLASFEO")
endif ()

add_subdirectory(src)
add_subdirectory(test)
