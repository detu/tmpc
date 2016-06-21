cmake_minimum_required (VERSION 2.8)
project (MPMC)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
#SET(WABS_ROOT ${PROJECT_SOURCE_DIR}/.. CACHE PATH "WABS root directory")

set(RTMC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
#set(EXTERN_DIR ${WABS_ROOT}/extern CACHE STRING "3rd-party software directory")
#set(EXTERN_INCLUDE_DIR ${EXTERN_DIR}/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#SET(MATLAB_ROOT ${EXTERN_DIR}/MATLAB CACHE PATH "Path to MATLAB")
#FIND_PACKAGE(MATLAB REQUIRED)

set(TMPC_DEPENDENCY_LIBS)

find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

find_package( qpDUNES )
if (qpDUNES_FOUND)
    include_directories(${qpDUNES_INCLUDE_DIRS})
    list(APPEND TMPC_DEPENDENCY_LIBS ${qpDUNES_STATIC_LIBRARIES})
endif()

find_package( qpOASES )
if (qpOASES_FOUND)
    include_directories(${qpOASES_INCLUDE_DIRS})
    list(APPEND TMPC_DEPENDENCY_LIBS ${qpOASES_STATIC_LIBRARIES})
endif()

find_package(HPMPC)
if (HPMPC_FOUND)
    include_directories(${HPMPC_INCLUDE_DIR})
    set(TMPC_DEPENDENCY_LIBS ${TMPC_DEPENDENCY_LIBS} ${HPMPC_STATIC_LIB})
endif()

message(STATUS "TMPC_DEPENDENCY_LIBS=" ${TMPC_DEPENDENCY_LIBS})

include_directories(${RTMC_INCLUDE_DIR})

include(CompilerOptions.cmake)

add_subdirectory(src)
add_subdirectory(test)