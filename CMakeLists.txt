cmake_minimum_required (VERSION 2.8)
project (tmpc)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(NOT TARGET gtest)
    add_subdirectory(test/googletest)
endif()

unset(PROJECT_DEPENDENCY_LIBS)

find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(qpDUNES)
if (qpDUNES_FOUND)
    include_directories(${qpDUNES_INCLUDE_DIRS})
    list(APPEND PROJECT_DEPENDENCY_LIBS ${qpDUNES_STATIC_LIBRARIES})
endif()

find_package(qpOASES)
if (qpOASES_FOUND)
    include_directories(${qpOASES_INCLUDE_DIRS})
    list(APPEND PROJECT_DEPENDENCY_LIBS ${qpOASES_STATIC_LIBRARIES})
endif()

set(HPMPC_BUILD_TESTS OFF)    # disable tests for hpmpc to avoid duplication of gtest
add_subdirectory(extern/hpmpc)
set(HPMPC_FOUND TRUE)    # notify the other CMake scripts that hpmpc is available
set(HPMPC_INCLUDE_DIR extern/hpmpc/include)
set(HPMPC_STATIC_LIB hpmpc)
include_directories(${HPMPC_INCLUDE_DIR})
set(PROJECT_DEPENDENCY_LIBS ${PROJECT_DEPENDENCY_LIBS} ${HPMPC_STATIC_LIB})

message(STATUS "PROJECT_DEPENDENCY_LIBS=" ${PROJECT_DEPENDENCY_LIBS})

include_directories(${PROJECT_INCLUDE_DIR})

include(CompilerOptions.cmake)

add_subdirectory(src)
add_subdirectory(test)