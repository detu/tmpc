cmake_minimum_required (VERSION 2.8)
project (tmpc)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

unset(PROJECT_DEPENDENCY_LIBS)

find_package( Eigen3 3.2.92 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

include(ExternalProject)
option(BUILD_QPOASES "Build qpOASES along with TMPC." ON)
option(BUILD_HPMPC "Build HPMPC along with TMPC." ON)
option(BUILD_QPDUNES "Build qpDUNES along with TMPC." ON)

if(BUILD_QPDUNES)
    ExternalProject_Add(qpDUNES
                        GIT_REPOSITORY https://github.com/jfrasch/qpDUNES.git
                        INSTALL_DIR "${CMAKE_BINARY_DIR}/qpdunes-build"
                      )
    set( qpDUNES_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/qpdunes-build/include/qpdunes )
    set( qpDUNES_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/qpdunes-build/lib/libqpDUNES.a )
    set( qpDUNES_FOUND TRUE )
elseif()
    find_package(qpDUNES)
endif()
if (qpDUNES_FOUND)
    include_directories(${qpDUNES_INCLUDE_DIRS})
    list(APPEND PROJECT_DEPENDENCY_LIBS ${qpDUNES_STATIC_LIBRARIES})
endif()

if(BUILD_QPOASES)
    ExternalProject_Add(qpOASES
                        SVN_REPOSITORY https://projects.coin-or.org/svn/qpOASES/stable/3.2
                        SVN_TRUST_CERT 1
                        INSTALL_DIR "${CMAKE_BINARY_DIR}/qpoases-build"
                      )
    set( qpOASES_STATIC_LIBRARIES ${CMAKE_BINARY_DIR}/qpoases-build/lib/libqpOASES.a )
    set( qpOASES_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/qpoases-build/include )
    set( qpOASES_FOUND TRUE )
elseif()
    find_package(qpOASES)
endif()
if (qpOASES_FOUND)
    include_directories(${qpOASES_INCLUDE_DIRS})
    list(APPEND PROJECT_DEPENDENCY_LIBS ${qpOASES_STATIC_LIBRARIES})
endif()

if(BUILD_HPMPC)
    ExternalProject_Add(hpmpc
                        GIT_REPOSITORY https://github.com/mkotlyar/hpmpc.git
                        BUILD_IN_SOURCE 1
                        CONFIGURE_COMMAND ""
                        BUILD_COMMAND make
                        INSTALL_COMMAND make PREFIX="${CMAKE_BINARY_DIR}/hpmpc-build" install_static
                      )
    set( HPMPC_STATIC_LIB ${CMAKE_BINARY_DIR}/hpmpc-build/hpmpc/lib/libhpmpc.a )
    set( HPMPC_INCLUDE_DIR ${CMAKE_BINARY_DIR}/hpmpc-build/hpmpc/include )
    set( HPMPC_FOUND TRUE )
elseif()
    find_package(HPMPC)
endif()
if (HPMPC_FOUND)
    include_directories(${HPMPC_INCLUDE_DIR})
    set(PROJECT_DEPENDENCY_LIBS ${PROJECT_DEPENDENCY_LIBS} ${HPMPC_STATIC_LIB})
endif()

message(STATUS "PROJECT_DEPENDENCY_LIBS=" ${PROJECT_DEPENDENCY_LIBS})

include_directories(${PROJECT_INCLUDE_DIR})

include(CompilerOptions.cmake)

add_subdirectory(src)
add_subdirectory(test)
