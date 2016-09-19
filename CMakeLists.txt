cmake_minimum_required (VERSION 2.8)
project (tmpc)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package( Eigen3 3.2.92 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

include(ExternalProject)
option(BUILD_QPOASES "Build qpOASES along with TMPC." ON)
option(BUILD_HPMPC "Build HPMPC along with TMPC." ON)
option(BUILD_QPDUNES "Build qpDUNES along with TMPC." ON)

set(EXTERN_DIR ${CMAKE_BINARY_DIR}/extern)
#set(EP_PREFIX ${CMAKE_BINARY_DIR}/extern)

if(BUILD_QPDUNES)
    ExternalProject_Add(qpDUNES
                        PREFIX ${EXTERN_DIR}
                        GIT_REPOSITORY https://github.com/jfrasch/qpDUNES.git
                        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${EXTERN_DIR}"
                      )
    set( qpDUNES_INCLUDE_DIRS ${EXTERN_DIR}/include/qpdunes )
    set( qpDUNES_LIBRARY_DIRS ${EXTERN_DIR}/lib/libqpdunes.a )
    set( qpDUNES_FOUND TRUE )
elseif()
    find_package(qpDUNES)
endif()

if(BUILD_QPOASES)
    ExternalProject_Add(qpOASES
                        PREFIX ${EXTERN_DIR}
                        SVN_REPOSITORY https://projects.coin-or.org/svn/qpOASES/stable/3.2
                        SVN_TRUST_CERT 1
                        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${EXTERN_DIR}"
                      )
    set( qpOASES_STATIC_LIBRARIES ${EXTERN_DIR}/lib/libqpOASES.a )
    set( qpOASES_INCLUDE_DIRS ${EXTERN_DIR}/include )
    set( qpOASES_FOUND TRUE )
elseif()
    find_package(qpOASES)
endif()

if(BUILD_HPMPC)
    ExternalProject_Add(hpmpc
                        PREFIX ${EXTERN_DIR}
                        GIT_REPOSITORY https://github.com/mkotlyar/hpmpc.git
                        BUILD_IN_SOURCE 1
                        CONFIGURE_COMMAND ""
                        BUILD_COMMAND make static_library
                        BUILD_BYPRODUCTS "${EXTERN_DIR}/src/hpmpc/libhpmpc.a"
                        INSTALL_COMMAND make PREFIX="${EXTERN_DIR}" install_static
                      )
    set( HPMPC_STATIC_LIB ${EXTERN_DIR}/hpmpc/lib/libhpmpc.a )
    set( HPMPC_INCLUDE_DIR ${EXTERN_DIR}/hpmpc/include )
    set( HPMPC_FOUND TRUE )
elseif()
    find_package(HPMPC)
endif()

include_directories(${PROJECT_INCLUDE_DIR})

include(CompilerOptions.cmake)

add_subdirectory(src)
add_subdirectory(test)
