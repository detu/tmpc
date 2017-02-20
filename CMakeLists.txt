cmake_minimum_required (VERSION 2.8)
project (tmpc)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include(ExternalProject)
option(BUILD_QPDUNES "Build qpDUNES along with TMPC." ON)

set(EXTERN_DIR "${CMAKE_CURRENT_BINARY_DIR}/extern")
#set(EP_PREFIX ${CMAKE_BINARY_DIR}/extern)

include_directories("${EXTERN_DIR}/include")
link_directories("${EXTERN_DIR}/lib")

if(BUILD_QPDUNES)
    ExternalProject_Add(qpDUNES
                        PREFIX ${EXTERN_DIR}
                        GIT_REPOSITORY https://github.com/mkotlyar/qpDUNES.git
                        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${EXTERN_DIR}" "-DQPDUNES_MAKE_EXAMPLES=OFF"
                      )
    set( qpDUNES_INCLUDE_DIRS "${EXTERN_DIR}/include/qpdunes" )
    set( qpDUNES_LIBRARY_DIRS "${EXTERN_DIR}/lib/libqpdunes.a" )
    set( qpDUNES_FOUND TRUE )
elseif()
    find_package(qpDUNES)
endif()

find_package(qpOASES)
find_package(HPMPC)

add_subdirectory(src)
add_subdirectory(test)
