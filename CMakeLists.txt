cmake_minimum_required (VERSION 2.8)
project (tmpc)

# Enable modern C++
set(CMAKE_CXX_STANDARD 14)

# Don't use compiler's C++ extensions
set(CMAKE_CXX_EXTENSIONS False)

# Include CMakeToolsHelpers for easier debugging from VS Code
include(CMakeToolsHelpers OPTIONAL)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH})
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

# List of QP solvers
set(qp_solvers "qpOASES" "HPMPC" "HPIPM")

# List of selectable extern libs
set(extern_libs ${qp_solvers} "BLASFEO")

# Create WITH_${lib} cached variables and find selected packages
foreach (lib ${extern_libs})
    set(WITH_${lib} OFF CACHE BOOL "Link to ${lib}")

    if (WITH_${lib})
        find_package(${lib} REQUIRED)
    endif ()

    message(STATUS "WITH_${lib}=${WITH_${lib}}")
endforeach ()

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

add_subdirectory(src)
add_subdirectory(test)
