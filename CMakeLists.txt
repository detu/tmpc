cmake_minimum_required (VERSION 3.10)   # Need at least 3.10 for gtest_discover_tests()
project (tmpc)

# Enable modern C++
set(CMAKE_CXX_STANDARD 17)

# Don't use compiler's C++ extensions
set(CMAKE_CXX_EXTENSIONS False)

# Enable position-independent code on all targets by default.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Include CMakeToolsHelpers for easier debugging from VS Code
include(CMakeToolsHelpers OPTIONAL)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_DEBUG_POSTFIX d)

# List of QP solvers
set(qp_solvers "qpOASES" "HPMPC" "HPIPM" "treeQP")

# List of selectable extern libs
set(extern_libs ${qp_solvers})

# Create TMPC_WITH_${lib} cached variables and find selected packages
foreach (lib ${extern_libs})
    set(TMPC_WITH_${lib} OFF CACHE BOOL "Link to ${lib}")

    if (TMPC_WITH_${lib})
        find_package(${lib} REQUIRED)
    endif ()

    message(STATUS "TMPC_WITH_${lib}=${TMPC_WITH_${lib}}")
endforeach ()

# TODO: Eigen3 is required now, because some older tests depend on it.
# This dependency must be removed in future.
find_package(Eigen3 3.2.10 REQUIRED)
find_package(blaze 3.7 REQUIRED)
find_package(LAPACK REQUIRED)

# Find boost.
find_package(Boost REQUIRED COMPONENTS graph exception)

# Find CasADi. Besides Used in some tests
find_package(CasADi 3.2.0)

if (TARGET casadi)
    message(STATUS "CASADI_INCLUDE_DIRS=${CASADI_INCLUDE_DIRS}")
else ()
    message(STATUS "CasADi not found; the tests depending on CasADi will not be built")
endif ()


# TMPC_WITH_CASADI
set(TMPC_WITH_CASADI ON CACHE BOOL "Build tmpc with CasADi generated function support classes")

if (TMPC_WITH_CASADI)
    if (TARGET casadi)
        include("CMakeModules/GenerateCasadiFunction.cmake")
    else ()
        message(SEND_ERROR, "You have selected TMPC_WITH_CASADI, but CasADi was not found")
    endif ()
endif ()


# TMPC_WITH_JSON
set(TMPC_WITH_JSON ON CACHE BOOL "Build tmpc with JSON support")

add_subdirectory(src)
add_subdirectory(examples)

# TMPC_WITH_TEST
set(TMPC_WITH_TEST ON CACHE BOOL "Build tmpc tests")

if (TMPC_WITH_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

# TMPC_WITH_BENCHMARK
set(TMPC_WITH_BENCHMARK ON CACHE BOOL "Build tmpc benchmarks")

if (TMPC_WITH_BENCHMARK)
    add_subdirectory(bench)
endif()
