cmake_minimum_required (VERSION 3.8)
project (tmpc)

# Enable modern C++
set(CMAKE_CXX_STANDARD 17)

# Don't use compiler's C++ extensions
set(CMAKE_CXX_EXTENSIONS False)

# Include CMakeToolsHelpers for easier debugging from VS Code
include(CMakeToolsHelpers OPTIONAL)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

# List of QP solvers
set(qp_solvers "qpOASES" "HPMPC" "HPIPM")

# List of selectable extern libs
set(extern_libs ${qp_solvers} "BLASFEO")

# Create TMPC_WITH_${lib} cached variables and find selected packages
foreach (lib ${extern_libs})
    set(TMPC_WITH_${lib} OFF CACHE BOOL "Link to ${lib}")

    if (TMPC_WITH_${lib})
        find_package(${lib} REQUIRED)
    endif ()

    message(STATUS "TMPC_WITH_${lib}=${TMPC_WITH_${lib}}")
endforeach ()

if (TMPC_WITH_HPIPM AND NOT TMPC_WITH_BLASFEO)
    message(SEND_ERROR "If you choose TMPC_WITH_HPIPM, you have to choose TMPC_WITH_BLASFEO as well.")
endif ()

find_package(Eigen3 3.2.10)

if (EIGEN3_FOUND)
    message(STATUS "Eigen3 found in ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen3 NOT FOUND")
endif()


find_package(blaze)


# TMPC_WITH_CASADI_INTERFACE
set(TMPC_WITH_CASADI_INTERFACE ON CACHE BOOL "Build tmpc with CasADi generated function support classes")

if (TMPC_WITH_CASADI_INTERFACE)
    find_package(CasADi 3.2.0 REQUIRED)
    message(STATUS "CASADI_INCLUDE_DIRS=${CASADI_INCLUDE_DIRS}")

    include ("CMakeModules/GenerateCasadiFunction.cmake")
endif ()


# TMPC_WITH_JSON
set(TMPC_WITH_JSON ON CACHE BOOL "Build tmpc with JSON support")

if (TMPC_WITH_JSON)
    # Find Nlohmann json library
    find_package(nlohmann_json REQUIRED)
endif ()

add_subdirectory(src)
add_subdirectory(examples)

# TMPC_WITH_TEST
set(TMPC_WITH_TEST ON CACHE BOOL "Build tmpc tests")

if (TMPC_WITH_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

# TMPC_WITH_BENCH
set(TMPC_WITH_BENCH ON CACHE BOOL "Build tmpc benchmarks")

if (TMPC_WITH_BENCH)
    add_subdirectory(bench)
endif()
