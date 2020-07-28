find_package(benchmark REQUIRED)
find_package(Threads REQUIRED)


#~~~~~~~~~~~~~~~~
add_executable(bench-tmpc
    Main.cpp
    
    math/Cholesky.cpp
    math/Trsv.cpp
    math/SyrkPotrf.cpp
    
    # CondensingBench.cpp
    qp/DynamicClassicalRiccatiBench.cpp
    qp/DynamicFactorizedRiccatiBench.cpp
    qp/StaticClassicalRiccatiBench.cpp
    qp/StaticFactorizedRiccatiBench.cpp
    OcpTreeBench.cpp
)


target_link_libraries(bench-tmpc
    tmpc
    benchmark
    blazefeo-blasfeo
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)


if (TMPC_WITH_CASADI)
    add_subdirectory("casadi")
endif ()


if (TMPC_WITH_HPIPM)
    add_subdirectory("qp-hpipm")
endif ()
