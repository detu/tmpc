find_package(benchmark REQUIRED)
find_package(Threads REQUIRED)


#~~~~~~~~~~~~~~~~
add_executable(bench-tmpc
    math/Cholesky.cpp
    math/Trsv.cpp
    math/SyrkPotrf.cpp
    
    CondensingBench.cpp
    ClassicalRiccatiBench.cpp
    FactorizedRiccatiBench.cpp
    StaticFactorizedRiccatiBench.cpp
    OcpGraphBench.cpp
)


target_link_libraries(bench-tmpc
    tmpc
    benchmark
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
