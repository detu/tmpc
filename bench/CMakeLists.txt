find_package(benchmark REQUIRED)
find_package(Threads REQUIRED)


#~~~~~~~~~~~~~~~~
add_executable(bench-tmpc
    CondensingBench.cpp
    math/ZeroMatrixAssign.cpp
    math/Cholesky.cpp
    math/Syrk.cpp
    math/Trmm.cpp
    math/Trmv.cpp
    math/Trsv.cpp
    math/Gemm.cpp
    math/GemmSubmatrix.cpp
    math/SqrNorm.cpp
    math/LowerMatrixScalarMultiply.cpp
    math/SyrkPotrf.cpp
    math/Column.cpp
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


add_subdirectory("blas")

if (TMPC_WITH_CASADI)
    add_subdirectory("casadi")
endif ()


if (TMPC_WITH_BLASFEO)
    add_subdirectory("blasfeo")
endif ()


if (TMPC_WITH_HPIPM)
    add_subdirectory("qp-hpipm")
endif ()
