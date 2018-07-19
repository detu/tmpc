find_package(benchmark REQUIRED)
find_package(Threads REQUIRED)


#~~~~~~~~~~~~~~~~
add_executable(tmpc_bench
    CondensingBench.cpp
    MathBench.cpp
)


if (TMPC_WITH_CASADI_INTERFACE)
    tmpc_generate_casadi_function(GENERATED_SRCS GENERATED_HDRS "casadi_math_bench.py")
    
    target_sources(tmpc_bench
        PRIVATE CasadiMathBench.cpp ${GENERATED_SRCS} ${GENERATED_HDRS}
    )

    target_include_directories(tmpc_bench
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
    )
endif ()


if (TMPC_WITH_BLASFEO)
    target_sources(tmpc_bench
        PRIVATE BlasfeoMathBench.cpp
    )
endif ()


target_link_libraries(tmpc_bench
    tmpc
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
