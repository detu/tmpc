find_package(benchmark REQUIRED)
find_package(Threads REQUIRED)

#~~~~~~~~~~~~~~~~
add_executable(tmpc_bench
    CondensingBench.cpp
)

target_link_libraries(tmpc_bench
    tmpc
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
