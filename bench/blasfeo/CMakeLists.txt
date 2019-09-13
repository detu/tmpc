add_executable(bench-tmpc-blasfeo
    Main.cpp
    BlasfeoMathBench.cpp
)


target_link_libraries(bench-tmpc-blasfeo
    tmpc-blasfeo
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
