add_executable(bench-tmpc-blasfeo
    Main.cpp
    Trmm.cpp
    Gemm.cpp
    GemmNoWrapper.cpp
    SyrkPotrf.cpp
)


target_link_libraries(bench-tmpc-blasfeo
    tmpc-blasfeo
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
