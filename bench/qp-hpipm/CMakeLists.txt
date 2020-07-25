add_executable(bench-tmpc-hpipm
    Main.cpp
    HpipmRiccatiBench.cpp
)


target_link_libraries(bench-tmpc-hpipm
    tmpc-hpipm
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
