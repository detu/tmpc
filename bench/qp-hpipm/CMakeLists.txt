add_executable(bench-tmpc-qp-hpipm
    Main.cpp
    HpipmRiccatiBench.cpp
)


target_link_libraries(bench-tmpc-qp-hpipm
    tmpc-qp-hpipm
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
