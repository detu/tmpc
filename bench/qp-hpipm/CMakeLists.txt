add_executable(bench-tmpc-hpipm
    Main.cpp
    Riccati.cpp
)


target_link_libraries(bench-tmpc-hpipm
    tmpc-hpipm
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
