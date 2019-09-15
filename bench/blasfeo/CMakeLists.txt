add_executable(bench-tmpc-blasfeo
    Main.cpp
    BlasfeoMath.cpp
    SyrkPotrf.cpp
)


target_link_libraries(bench-tmpc-blasfeo
    tmpc-blasfeo
    benchmark
    ${CMAKE_THREAD_LIBS_INIT}
)
