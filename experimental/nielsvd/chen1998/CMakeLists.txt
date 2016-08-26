add_executable(example_chen1998 chen1998.cpp chen1998_ode_generated.c)
target_link_libraries(example_chen1998 tmpc ${qpOASES_STATIC_LIBRARIES})
