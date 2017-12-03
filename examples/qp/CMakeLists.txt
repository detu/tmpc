add_executable(create_qp create_qp.cpp)
target_link_libraries(create_qp PRIVATE tmpc)

add_executable(solve_qp solve_qp.cpp)
target_link_libraries(solve_qp PRIVATE tmpc)

if (TMPC_WITH_HPIPM)
    add_executable(soft_constraints soft_constraints.cpp)
    target_link_libraries(soft_constraints PRIVATE tmpc)
endif ()