add_executable(create_qp
    create_qp.cpp
)

target_link_libraries(create_qp
    PRIVATE tmpc
)

add_executable(solve_qp
    solve_qp.cpp
)

target_link_libraries(solve_qp
    PRIVATE tmpc
)