add_executable(test-tmpc-qp-treeqp
    DualNewtonTreeWorkspaceTest.cpp
)

target_include_directories(test-tmpc-qp-treeqp
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-qp-treeqp
    tmpc-qp-treeqp
    GTest::GTest
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_qp "" AUTO)
add_test(NAME test-tmpc-qp-treeqp COMMAND test-tmpc-qp-treeqp)
