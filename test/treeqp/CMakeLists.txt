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

gtest_discover_tests(test-tmpc-qp-treeqp)
