add_executable(test-tmpc-qp-qpoases
    QpOasesWorkspaceTest.cpp
)

target_include_directories(test-tmpc-qp-qpoases
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-qp-qpoases
    tmpc-qp-qpoases
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-qp-qpoases)
