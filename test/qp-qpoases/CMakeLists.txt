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

#GTEST_ADD_TESTS(tmpc_test_qp "" AUTO)
add_test(NAME test-tmpc-qp-qpoases COMMAND test-tmpc-qp-qpoases)
