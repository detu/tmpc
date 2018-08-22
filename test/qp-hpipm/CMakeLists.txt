add_executable(test-tmpc-qp-hpipm
    HpipmWorkspaceTest.cpp
)

target_include_directories(test-tmpc-qp-hpipm
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-qp-hpipm
    tmpc-qp-hpipm
    GTest::GTest
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_qp "" AUTO)
add_test(NAME test-tmpc-qp-hpipm COMMAND test-tmpc-qp-hpipm)
