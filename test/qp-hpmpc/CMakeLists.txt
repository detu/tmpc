add_executable(test-tmpc-qp-hpmpc
    HpmpcWorkspaceTest.cpp
)

target_include_directories(test-tmpc-qp-hpmpc
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-qp-hpmpc
    tmpc-qp-hpmpc
    GTest::GTest
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_qp "" AUTO)
add_test(NAME test-tmpc-qp-hpmpc COMMAND test-tmpc-qp-hpmpc)
