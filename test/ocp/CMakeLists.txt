add_executable(tmpc_test_ocp
    SensitivityTupleTest.cpp
)

target_link_libraries(tmpc_test_ocp
    tmpc
    GTest::GTest
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_ocp "" AUTO)
add_test(NAME tmpc_test_ocp COMMAND tmpc_test_ocp)