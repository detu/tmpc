add_executable(tmpc_test_matrix
    #EigenAdaptorTest.cpp
    DynamicSubmatrixTest.cpp
    DynamicSubvectorTest.cpp
    KernelTest.cpp
)

target_link_libraries(tmpc_test_matrix PRIVATE 
    tmpc 
    GTest::GTest    
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_matrix "" AUTO)
add_test(NAME tmpc_test_matrix COMMAND tmpc_test_matrix)
