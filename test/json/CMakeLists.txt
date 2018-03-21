add_executable(tmpc_test_json
    JsonQpTest.cpp
)

target_link_libraries(tmpc_test_json PRIVATE 
    tmpc 
    GTest::GTest    
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_matrix "" AUTO)
add_test(NAME tmpc_test_json COMMAND tmpc_test_json)
