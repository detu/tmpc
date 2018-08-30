add_executable(test-tmpc-json
    JsonQpTest.cpp
)

target_link_libraries(test-tmpc-json PRIVATE 
    tmpc 
    GTest::GTest    
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_matrix "" AUTO)
add_test(NAME test-tmpc-json COMMAND test-tmpc-json)
