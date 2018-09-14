add_executable(test-tmpc-json
    JsonQpTest.cpp
)

target_link_libraries(test-tmpc-json PRIVATE 
    tmpc 
    GTest::GTest    
    GTest::Main
)

gtest_discover_tests(test-tmpc-json)
