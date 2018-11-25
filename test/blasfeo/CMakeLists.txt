add_executable(test-tmpc-blasfeo
    BlasfeoTest.cpp
)

target_link_libraries(test-tmpc-blasfeo
    tmpc-blasfeo
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-blasfeo)
