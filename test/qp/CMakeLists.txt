add_executable(test-tmpc-qp 
    CondensingN3Test.cpp
    CondensingN2Test.cpp
    VectorPropertyMapTest.cpp
)

target_link_libraries(test-tmpc-qp
    tmpc
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-qp)