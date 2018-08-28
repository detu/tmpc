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

#GTEST_ADD_TESTS(test-tmpc-qp "" AUTO)
add_test(NAME test-tmpc-qp COMMAND test-tmpc-qp)