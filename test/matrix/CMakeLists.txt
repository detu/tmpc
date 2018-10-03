add_executable(test-tmpc-matrix
    #EigenAdaptorTest.cpp
    DynamicSubmatrixTest.cpp
    DynamicSubvectorTest.cpp
    KernelTest.cpp
    BlockMatrixViewTest.cpp
)

target_link_libraries(test-tmpc-matrix PRIVATE 
    tmpc 
    GTest::GTest    
    GTest::Main
)

gtest_discover_tests(test-tmpc-matrix)
