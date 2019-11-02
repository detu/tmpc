add_executable(test-tmpc-blasfeo
    CustomMatrixTest.cpp
    DynamicMatrixTest.cpp
    SyrkPotrfTest.cpp
    PotrfTest.cpp
    GemmTest.cpp
)

target_link_libraries(test-tmpc-blasfeo
    tmpc-blasfeo
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-blasfeo)
