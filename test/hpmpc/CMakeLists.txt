add_executable(test-tmpc-qp-hpmpc
    HpmpcWorkspaceTest.cpp
)

target_include_directories(test-tmpc-qp-hpmpc
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-qp-hpmpc
    tmpc-qp-hpmpc
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-qp-hpmpc)
