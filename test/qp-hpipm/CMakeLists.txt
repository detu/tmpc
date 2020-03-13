add_executable(test-tmpc-qp-hpipm
    HpipmWorkspaceTest.cpp
    HpipmVsReferenceTest.cpp
)

target_include_directories(test-tmpc-qp-hpipm
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-qp-hpipm
    tmpc-qp-hpipm
    ${LAPACK_LIBRARIES}
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-qp-hpipm)
