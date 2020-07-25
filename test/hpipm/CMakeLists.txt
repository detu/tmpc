add_executable(test-tmpc-hpipm
    TreeTest.cpp
    NominalSolverTest.cpp
    TreeSolverTest.cpp
    ConvertTreeOcpQpTest.cpp
)

target_include_directories(test-tmpc-hpipm
    PRIVATE "../qp"
)

target_link_libraries(test-tmpc-hpipm
    tmpc-hpipm
    ${LAPACK_LIBRARIES}
    GTest::GTest
    GTest::Main
)

gtest_discover_tests(test-tmpc-hpipm)
