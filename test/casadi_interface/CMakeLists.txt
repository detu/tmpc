find_package(OpenMP REQUIRED)

# Generate code for casadi_generated_function test
tmpc_generate_casadi_function(SRCS HDRS "test_functions.py")

add_executable(test-tmpc-casadi
    GeneratedFunctionTest.cpp
    ${SRCS} ${HDRS}
)

target_include_directories(test-tmpc-casadi
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(test-tmpc-casadi
    tmpc
    GTest::GTest
    GTest::Main
    OpenMP::OpenMP_CXX
)

gtest_discover_tests(test-tmpc-casadi)
