find_package(OpenMP REQUIRED)

# Generate code for casadi_generated_function test
tmpc_generate_casadi_function(SRCS HDRS "test_functions.py")

add_executable(tmpc_test_casadi
    GeneratedFunctionTest.cpp
    ${SRCS} ${HDRS}
)

target_include_directories(tmpc_test_casadi
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(tmpc_test_casadi
    tmpc
    GTest::GTest
    GTest::Main
    OpenMP::OpenMP_CXX
)

#GTEST_ADD_TESTS(tmpc_test_casadi "" AUTO)
add_test(NAME tmpc_test_casadi COMMAND tmpc_test_casadi)
