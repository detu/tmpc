# Generate code for casadi_generated_function test
add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.h
        ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.c
    COMMAND
        python
    ARGS
        ${CMAKE_CURRENT_SOURCE_DIR}/casadi_function_test_generate_code.py casadi_function_test_generated.c
    #[COMMAND command2 [ARGS] [args2...] ...]
    #[MAIN_DEPENDENCY depend]
    DEPENDS
        casadi_function_test_generate_code.py
    #[IMPLICIT_DEPENDS <lang1> depend1
    #                 [<lang2> depend2] ...]
    WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT
        "Generating code for casadi function test..."
)

add_executable(tmpc_test_casadi
    casadi_function_test.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.h
    ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.c
)

target_include_directories(tmpc_test_casadi
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(tmpc_test_casadi
    tmpc
    GTest::GTest
    GTest::Main
)