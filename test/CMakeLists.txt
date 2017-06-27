#
# test
#
find_package(GTest REQUIRED)

find_package(Threads REQUIRED)    # required for gtest

set(TEST_DATA_PATH ${CMAKE_CURRENT_BINARY_DIR}/data)

# Generate code and test data for rk4 test
add_custom_command(
   OUTPUT
        pendulum_ode_generated.h
        pendulum_ode_generated.c
        ${TEST_DATA_PATH}/rk4/pendulum.txt
   COMMAND
       python
   ARGS
       ${CMAKE_CURRENT_SOURCE_DIR}/rk4_test_generate_data.py ${TEST_DATA_PATH}
   #[COMMAND command2 [ARGS] [args2...] ...]
   #[MAIN_DEPENDENCY depend]
   DEPENDS
       rk4_test_generate_data.py
   #[IMPLICIT_DEPENDS <lang1> depend1
   #                 [<lang2> depend2] ...]
   WORKING_DIRECTORY
       ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT
       "Generating model code and test data for rk4..."
)

# Generate code for casadi_generated_function test
add_custom_command(
   OUTPUT
        casadi_function_test_generated.h
        casadi_function_test_generated.c
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

# ~~~~~~~~~~ Our tests ~~~~~~~~~~~~~~
add_executable(tmpc_test
    #condensing_test.cpp
    rk4_test.cpp
    RealtimeIterationTest.cpp
    #mhe_test.cpp
    QpWorkspaceTest.cpp
    QpOasesWorkspaceTest.cpp
    QpWorkspaceSolveTest.cpp
    casadi_function_test.cpp
    gtest_tools_eigen.hpp

    ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.h
    ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.c
    ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.h
    ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.c
)

target_compile_definitions(tmpc_test
    PRIVATE TEST_DATA_PATH="${TEST_DATA_PATH}")

if (HPMPC_FOUND)
    #target_sources(tmpc_test
    #    PRIVATE hpmpc_test.cpp
    #)
endif()

target_include_directories(tmpc_test
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
)

target_link_libraries(tmpc_test
    tmpc
    GTest::GTest
    GTest::Main
    gmock
    ${CMAKE_THREAD_LIBS_INIT}    # Needed for gtest
)

add_subdirectory(matrix)

#enable_testing()
#add_test(AllTests tmpc_test)
