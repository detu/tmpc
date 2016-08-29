#
# test
#

add_subdirectory(googletest)

# Generate code and test data for rk4 test
add_custom_command(
   OUTPUT
        pendulum_ode_generated.h
        pendulum_ode_generated.c
        data/rk4/pendulum.txt
   COMMAND 
       python 
   ARGS 
       ${CMAKE_CURRENT_SOURCE_DIR}/rk4_test_generate_data.py
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
set(SRC
    condensing_test.cpp
    rk4_test.cpp
    realtime_iteration_test.cpp
    qp_test.cpp
    qp_solver_test.cpp
    casadi_function_test.cpp
    gtest_tools_eigen.hpp

    ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.h    
    ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.c
    ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.h    
    ${CMAKE_CURRENT_BINARY_DIR}/casadi_function_test_generated.c
   )

if (HPMPC_FOUND)
    set(SRC ${SRC} hpmpc_test.cpp)
endif()

add_executable(tmpc_test ${SRC})

target_include_directories(tmpc_test
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(tmpc_test
    tmpc
    gtest
    gtest_main
    gmock
    )

if (qpOASES_FOUND)
    target_link_libraries(tmpc_test ${qpOASES_STATIC_LIBRARIES})
endif ()

if (HPMPC_FOUND)
    target_link_libraries(tmpc_test ${HPMPC_STATIC_LIB})
endif ()

#target_link_libraries(tmpc_test ${PROJECT_DEPENDENCY_LIBS})

#enable_testing()
#add_test(AllTests tmpc_test)
