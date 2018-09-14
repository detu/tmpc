#
# test
#
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)    # required for gtest
find_package(PythonInterp REQUIRED)

# Generate and build models
foreach (model pendulum)
    set(generated_name "generated_${model}")
    add_custom_command(
        OUTPUT
            ${generated_name}.h
            ${generated_name}.c
        COMMAND
            ${PYTHON_EXECUTABLE}
        ARGS
            ${CMAKE_CURRENT_SOURCE_DIR}/model/${model}.py ${generated_name}.c
        DEPENDS
            model/${model}.py
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT
            "Generating model code and test data for model ${model}..."
    )

    add_library(model_${model} STATIC
        model/${model}.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/${generated_name}.h
        ${CMAKE_CURRENT_BINARY_DIR}/${generated_name}.c
    )

    target_include_directories(model_${model}
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
    )

    target_link_libraries(model_${model} PRIVATE tmpc)
endforeach ()

set(TEST_DATA_PATH ${CMAKE_CURRENT_BINARY_DIR}/data)

# Generate code and test data for rk4 test
add_custom_command(
    OUTPUT
        pendulum_ode_generated.h
        pendulum_ode_generated.c
        ${TEST_DATA_PATH}/rk4/pendulum.txt
    COMMAND
        ${PYTHON_EXECUTABLE}
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

# ~~~~~~~~~~ Our tests ~~~~~~~~~~~~~~
add_executable(test-tmpc
    #rk4_test.cpp
    #mhe_test.cpp
    ocp/SensitivityTupleTest.cpp
    ocp/CondensedOcpSizeTest.cpp
    ocp/OcpGraphTest.cpp
    ocp/OcpSizePropertiesTest.cpp

    ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.h
    ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.c
)

target_compile_definitions(test-tmpc
    PRIVATE TEST_DATA_PATH="${TEST_DATA_PATH}")

#if (TMPC_WITH_QPOASES)
#    target_sources(test-tmpc
#        PRIVATE RealtimeIterationTest.cpp
#    )
#endif()

target_include_directories(test-tmpc
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
)

target_link_libraries(test-tmpc PRIVATE
    tmpc
    GTest::GTest
    GTest::Main
    gmock
    ${CMAKE_THREAD_LIBS_INIT}    # Needed for gtest
)

add_subdirectory(matrix)
add_subdirectory(qp)
add_subdirectory(integrator)

if (TMPC_WITH_CASADI_INTERFACE)
    add_subdirectory(casadi_interface)
endif ()

if (TMPC_WITH_JSON)
    add_subdirectory(json)
endif ()

gtest_discover_tests(test-tmpc)


if (TMPC_WITH_qpOASES)
    add_subdirectory("qp-qpoases")
endif ()


if (TMPC_WITH_treeQP)
    add_subdirectory("qp-treeqp")
endif ()


if (TMPC_WITH_HPMPC)
    add_subdirectory("qp-hpmpc")
endif ()


if (TMPC_WITH_HPIPM)
    add_subdirectory("qp-hpipm")
endif ()