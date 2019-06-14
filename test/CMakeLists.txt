#
# test
#
find_package(GTest 1.9.0 REQUIRED)
find_package(Threads REQUIRED)    # required for gtest
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Generate and build models
foreach (model pendulum)
    set(generated_name "generated_${model}")
    add_custom_command(
        OUTPUT
            ${generated_name}.h
            ${generated_name}.c
        COMMAND
            ${Python3_EXECUTABLE}
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

# ~~~~~~~~~~ Our tests ~~~~~~~~~~~~~~
add_executable(test-tmpc
    #mhe_test.cpp
    ocp/SensitivityTupleTest.cpp
    ocp/CondensedOcpSizeTest.cpp
    ocp/OcpGraphTest.cpp
    ocp/OcpSizePropertiesTest.cpp

    qp/CondensingN3Test.cpp
    qp/CondensingN2Test.cpp
    qp/VectorPropertyMapTest.cpp
    qp/MpipmWorkspaceTest.cpp
    qp/RiccatiTest.cpp

    integrator/ImplicitRungeKuttaTest.cpp

    matrix/MatrixTest.cpp
    matrix/DynamicSubmatrixTest.cpp
    matrix/DynamicSubvectorTest.cpp
    matrix/KernelTest.cpp
    matrix/BlockMatrixViewTest.cpp

    math/ExpmTest.cpp

    numeric/CentralDifferenceDerivativeTest.cpp
    numeric/NewtonSolverTest.cpp

    system/LtiContinuousToDiscreteTest.cpp

    random/MultivariateNormalDistributionTest.cpp

    estimation/KalmanFilterTest.cpp

    util/UnwrapTest.cpp
)

if (TMPC_WITH_CASADI_INTERFACE)
    # Generate code and test data for rk4 test
    add_custom_command(
        OUTPUT
            pendulum_ode_generated.h
            pendulum_ode_generated.c
            ${TEST_DATA_PATH}/rk4/pendulum.txt
        COMMAND
            ${Python3_EXECUTABLE}
        ARGS
            ${CMAKE_CURRENT_SOURCE_DIR}/integrator/rk4_test_generate_data.py ${TEST_DATA_PATH}
        #[COMMAND command2 [ARGS] [args2...] ...]
        #[MAIN_DEPENDENCY depend]
        DEPENDS
            integrator/rk4_test_generate_data.py
        #[IMPLICIT_DEPENDS <lang1> depend1
        #                 [<lang2> depend2] ...]
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT
            "Generating model code and test data for ExplicitRungeKutta4..."
    )

    target_sources(test-tmpc
        PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.h
        ${CMAKE_CURRENT_BINARY_DIR}/pendulum_ode_generated.c
        
        # TODO: Integrator tests depend on CasADi.
        # Remove this dependency.
        integrator/PendulumData.cpp
        integrator/PendulumOdeTest.cpp
        integrator/ExplicitRungeKutta4Test.cpp
        integrator/MultiStepIntegratorTest.cpp
    )

    target_link_libraries(test-tmpc
        PUBLIC tmpc-casadi
    )
endif ()

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

if (TARGET casadi)
    target_sources(test-tmpc PRIVATE
        integrator/PendulumTest.cpp
    )

    target_link_libraries(test-tmpc PRIVATE casadi)
endif ()


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