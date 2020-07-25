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


# ~~~~~~~~~~ Our tests ~~~~~~~~~~~~~~
add_executable(test-tmpc
    ElementGroupTest.cpp

    # property_map/VectorPropertyMapTest.cpp
    
    ocp/CondensedSizeTest.cpp
    ocp/OcpTreeTest.cpp
    ocp/DynamicOcpSizeTest.cpp
    ocp/StaticOcpSizeTest.cpp
    ocp/StaticOcpSolutionTest.cpp
    ocp/DynamicOcpKktValueTest.cpp

    # qp/CondensingN3Test.cpp
    # qp/CondensingN2Test.cpp
    qp/MpipmSolverTest.cpp
    qp/StaticOcpQpTest.cpp
    qp/DynamicOcpQpTest.cpp
    qp/DynamicRiccatiTest.cpp
    qp/StaticRiccatiTest.cpp

    integrator/ErkMassSpringDamperTest.cpp
    integrator/ErkDecayTest.cpp
    integrator/IrkTest.cpp
    integrator/IrkSimpleLinearOdeTest.cpp
    integrator/IrkSimpleTimeDependentOdeTest.cpp
    integrator/IrkSimpleDaeTest.cpp
    integrator/IrkMassSpringDamperTest.cpp
    integrator/IrkDecayTest.cpp

    matrix/MatrixTest.cpp
    matrix/DynamicSubmatrixTest.cpp
    matrix/KernelTest.cpp
    matrix/BlockMatrixViewTest.cpp

    math/ExpmTest.cpp
    math/RankTest.cpp
    math/LlhTest.cpp
    math/SyrkPotrfTest.cpp
    math/TrsvTest.cpp
    math/GgesTest.cpp

    numeric/CentralDifferenceDerivativeTest.cpp
    numeric/NewtonSolverTest.cpp

    lti/LtiContinuousToDiscreteTest.cpp
    lti/DareTest.cpp

    random/MultivariateNormalDistributionTest.cpp

    estimation/KalmanFilterTest.cpp

    util/UnwrapTest.cpp
)


target_link_libraries(test-tmpc PRIVATE
    tmpc
    GTest::GTest
    GTest::Main
    gmock
    ${LAPACK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}    # Needed for gtest
)


if (TARGET casadi)
    target_sources(test-tmpc PRIVATE
        integrator/IrkPendulumTest.cpp
    )

    target_link_libraries(test-tmpc PRIVATE casadi)
endif ()


if (TMPC_WITH_CASADI)
    add_subdirectory(casadi)
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
    add_subdirectory("hpipm")
endif ()
