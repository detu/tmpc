set(src
    CondensingTest.cpp
)

if (TMPC_WITH_qpOASES)
    list(APPEND src QpOasesWorkspaceTest.cpp)
endif ()

if (TMPC_WITH_HPMPC)
    list(APPEND src HpmpcWorkspaceTest.cpp)
endif ()

if (TMPC_WITH_HPIPM)
    list(APPEND src HpipmWorkspaceTest.cpp)
endif ()

message("src=${src}")

add_executable(tmpc_test_qp ${src})

target_link_libraries(tmpc_test_qp
    tmpc
    GTest::GTest
    GTest::Main
)

#GTEST_ADD_TESTS(tmpc_test_qp "" AUTO)
add_test(NAME tmpc_test_qp COMMAND tmpc_test_qp)