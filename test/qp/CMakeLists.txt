if (WITH_qpOASES)
    target_sources(tmpc_test
        PRIVATE qp/QpOasesWorkspaceTest.cpp
    )
endif ()

if (WITH_HPMPC)
    target_sources(tmpc_test
        PRIVATE qp/HpmpcWorkspaceTest.cpp
    )
endif ()

if (WITH_HPIPM)
    target_sources(tmpc_test
        PRIVATE qp/HpipmWorkspaceTest.cpp
    )
endif ()

target_sources(tmpc_test
    PRIVATE qp/CondensingTest.cpp
    PRIVATE qp/QpTest.cpp
)
