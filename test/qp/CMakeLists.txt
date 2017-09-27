if (TMPC_WITH_qpOASES)
    target_sources(tmpc_test
        PRIVATE qp/QpOasesWorkspaceTest.cpp
    )
endif ()

if (TMPC_WITH_HPMPC)
    target_sources(tmpc_test
        PRIVATE qp/HpmpcWorkspaceTest.cpp
    )
endif ()

if (TMPC_WITH_HPIPM)
    target_sources(tmpc_test
        PRIVATE qp/HpipmWorkspaceTest.cpp
    )
endif ()

target_sources(tmpc_test
    PRIVATE qp/CondensingTest.cpp
)
