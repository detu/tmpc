add_library(tmpc-qp-qpoases
    ${TMPC_INCLUDE_DIR}/qp/QpOasesWorkspace.hpp
    QpOasesWorkspace.cpp
)

target_include_directories(tmpc-qp-qpoases
    PRIVATE ${qpOASES_INCLUDE_DIRS}
)

target_link_libraries(tmpc-qp-qpoases
    PUBLIC tmpc
    PUBLIC ${qpOASES_STATIC_LIBRARIES}
)