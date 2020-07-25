# find_package(qpOASES_e REQUIRED)
# find_package(BLASFEO REQUIRED)

add_library(tmpc-qp-treeqp
    ${TMPC_INCLUDE_DIR}/qp/DualNewtonTreeWorkspace.hpp
    DualNewtonTreeWorkspace.cpp
)

target_include_directories(tmpc-qp-treeqp
    PUBLIC ${treeQP_INCLUDE_DIRS}
    #PUBLIC ${qpOASES_e_INCLUDE_DIR}
    #PUBLIC ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-qp-treeqp
    PUBLIC tmpc
    PUBLIC ${treeQP_STATIC_LIBS}
    #PUBLIC ${treeQP_STATIC_LIB}
    #PUBLIC ${qpOASES_e_STATIC_LIBRARIES}
    #PUBLIC ${BLASFEO_STATIC_LIB}
)
