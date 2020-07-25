find_package(BLASFEO REQUIRED)

add_library(tmpc-qp-hpmpc INTERFACE)

target_include_directories(tmpc-qp-hpmpc
    INTERFACE ${HPMPC_INCLUDE_DIR}
    INTERFACE ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-qp-hpmpc
    INTERFACE tmpc
    INTERFACE ${HPMPC_STATIC_LIB}
    INTERFACE ${BLASFEO_STATIC_LIB}
)
