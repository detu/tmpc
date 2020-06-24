find_package(BLASFEO REQUIRED)

add_library(tmpc-qp-hpipm INTERFACE)

target_include_directories(tmpc-qp-hpipm
    INTERFACE ${HPIPM_INCLUDE_DIR}
    INTERFACE ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-qp-hpipm
    INTERFACE tmpc
    INTERFACE ${HPIPM_STATIC_LIB}
    INTERFACE ${BLASFEO_STATIC_LIB}
)
