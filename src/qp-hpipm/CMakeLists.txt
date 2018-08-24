find_package(BLASFEO REQUIRED)

add_library(tmpc-qp-hpipm
    ${TMPC_INCLUDE_DIR}/qp/HpipmWorkspace.hpp
    HpipmWorkspace.cpp
)

target_include_directories(tmpc-qp-hpipm
    PUBLIC ${HPIPM_INCLUDE_DIR}
    PUBLIC ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-qp-hpipm
    PUBLIC tmpc
    PUBLIC ${HPIPM_STATIC_LIB}
    PUBLIC ${BLASFEO_STATIC_LIB}
)
