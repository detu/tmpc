find_package(BLASFEO REQUIRED)

add_library(tmpc-qp-hpmpc
    ${TMPC_INCLUDE_DIR}/qp/HpmpcWorkspace.hpp
    HpmpcWorkspace.cpp
)

target_include_directories(tmpc-qp-hpmpc
    PUBLIC ${HPMPC_INCLUDE_DIR}
    PUBLIC ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-qp-hpmpc
    PUBLIC tmpc
    PUBLIC ${HPMPC_STATIC_LIB}
    PUBLIC ${BLASFEO_STATIC_LIB}
)
