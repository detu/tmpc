find_package(BLASFEO REQUIRED)

add_library(tmpc-hpipm STATIC
    Tree.cpp
)

target_include_directories(tmpc-hpipm
    PUBLIC ${HPIPM_INCLUDE_DIR}
    PUBLIC ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-hpipm
    PUBLIC tmpc
    PUBLIC ${HPIPM_STATIC_LIB}
    PUBLIC blazefeo-blasfeo
)
