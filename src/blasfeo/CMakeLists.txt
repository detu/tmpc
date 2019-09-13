find_package(BLASFEO REQUIRED)

add_library(tmpc-blasfeo
    Blasfeo.cpp
)

target_include_directories(tmpc-blasfeo
    PUBLIC ${BLASFEO_INCLUDE_DIR}
)

target_link_libraries(tmpc-blasfeo
    PUBLIC tmpc
    PUBLIC ${BLASFEO_STATIC_LIB}
)
