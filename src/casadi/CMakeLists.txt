find_package(casadi REQUIRED 3.5.0)

add_library(tmpc-casadi STATIC 
    GeneratedFunction.cpp
)

target_include_directories(tmpc-casadi
    PUBLIC ${CASADI_INCLUDE_DIRS}
)

target_link_libraries(tmpc-casadi 
    PUBLIC tmpc
    PUBLIC casadi
)