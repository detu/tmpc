find_package(casadi REQUIRED 3.5.0)

message(STATUS "CASADI_INCLUDE_DIRS=${CASADI_INCLUDE_DIRS}")

include("cmake/GenerateCasadiFunction.cmake")

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
