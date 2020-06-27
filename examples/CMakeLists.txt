add_subdirectory(qp)
add_subdirectory(numeric)

if (TMPC_WITH_JSON)
    add_subdirectory(json)
endif ()


if (TMPC_WITH_CASADI)
    add_subdirectory(casadi)
endif ()