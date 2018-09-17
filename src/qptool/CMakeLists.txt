if (TMPC_WITH_JSON AND TMPC_WITH_treeQP)
    add_executable(qptool qptool.cpp)
    target_link_libraries(qptool 
        PRIVATE tmpc-qp-treeqp
        PRIVATE tmpc-json
    )
endif ()