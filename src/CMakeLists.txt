#
# tmpc/src
#
# TODO: add all sources from the root CMakeLists.txt
#
set(TMPC_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/tmpc)

set (SRC    
    ocp/OcpTree.cpp
    ocp/DynamicOcpSize.cpp

    print/Exception.cpp
)

add_library(tmpc STATIC ${SRC})

target_include_directories(tmpc
    PUBLIC ${PROJECT_INCLUDE_DIR}
    PUBLIC ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(tmpc
    PUBLIC Boost::exception
    PUBLIC blaze::blaze
    PUBLIC blazefeo
)


#
# Add CasADi generated functions interface
#
if (TMPC_WITH_CASADI)
    add_subdirectory("casadi")
endif ()


#
# Add qpOASES interface
#
if (TMPC_WITH_qpOASES)
    add_subdirectory("qpoases")
endif()


#
# Add HPMPC interface
#
if (TMPC_WITH_HPMPC)
    add_subdirectory("hpmpc")
endif()

#
# Add HPIPM interface
#
if (TMPC_WITH_HPIPM)
    add_subdirectory("hpipm")
endif()

#
# Add treeQP interface
#
if (TMPC_WITH_treeQP)
    add_subdirectory("treeqp")
endif()

#
# Link to nlohmann_json
#
if (TMPC_WITH_JSON)
    add_subdirectory("json")
endif ()


if (TMPC_WITH_treeQP)
    add_subdirectory("qptool")
endif ()