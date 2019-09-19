#
# tmpc/src
#
# TODO: add all sources from the root CMakeLists.txt
#
set(TMPC_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/tmpc)

set (SRC
    ${TMPC_INCLUDE_DIR}/Matrix.hpp
    
	${TMPC_INCLUDE_DIR}/qp/CondensingN3.hpp
	${TMPC_INCLUDE_DIR}/qp/OcpQp.hpp

	${TMPC_INCLUDE_DIR}/qp/diagnostics.hpp
	${TMPC_INCLUDE_DIR}/qp/QpSolverException.hpp

    ocp/OcpGraph.cpp

    mpc/MpcOcpSize.cpp

    print/Exception.cpp
)

add_library(tmpc STATIC ${SRC})

target_include_directories(tmpc
    PUBLIC ${PROJECT_INCLUDE_DIR}
    PUBLIC ${EIGEN3_INCLUDE_DIR}
)


target_link_libraries(tmpc
    # TODO: this should be unnecessary since we link to Blaze, and Blaze already does this:
    # target_link_libraries(blaze INTERFACE $<BUILD_INTERFACE:${LAPACK_LIBRARIES}>)
    # But for some reason it does not work.
    # PUBLIC ${LAPACK_openblas_LIBRARY}
    # PUBLIC ${BLAS_mkl_core_LIBRARY}
    # PUBLIC ${BLAS_mkl_sequential_LIBRARY}
    PUBLIC "/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64_lin/libmkl_rt.so"
    PUBLIC Boost::graph
    PUBLIC Boost::exception
)

#
# Add CasADi interface
#
if (TMPC_WITH_CASADI)
    add_library(tmpc-casadi STATIC 
        casadi/GeneratedFunction.cpp
    )

    target_include_directories(tmpc-casadi
        PUBLIC ${CASADI_INCLUDE_DIRS}
    )

    target_link_libraries(tmpc-casadi 
        PUBLIC tmpc
    )
endif ()

#
# Add qpOASES interface
#
if (TMPC_WITH_qpOASES)
    add_subdirectory("qp-qpoases")
endif()


if (blaze_FOUND)
    target_link_libraries(tmpc
        PUBLIC blaze::blaze
    )
endif ()

#
# Add HPMPC interface
#
if (TMPC_WITH_HPMPC)
    add_subdirectory("qp-hpmpc")
endif()

#
# Add HPIPM interface
#
if (TMPC_WITH_HPIPM)
    add_subdirectory("qp-hpipm")
endif()

#
# Add treeQP interface
#
if (TMPC_WITH_treeQP)
    add_subdirectory("qp-treeqp")
endif()

#
# Build tmpc-blasfeo module
#
if (TMPC_WITH_BLASFEO)
    add_subdirectory("blasfeo")
endif()

#
# Link to nlohmann_json
#
if (TMPC_WITH_JSON)
    # Find Nlohmann json library
    find_package(nlohmann_json REQUIRED)

    add_library(tmpc-json INTERFACE)

    target_link_libraries(tmpc-json
        INTERFACE tmpc
        INTERFACE nlohmann_json
    )
endif ()


if (TMPC_WITH_treeQP)
    add_subdirectory(qptool)
endif ()