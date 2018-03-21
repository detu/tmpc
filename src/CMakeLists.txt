#
# tmpc/src
#
# TODO: add all sources from the root CMakeLists.txt
#
set(TMPC_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/tmpc)

set (SRC
    ${TMPC_INCLUDE_DIR}/Matrix.hpp
    ${TMPC_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp

	${TMPC_INCLUDE_DIR}/qp/Condensing.hpp    qp/Condensing.cpp
	${TMPC_INCLUDE_DIR}/qp/OcpQp.hpp

	${TMPC_INCLUDE_DIR}/qp/diagnostics.hpp
	${TMPC_INCLUDE_DIR}/qp/QpSolverException.hpp

    mpc/MpcOcpSize.cpp

    print/Exception.cpp
)

add_library(tmpc STATIC ${SRC})


# Enable position-independent code in tmpc.
# Otherwise, when trying to link it as a part of a shared library,
# a linker can complain.
set_property(TARGET tmpc PROPERTY POSITION_INDEPENDENT_CODE True)

target_include_directories(tmpc
    PUBLIC ${PROJECT_INCLUDE_DIR}
)

#
# Add CasADi interface
#
if (TMPC_WITH_CASADI_INTERFACE)
    target_include_directories(tmpc
        PUBLIC ${CASADI_INCLUDE_DIRS}
    )
endif ()

#
# Add qpOASES interface
#
if (TMPC_WITH_qpOASES)
    target_sources(tmpc
        PRIVATE ${TMPC_INCLUDE_DIR}/qp/QpOasesWorkspace.hpp
        PRIVATE qp/QpOasesWorkspace.cpp
	)

	target_include_directories(tmpc
        PRIVATE ${qpOASES_INCLUDE_DIRS}
    )

    target_link_libraries(tmpc
        PUBLIC ${qpOASES_STATIC_LIBRARIES}
    )
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
    target_sources(tmpc PRIVATE 
        ${TMPC_INCLUDE_DIR}/qp/HpmpcWorkspace.hpp
        qp/HpmpcWorkspace.cpp
    )

    target_include_directories(tmpc
        PRIVATE ${HPMPC_INCLUDE_DIR}
    )

    target_link_libraries(tmpc
        PUBLIC ${HPMPC_STATIC_LIB}
    )
endif()

#
# Add HPIPM interface
#
if (TMPC_WITH_HPIPM)
    target_sources(tmpc PRIVATE 
        ${TMPC_INCLUDE_DIR}/qp/HpipmWorkspace.hpp
        qp/HpipmWorkspace.cpp
    )

    target_include_directories(tmpc
        PUBLIC ${HPIPM_INCLUDE_DIR}
    )

    target_link_libraries(tmpc
        PUBLIC ${HPIPM_STATIC_LIB}
    )
endif()

#
# Link to BLASFEO
#
if (TMPC_WITH_BLASFEO)
    target_include_directories(tmpc
        PUBLIC ${BLASFEO_INCLUDE_DIR}
    )

    target_link_libraries(tmpc
        PUBLIC ${BLASFEO_STATIC_LIB}
    )
endif()

#
# Link to nlohmann_json
#
if (TMPC_WITH_JSON)
    target_link_libraries(tmpc
        PUBLIC nlohmann_json
    )
endif ()