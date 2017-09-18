#
# tmpc/src
#
# TODO: add all sources from the root CMakeLists.txt
#
set(TMPC_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/tmpc)

set (SRC
    ${TMPC_INCLUDE_DIR}/Matrix.hpp
    
    ${TMPC_INCLUDE_DIR}/util/problem_specific.hpp

    ${TMPC_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp

	${TMPC_INCLUDE_DIR}/integrator/rk4.hpp

    ${TMPC_INCLUDE_DIR}/qp/QpSize.hpp    qp/QpSize.cpp
	${TMPC_INCLUDE_DIR}/qp/Condensing.hpp    qp/Condensing.cpp
	${TMPC_INCLUDE_DIR}/qp/QuadraticProblem.hpp

	${TMPC_INCLUDE_DIR}/qp/diagnostics.hpp
	${TMPC_INCLUDE_DIR}/qp/QpSolverException.hpp

    ${TMPC_INCLUDE_DIR}/mpc/MpcQpSize.hpp mpc/MpcQpSize.cpp
    ${TMPC_INCLUDE_DIR}/mpc/MpcRealtimeIteration.hpp
    ${TMPC_INCLUDE_DIR}/mpc/MpcTrajectoryPoint.hpp
)

add_library(tmpc STATIC ${SRC})

# Enable modern C++ for tmpc.
set_property(TARGET tmpc PROPERTY CXX_STANDARD 14)

# Don't use compiler's C++ extensions
set_property(TARGET tmpc PROPERTY CXX_EXTENSIONS False)

# Enable position-independent code in tmpc.
# Otherwise, when trying to link it as a part of a shared library,
# a linker can complain.
set_property(TARGET tmpc PROPERTY POSITION_INDEPENDENT_CODE True)

target_include_directories(tmpc
    PUBLIC ${PROJECT_INCLUDE_DIR}
)

#
# Add qpOASES interface
#
if (WITH_qpOASES)
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

#
# Add HPMPC interface
#
if (WITH_HPMPC)
    target_sources(tmpc PRIVATE 
        ${TMPC_INCLUDE_DIR}/qp/HpmpcWorkspace.hpp
        qp/HpmpcWorkspace.cpp
        qp/HPMPCProblemExport.hpp
        qp/HPMPCProblemExport.cpp
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
if (WITH_HPIPM)
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
if (WITH_BLASFEO)
    target_include_directories(tmpc
        PUBLIC ${BLASFEO_INCLUDE_DIR}
    )

    target_link_libraries(tmpc
        PUBLIC ${BLASFEO_STATIC_LIB}
    )
endif()