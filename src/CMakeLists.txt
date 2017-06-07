#
# tmpc/src
#
# TODO: add all sources from the root CMakeLists.txt
#
set(TMPC_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/tmpc)

set (SRC
    ${TMPC_INCLUDE_DIR}/Matrix.hpp
    
    ${TMPC_INCLUDE_DIR}/util/TrajectoryPoint.hpp
    ${TMPC_INCLUDE_DIR}/util/gauss_newton.hpp
    ${TMPC_INCLUDE_DIR}/util/problem_specific.hpp

    ${TMPC_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp

	${TMPC_INCLUDE_DIR}/integrator/rk4.hpp

    ${TMPC_INCLUDE_DIR}/qp/QpSize.hpp    qp/QpSize.cpp
	${TMPC_INCLUDE_DIR}/qp/Condensing.hpp    qp/Condensing.cpp
	${TMPC_INCLUDE_DIR}/qp/MultiStageQPSolution.hpp
	${TMPC_INCLUDE_DIR}/qp/QuadraticProblem.hpp

	${TMPC_INCLUDE_DIR}/qp/diagnostics.hpp
	${TMPC_INCLUDE_DIR}/qp/UnsolvedQpException.hpp

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

if (BLAZE_DIR)
    target_include_directories(tmpc PUBLIC ${BLAZE_DIR})
    message(STATUS "BLAZE_DIR=${BLAZE_DIR}")
endif()

#
# Add qpOASES interface
#
if (qpOASES_FOUND)
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
if (HPMPC_FOUND)
    target_sources(tmpc PRIVATE 
        ${TMPC_INCLUDE_DIR}/qp/HPMPCSolver.hpp
        qp/HPMPCSolver.cpp
        ${TMPC_INCLUDE_DIR}/qp/HPMPCProblem.hpp
        ${TMPC_INCLUDE_DIR}/qp/HPMPCSolution.hpp
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
