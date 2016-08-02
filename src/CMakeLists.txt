#
# tmpc/src
#

# Disable alignment of fixed-sized Eigen matrices.
# With alignment enabled, this error appears in the 32-bit version:
#
# ******************************************************************************************************************
# Expression: (reinterpret_cast<size_t>(array) & 0xf) == 0 && "this assertion is explained here:
# " "http://eigen.tuxfamily.org/dox-devel/group__TopicUnalignedArrayAssert.html" " **** READ THIS WEB PAGE !!! ****"
# ******************************************************************************************************************
#
# See this link for more details: http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
#add_definitions(-DEIGEN_DONT_ALIGN_STATICALLY)

set (SRC
    ${PROJECT_INCLUDE_DIR}/core/RealtimeIteration.hpp
    ${PROJECT_INCLUDE_DIR}/core/Trajectory.hpp
    ${PROJECT_INCLUDE_DIR}/core/gauss_newton.hpp
    ${PROJECT_INCLUDE_DIR}/core/matrix.hpp

	${PROJECT_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp casadi_interface/CasADiGeneratedFunction.cpp

	${PROJECT_INCLUDE_DIR}/integrator/RK4.hpp

	${PROJECT_INCLUDE_DIR}/qp/Condensing.hpp
	${PROJECT_INCLUDE_DIR}/qp/CondensingSolver.hpp  qp/CondensingSolver.cpp
	${PROJECT_INCLUDE_DIR}/qp/qpDUNESProgram.hpp
	${PROJECT_INCLUDE_DIR}/qp/qpOASESProgram.hpp    qp/qpOASESProgram.cpp
	${PROJECT_INCLUDE_DIR}/qp/MultiStageQPSolution.hpp
)

if (HPMPC_FOUND)
    set(SRC ${SRC}
        ${PROJECT_INCLUDE_DIR}/qp/HPMPCSolver.hpp     qp/HPMPCSolver.cpp
        ${PROJECT_INCLUDE_DIR}/qp/HPMPCProblem.hpp
        ${PROJECT_INCLUDE_DIR}/qp/HPMPCSolution.hpp
        qp/HPMPCProblemExport.hpp                  qp/HPMPCProblemExport.cpp
    )
endif()

add_library(tmpc STATIC ${SRC})
