#
# RTMC/src
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
    ${RTMC_INCLUDE_DIR}/core/RealtimeIteration.hpp
	${RTMC_INCLUDE_DIR}/core/OptimalControlProblem.hpp
    ${RTMC_INCLUDE_DIR}/core/ODEModel.hpp
    ${RTMC_INCLUDE_DIR}/core/Trajectory.hpp
        
	${RTMC_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp casadi_interface/CasADiGeneratedFunction.cpp
	
	${RTMC_INCLUDE_DIR}/integrator/RK4.hpp
	
	${RTMC_INCLUDE_DIR}/qp/Condensing.hpp
	${RTMC_INCLUDE_DIR}/qp/CondensingSolver.hpp  qp/CondensingSolver.cpp
	${RTMC_INCLUDE_DIR}/qp/qpDUNESProgram.hpp
	${RTMC_INCLUDE_DIR}/qp/qpOASESProgram.hpp    qp/qpOASESProgram.cpp
)

if (HPMPC_FOUND)
    set(SRC ${SRC} ${RTMC_INCLUDE_DIR}/qp/HPMPCSolver.hpp qp/HPMPCSolver.cpp)
endif()

add_library(mpmc STATIC ${HDR} ${SRC})
#target_link_libraries(mpmc ${qpDUNES_STATIC_LIBRARIES})
target_link_libraries(mpmc ${qpOASES_STATIC_LIBRARIES})


