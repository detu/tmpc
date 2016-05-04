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

set(HDR 
	${RTMC_INCLUDE_DIR}/MotionPlatform.hpp
	${RTMC_INCLUDE_DIR}/MotionPlatformX.hpp
	${RTMC_INCLUDE_DIR}/CyberMotion.hpp
	${RTMC_INCLUDE_DIR}/CyberMotion1D.hpp
	${RTMC_INCLUDE_DIR}/ModelPredictiveController.hpp
	${RTMC_INCLUDE_DIR}/OptimalControlProblem.hpp
	${RTMC_INCLUDE_DIR}/CyberMotionOCP.hpp
	${RTMC_INCLUDE_DIR}/CableRobotOCP.hpp
)

set (SRC
	MotionPlatformX.cpp
	CyberMotion.cpp
	CyberMotion1D.cpp
	
	CyberMotionOCP.cpp
	CableRobotOCP.cpp
	
	cablerobot_generated.c
	cablerobot_generated.h
	cybermotion_generated.c
	cybermotion_generated.h
	
	${RTMC_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp casadi_interface/CasADiGeneratedFunction.cpp
	
	${RTMC_INCLUDE_DIR}/qp/Condensing.hpp
	${RTMC_INCLUDE_DIR}/qp/CondensingSolver.hpp
	${RTMC_INCLUDE_DIR}/qp/MultiStageQPSize.hpp    qp/MultiStageQPSize.cpp
	${RTMC_INCLUDE_DIR}/qp/qpDUNESProgram.hpp
	${RTMC_INCLUDE_DIR}/qp/qpOASESProgram.hpp    qp/qpOASESProgram.cpp
)

add_library(mpmc STATIC ${HDR} ${SRC})
#target_link_libraries(mpmc ${qpDUNES_STATIC_LIBRARIES})
target_link_libraries(mpmc ${qpOASES_STATIC_LIBRARIES})


