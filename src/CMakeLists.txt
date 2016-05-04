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
	${RTMC_INCLUDE_DIR}/qpDUNESProgram.hpp
	${RTMC_INCLUDE_DIR}/MultiStageQPSize.hpp
	${RTMC_INCLUDE_DIR}/qpOASESProgram.hpp
	${RTMC_INCLUDE_DIR}/CondensingSolver.hpp
	${RTMC_INCLUDE_DIR}/OptimalControlProblem.hpp
	${RTMC_INCLUDE_DIR}/CyberMotionOCP.hpp
	${RTMC_INCLUDE_DIR}/CableRobotOCP.hpp
	${RTMC_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp
	${RTMC_INCLUDE_DIR}/Condensing.hpp
)

set (SRC
	MotionPlatformX.cpp
	CyberMotion.cpp
	CyberMotion1D.cpp
	MultiStageQPSize.cpp
	CyberMotionOCP.cpp
	CableRobotOCP.cpp
	casadi_interface/CasADiGeneratedFunction.cpp
	cablerobot_generated.c
	cablerobot_generated.h
	cybermotion_generated.c
	cybermotion_generated.h
	qpOASESProgram.cpp
)

add_library(mpmc STATIC ${HDR} ${SRC})
#target_link_libraries(mpmc ${qpDUNES_STATIC_LIBRARIES})
target_link_libraries(mpmc ${qpOASES_STATIC_LIBRARIES})


