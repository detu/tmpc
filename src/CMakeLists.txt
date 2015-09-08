#
# RTMC/src
#

enable_testing()
find_package(GTest REQUIRED)

include(CompilerOptions.cmake)

# Disable alignment of fixed-sized Eigen matrices.
# With alignment enabled, this error appears in the 32-bit version:
#
# ******************************************************************************************************************
# Expression: (reinterpret_cast<size_t>(array) & 0xf) == 0 && "this assertion is explained here: 
# " "http://eigen.tuxfamily.org/dox-devel/group__TopicUnalignedArrayAssert.html" " **** READ THIS WEB PAGE !!! ****"
# ******************************************************************************************************************
#
# See this link for more details: http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
add_definitions(-DEIGEN_DONT_ALIGN_STATICALLY)

set(HDR 
	${RTMC_INCLUDE_DIR}/MotionPlatform.hpp
	${RTMC_INCLUDE_DIR}/MotionPlatformX.hpp
	${RTMC_INCLUDE_DIR}/CyberMotion.hpp
	${RTMC_INCLUDE_DIR}/MPC_Controller.hpp
	${RTMC_INCLUDE_DIR}/MotionPlatformModelPredictiveController.hpp
	${RTMC_INCLUDE_DIR}/MultiStageQP.hpp
	${RTMC_INCLUDE_DIR}/MultiStageQPSize.hpp
	${RTMC_INCLUDE_DIR}/QuadraticProgram.hpp
	${RTMC_INCLUDE_DIR}/CondensingSolver.hpp
)

set (SRC
	MotionPlatformX.cpp
	CyberMotion.cpp
	MPC_Controller.cpp
	MotionPlatformModelPredictiveController.cpp
	CMS_output_jacobian.cpp
	MultiStageQP.cpp
	MultiStageQPSize.cpp
	CondensingSolver.cpp
)

add_library(mpmc STATIC ${HDR} ${SRC})
#target_link_libraries(mpmc ${qpDUNES_STATIC_LIBRARIES})
target_link_libraries(mpmc ${qpOASES_LIBRARIES})

add_executable(perftest perftest.cpp)
target_link_libraries(perftest mpmc)

include_directories(${GTEST_INCLUDE_DIRS})
add_executable(mpc_test mpc_test.cpp condensing_test.cpp)
#add_executable(mpc_test condensing_test.cpp)

set (LIB mpmc ${GTEST_BOTH_LIBRARIES} ${qpOASES_LIBRARIES})
if (UNIX)
	# On UNIX, link to pthread -- gtest requires it.
	set(LIB ${LIB} pthread)
endif (UNIX)

target_link_libraries(mpc_test ${LIB})

add_test(AllTests mpc_test)