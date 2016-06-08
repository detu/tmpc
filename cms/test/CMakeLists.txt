#
# RTMC/test
#

enable_testing()
find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})
add_executable(mpc_test 
    mpc_test.cpp 
    condensing_test.cpp
    cybermotion_ocp_test.cpp
    cybermotion_generated_test.cpp
    rk4_test.cpp
    pendulum_ode_generated.h    pendulum_ode_generated.c
)
#add_executable(mpc_test condensing_test.cpp)

set (LIB mpmc ${GTEST_BOTH_LIBRARIES})
if (UNIX)
	# On UNIX, link to pthread -- gtest requires it.
	set(LIB ${LIB} pthread)
endif (UNIX)

target_link_libraries(mpc_test ${LIB})

add_executable(perftest perftest.cpp)
target_link_libraries(perftest mpmc)

add_test(AllTests mpc_test)