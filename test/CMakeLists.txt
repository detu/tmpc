#
# RTMC/test
#

enable_testing()
find_package(GTest REQUIRED)

set(SRC 
    condensing_test.cpp
    rk4_test.cpp
    pendulum_ode_generated.h    pendulum_ode_generated.c
   )
   
if (HPMPC_FOUND)
    set(SRC ${SRC} hpmpc_test.cpp)
endif()
   
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(mpc_test ${SRC})

set (LIB mpmc ${GTEST_BOTH_LIBRARIES})
if (UNIX)
	# On UNIX, link to pthread -- gtest requires it.
	set(LIB ${LIB} pthread)
endif (UNIX)

target_link_libraries(mpc_test ${LIB})

add_test(AllTests mpc_test)