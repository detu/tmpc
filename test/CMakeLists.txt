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

target_link_libraries(mpc_test mpmc ${GTEST_BOTH_LIBRARIES})

if (UNIX)
	# On UNIX, link to pthread -- gtest requires it.
	target_link_libraries(mpc_test pthread)
endif (UNIX)

if (qpOASES_FOUND)
    target_link_libraries(mpc_test ${qpOASES_STATIC_LIBRARIES})
endif ()

if (HPMPC_FOUND)
    target_link_libraries(mpc_test ${HPMPC_STATIC_LIB})
endif ()

add_test(AllTests mpc_test)