#
# test
#

add_subdirectory(googletest)

# ~~~~~~~~~~ Our tests ~~~~~~~~~~~~~~
set(SRC 
    condensing_test.cpp
    rk4_test.cpp
    realtime_iteration_test.cpp
    qp_test.cpp
    qp_solver_test.cpp
    eigen_print_wrap.hpp
    
    pendulum_ode_generated.h    pendulum_ode_generated.c
   )
   
if (HPMPC_FOUND)
    set(SRC ${SRC} hpmpc_test.cpp)
endif()
   
add_executable(tmpc_test ${SRC})

target_link_libraries(tmpc_test 
    tmpc 
    gtest
    gtest_main
    )

if (qpOASES_FOUND)
    target_link_libraries(tmpc_test ${qpOASES_STATIC_LIBRARIES})
endif ()

if (HPMPC_FOUND)
    target_link_libraries(tmpc_test ${HPMPC_STATIC_LIB})
endif ()

#target_link_libraries(tmpc_test ${PROJECT_DEPENDENCY_LIBS})

#enable_testing()
#add_test(AllTests tmpc_test)
