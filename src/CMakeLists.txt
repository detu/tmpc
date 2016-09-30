#
# tmpc/src
#
# TODO: add all sources from the root CMakeLists.txt
#

set (SRC
    ${PROJECT_INCLUDE_DIR}/core/RealtimeIteration.hpp
    ${PROJECT_INCLUDE_DIR}/core/Trajectory.hpp
    ${PROJECT_INCLUDE_DIR}/core/gauss_newton.hpp
    ${PROJECT_INCLUDE_DIR}/core/matrix.hpp

	${PROJECT_INCLUDE_DIR}/casadi_interface/GeneratedFunction.hpp

	${PROJECT_INCLUDE_DIR}/integrator/rk4.hpp

	${PROJECT_INCLUDE_DIR}/qp/Condensing.hpp
	${PROJECT_INCLUDE_DIR}/qp/MultiStageQPSolution.hpp
	
	${PROJECT_INCLUDE_DIR}/qp/diagnostics.hpp
)

add_library(tmpc STATIC ${SRC})

target_include_directories(tmpc
    PUBLIC ${PROJECT_INCLUDE_DIR}
)

if (EIGEN3_FOUND)
    target_sources(tmpc
        PUBLIC ${PROJECT_INCLUDE_DIR}/kernel/eigen.hpp
    )
endif()

#
# Add qpDUNES interface 
#
if (qpDUNES_FOUND)
    add_dependencies(tmpc qpDUNES)
    
    target_sources(tmpc
        PRIVATE ${PROJECT_INCLUDE_DIR}/qp/qpDUNESProgram.hpp
	)
	
	target_include_directories(tmpc
        PRIVATE ${qpDUNES_INCLUDE_DIRS}
    )
    
    target_link_libraries(tmpc
        PUBLIC ${qpDUNES_LIBRARY_DIRS}
    )
endif()

#
# Add qpOASES interface 
#
if (qpOASES_FOUND)
    add_dependencies(tmpc qpOASES)

    target_sources(tmpc
        PRIVATE ${PROJECT_INCLUDE_DIR}/qp/CondensingSolver.hpp  
        PRIVATE qp/CondensingSolver.cpp
        PRIVATE ${PROJECT_INCLUDE_DIR}/qp/qpOASESProgram.hpp    
        PRIVATE qp/qpOASESProgram.cpp	
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
    add_dependencies(tmpc hpmpc)
    
    target_sources(tmpc
        PRIVATE ${PROJECT_INCLUDE_DIR}/qp/HPMPCSolver.hpp
        PRIVATE qp/HPMPCSolver.cpp
        PRIVATE ${PROJECT_INCLUDE_DIR}/qp/HPMPCProblem.hpp
        PRIVATE ${PROJECT_INCLUDE_DIR}/qp/HPMPCSolution.hpp
        PRIVATE qp/HPMPCProblemExport.hpp
        PRIVATE qp/HPMPCProblemExport.cpp
    )
    
    target_include_directories(tmpc
        PRIVATE ${HPMPC_INCLUDE_DIR}
    )
    
    target_link_libraries(tmpc
        PUBLIC ${HPMPC_STATIC_LIB}
    )
endif()


