MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for treeQP package: \n" )
#
# Include folders
#
MESSAGE( STATUS "Looking for treeQP include directories" )

FIND_PATH(treeQP_INCLUDE_DIR "treeqp/src/tree_ocp_qp_common.h"
	HINTS "${treeQP_DIR}" "$ENV{treeQP_DIR}"				
)
IF( treeQP_INCLUDE_DIR )
	SET( treeQP_INCLUDE_DIRS_FOUND TRUE )
	list(APPEND treeQP_INCLUDE_DIR "${treeQP_INCLUDE_DIR}/external")
	MESSAGE( STATUS "Found treeQP include directories: ${treeQP_INCLUDE_DIR} \n" )
ELSE( treeQP_INCLUDE_DIR )
	MESSAGE( STATUS "Could not find treeQP include directories \n" )
ENDIF( treeQP_INCLUDE_DIR )

#
# Libraries
#
FIND_LIBRARY( treeQP_STATIC_LIB treeqp 
	HINTS "${treeQP_DIR}/lib" "$ENV{treeQP_DIR}/lib" "/opt/treeqp/lib"
)

IF( treeQP_STATIC_LIB )
	MESSAGE( STATUS "Found treeQP static library: ${treeQP_STATIC_LIB} \n" )
	SET( treeQP_STATIC_LIBS_FOUND TRUE )
ELSE( treeQP_STATIC_LIB )
	MESSAGE( STATUS "Could not find treeQP static library.\n" )
	SET( treeQP_STATIC_LIBS_FOUND FALSE )
ENDIF( treeQP_STATIC_LIB )

#
# And finally set found flag...
#
IF( treeQP_INCLUDE_DIRS_FOUND AND treeQP_STATIC_LIBS_FOUND )
	SET( treeQP_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )
