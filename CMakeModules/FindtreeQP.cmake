MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for treeQP package:" )
#
# Include folders
#
MESSAGE( STATUS "Looking for treeQP include directories" )

FIND_PATH(treeQP_INCLUDE_DIR "treeqp/src/tree_qp_common.h"
	HINTS "${treeQP_DIR}" "$ENV{treeQP_DIR}"				
)

IF (treeQP_INCLUDE_DIR)
	SET( treeQP_INCLUDE_DIRS_FOUND TRUE )
	# list(APPEND treeQP_INCLUDE_DIR "${treeQP_INCLUDE_DIR}/external")
	MESSAGE( STATUS "Found treeQP include directories: ${treeQP_INCLUDE_DIR}" )
ELSE ()
	MESSAGE( STATUS "Could not find treeQP include directories" )
ENDIF ()


set (treeQP_INCLUDE_DIRS ${treeQP_INCLUDE_DIR})

foreach (LIB "blasfeo" "qpoases" "hpmpc")
	list (APPEND treeQP_INCLUDE_DIRS "${treeQP_INCLUDE_DIR}/external/${LIB}/include")
endforeach ()


#
# Libraries
#
set (treeQP_STATIC_LIBS)

foreach (LIB "treeqp" "blasfeo" "qpoases" "hpmpc")
	FIND_LIBRARY( treeQP_${LIB}_STATIC_LIB ${LIB} 
		HINTS "${treeQP_DIR}/lib" "$ENV{treeQP_DIR}/lib" "/opt/treeqp/lib"
	)

	IF (treeQP_${LIB}_STATIC_LIB)
		MESSAGE( STATUS "Found treeQP ${LIB} static library: ${treeQP_${LIB}_STATIC_LIB}" )
		SET( treeQP_${LIB}_STATIC_LIB_FOUND TRUE )
		list (APPEND treeQP_STATIC_LIBS ${treeQP_${LIB}_STATIC_LIB})
		#list (APPEND treeQP_INCLUDE_DIRS "${treeQP_INCLUDE_DIR}/external/")
	ELSE ()
		MESSAGE( STATUS "Could not find treeQP ${LIB} static library." )
		SET( treeQP_${LIB}_STATIC_LIB_FOUND FALSE )
	ENDIF ()
endforeach ()


#
# And finally set found flag...
#
IF (treeQP_INCLUDE_DIRS_FOUND AND 
	treeQP_treeqp_STATIC_LIB_FOUND AND
	treeQP_blasfeo_STATIC_LIB_FOUND AND
	treeQP_qpoases_STATIC_LIB_FOUND AND
	treeQP_hpmpc_STATIC_LIB_FOUND)
	set(treeQP_FOUND TRUE)
	message(STATUS "treeQP_STATIC_LIBS=${treeQP_STATIC_LIBS}")
	message(STATUS "treeQP_INCLUDE_DIRS=${treeQP_INCLUDE_DIRS}")
ELSE ()
	set(treeQP_FOUND FALSE)
ENDIF ()

MESSAGE( STATUS "********************************************************************************" )
