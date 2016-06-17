MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for HPMPC package: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for HPMPC include directories" )

FIND_PATH(HPMPC_INCLUDE_DIRS "c_interface.h" ${HPMPC_DIR}/include)
IF( HPMPC_INCLUDE_DIRS )
	MESSAGE( STATUS "Found HPMPC include directories: ${HPMPC_INCLUDE_DIRS} \n" )
	SET( HPMPC_INCLUDE_DIRS_FOUND TRUE )
ELSE( HPMPC_INCLUDE_DIRS )
	MESSAGE( STATUS "Could not find HPMPC include directories \n" )
ENDIF( HPMPC_INCLUDE_DIRS )

#
# Libraries
#
FIND_LIBRARY( HPMPC_STATIC_LIBRARIES hpmpc ${HPMPC_DIR})

IF( HPMPC_STATIC_LIBRARIES )
	MESSAGE( STATUS "Found HPMPC static library: ${HPMPC_STATIC_LIBRARIES} \n" )
	SET( HPMPC_STATIC_LIBS_FOUND TRUE )
ELSE( HPMPC_STATIC_LIBRARIES )
	MESSAGE( STATUS "Could not find HPMPC static library.\n" )
	SET( HPMPC_STATIC_LIBS_FOUND FALSE )
ENDIF( HPMPC_STATIC_LIBRARIES )

#
# And finally set found flag...
#
IF( HPMPC_INCLUDE_DIRS_FOUND AND HPMPC_STATIC_LIBS_FOUND )
	SET( HPMPC_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )
