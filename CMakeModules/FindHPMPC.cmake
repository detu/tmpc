MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for HPMPC package: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for HPMPC include directories" )

FIND_PATH(HPMPC_INCLUDE_DIR "c_interface.h" ${HPMPC_DIR}/include)
IF( HPMPC_INCLUDE_DIR )
	MESSAGE( STATUS "Found HPMPC include directories: ${HPMPC_INCLUDE_DIR} \n" )
	SET( HPMPC_INCLUDE_DIRS_FOUND TRUE )
ELSE( HPMPC_INCLUDE_DIR )
	MESSAGE( STATUS "Could not find HPMPC include directories \n" )
ENDIF( HPMPC_INCLUDE_DIR )

#
# Libraries
#
FIND_LIBRARY( HPMPC_STATIC_LIB hpmpc ${HPMPC_DIR})

IF( HPMPC_STATIC_LIB )
	MESSAGE( STATUS "Found HPMPC static library: ${HPMPC_STATIC_LIB} \n" )
	SET( HPMPC_STATIC_LIBS_FOUND TRUE )
ELSE( HPMPC_STATIC_LIB )
	MESSAGE( STATUS "Could not find HPMPC static library.\n" )
	SET( HPMPC_STATIC_LIBS_FOUND FALSE )
ENDIF( HPMPC_STATIC_LIB )

#
# And finally set found flag...
#
IF( HPMPC_INCLUDE_DIRS_FOUND AND HPMPC_STATIC_LIBS_FOUND )
	SET( HPMPC_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )
