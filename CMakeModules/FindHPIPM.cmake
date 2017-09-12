MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for HPIPM package: \n" )
#
# Include folders
#
MESSAGE( STATUS "Looking for HPIPM include directories" )

FIND_PATH(HPIPM_INCLUDE_DIR "hpipm_target.h"
	HINTS ${HPIPM_DIR}/include $ENV{HPIPM_DIR}/include "/opt/hpipm/include"
)
IF( HPIPM_INCLUDE_DIR )
	MESSAGE( STATUS "Found HPIPM include directories: ${HPIPM_INCLUDE_DIR} \n" )
	SET( HPIPM_INCLUDE_DIRS_FOUND TRUE )
ELSE( HPIPM_INCLUDE_DIR )
	MESSAGE( STATUS "Could not find HPIPM include directories \n" )
ENDIF( HPIPM_INCLUDE_DIR )

#
# Libraries
#
FIND_LIBRARY( HPIPM_STATIC_LIB hpipm 
	HINTS ${HPIPM_DIR} $ENV{HPIPM_DIR} "/opt/hpipm/lib"
)

IF( HPIPM_STATIC_LIB )
	MESSAGE( STATUS "Found HPIPM static library: ${HPIPM_STATIC_LIB} \n" )
	SET( HPIPM_STATIC_LIBS_FOUND TRUE )
ELSE( HPIPM_STATIC_LIB )
	MESSAGE( STATUS "Could not find HPIPM static library.\n" )
	SET( HPIPM_STATIC_LIBS_FOUND FALSE )
ENDIF( HPIPM_STATIC_LIB )

#
# And finally set found flag...
#
IF( HPIPM_INCLUDE_DIRS_FOUND AND HPIPM_STATIC_LIBS_FOUND )
	SET( HPIPM_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )
