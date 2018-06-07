MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for qpOASES_e package: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for qpOASES_e include directories" )
FIND_PATH(qpOASES_e_INCLUDE_DIRS "qpOASES_e.h"
	HINTS ${qpOASES_e_DIR}/include $ENV{qpOASES_e_DIR}/include
)
IF( qpOASES_e_INCLUDE_DIRS )
	MESSAGE( STATUS "Found qpOASES_e include directories: ${qpOASES_e_INCLUDE_DIRS} \n" )
	SET( qpOASES_e_INCLUDE_DIRS_FOUND TRUE )
ELSE( qpOASES_e_INCLUDE_DIRS )
	MESSAGE( STATUS "Could not find qpOASES_e include directories \n" )
ENDIF( qpOASES_e_INCLUDE_DIRS )

#
# Libraries
#
FIND_LIBRARY( qpOASES_e_STATIC_LIBRARIES
	NAMES libqpOASES_e.a qpOASES_e
	HINTS ${qpOASES_e_DIR}/lib ${qpOASES_e_DIR}/bin $ENV{qpOASES_e_DIR}/build/libs $ENV{qpOASES_e_DIR}/bin
)
IF( qpOASES_e_STATIC_LIBRARIES )
	MESSAGE( STATUS "Found qpOASES_e static library: ${qpOASES_e_STATIC_LIBRARIES} \n" )
	SET( qpOASES_e_STATIC_LIBS_FOUND TRUE )
ELSE( qpOASES_e_STATIC_LIBRARIES )
	MESSAGE( STATUS "Could not find qpOASES_e static library.\n" )
	SET( qpOASES_e_STATIC_LIBS_FOUND FALSE )
ENDIF( qpOASES_e_STATIC_LIBRARIES )

#
# And finally set found flag...
#
IF( qpOASES_e_INCLUDE_DIRS_FOUND AND qpOASES_e_STATIC_LIBS_FOUND )
	SET( qpOASES_e_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )
