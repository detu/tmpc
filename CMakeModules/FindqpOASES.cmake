MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for qpOASES package: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for qpOASES include directories" )
FIND_PATH(qpOASES_INCLUDE_DIRS "qpOASES.hpp")
IF( qpOASES_INCLUDE_DIRS )
	MESSAGE( STATUS "Found qpOASES include directories: ${qpOASES_INCLUDE_DIRS} \n" )
	SET( qpOASES_INCLUDE_DIRS_FOUND TRUE )
ELSE( qpOASES_INCLUDE_DIRS )
	MESSAGE( STATUS "Could not find qpOASES include directories \n" )
ENDIF( qpOASES_INCLUDE_DIRS )

#
# Libraries
#
FIND_LIBRARY( qpOASES_STATIC_LIBRARIES
	NAMES qpOASES
)
IF( qpOASES_STATIC_LIBRARIES )
	MESSAGE( STATUS "Found qpOASES static library: ${qpOASES_STATIC_LIBRARIES} \n" )
	SET( qpOASES_STATIC_LIBS_FOUND TRUE )
ELSE( qpOASES_STATIC_LIBRARIES )
	MESSAGE( STATUS "Could not find qpOASES static library.\n" )
	SET( qpOASES_STATIC_LIBS_FOUND FALSE )
ENDIF( qpOASES_STATIC_LIBRARIES )

#
# And finally set found flag...
#
IF( qpOASES_INCLUDE_DIRS_FOUND AND qpOASES_STATIC_LIBS_FOUND )
	SET( qpOASES_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )
