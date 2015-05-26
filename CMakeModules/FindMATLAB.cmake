# - this module looks for Matlab
# Defines:
#  MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MATLAB_LIBRARIES:   required libraries: libmex, etc
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY:  path to libmx.lib
#  MATLAB_ENG_LIBRARY: path to libeng.lib

#=============================================================================
# Copyright 2005-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

MESSAGE(STATUS "CMAKE_GENERATOR=" ${CMAKE_GENERATOR})

set(MATLAB_FOUND 0)
if(WIN32)
	#set(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.14;MATLABROOT]")
	if(${CMAKE_GENERATOR} MATCHES "Visual Studio 10 Win64" 
	    OR ${CMAKE_GENERATOR} MATCHES "Visual Studio 11 Win64"
		OR ${CMAKE_GENERATOR} MATCHES "Visual Studio 12 2013 Win64"
	    OR ${CMAKE_GENERATOR} MATCHES "Visual Studio 11 2012 Win64")
		set(MATLAB_LIB_DIR ${MATLAB_ROOT}/extern/lib/win64/microsoft)
	elseif(${CMAKE_GENERATOR} MATCHES "Visual Studio 10" OR ${CMAKE_GENERATOR} MATCHES "Visual Studio 11")
		set(MATLAB_LIB_DIR ${MATLAB_ROOT}/extern/lib/win32/microsoft)
	else()
		if(MATLAB_FIND_REQUIRED)
			message(FATAL_ERROR "Generator not compatible: ${CMAKE_GENERATOR}")
		endif()
    endif()
	
	MESSAGE(STATUS "MATLAB_LIB_DIR=" ${MATLAB_LIB_DIR})
  find_library(MATLAB_MEX_LIBRARY
    libmex
    ${MATLAB_LIB_DIR}
    )
  find_library(MATLAB_MX_LIBRARY
    libmx
    ${MATLAB_LIB_DIR}
    )
  find_library(MATLAB_ENG_LIBRARY
    libeng
    ${MATLAB_LIB_DIR}
    )

	find_path(MATLAB_INCLUDE_DIR
		"mex.h"
		${MATLAB_ROOT}/extern/include
    )
	
	find_path(SIMULINK_INCLUDE_DIR
		"simulink.h"
		${MATLAB_ROOT}/simulink/include
	)
else()
  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Regular x86
    set(MATLAB_LIB_DIR
      /usr/local/matlab-7sp1/bin/glnx86/
      /opt/matlab-7sp1/bin/glnx86/
      $ENV{HOME}/matlab-7sp1/bin/glnx86/
      $ENV{HOME}/redhat-matlab/bin/glnx86/
      )
  else()
    # AMD64:
    set(MATLAB_LIB_DIR
      /usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/glnxa64/
      /usr/local/MATLAB/R2014a/bin/glnxa64/
      /usr/local/matlab-7sp1/bin/glnxa64/
      /opt/matlab-7sp1/bin/glnxa64/
      $ENV{HOME}/matlab7_64/bin/glnxa64/
      $ENV{HOME}/matlab-7sp1/bin/glnxa64/
      $ENV{HOME}/redhat-matlab/bin/glnxa64/
      )
  endif()

  find_library(MATLAB_MEX_LIBRARY
    mex
    ${MATLAB_LIB_DIR}
    )
  find_library(MATLAB_MX_LIBRARY
    mx
    ${MATLAB_LIB_DIR}
    )
  find_library(MATLAB_ENG_LIBRARY
    eng
    ${MATLAB_LIB_DIR}
    )
  find_path(MATLAB_INCLUDE_DIR
    "mex.h"
    "/usr/local/MATLAB/MATLAB_Production_Server/R2013a/extern/include/"
    "/usr/local/MATLAB/R2014a/extern/include/"
    "/usr/local/matlab-7sp1/extern/include/"
    "/opt/matlab-7sp1/extern/include/"
    "$ENV{HOME}/matlab-7sp1/extern/include/"
    "$ENV{HOME}/redhat-matlab/extern/include/"
    )

endif()

# This is common to UNIX and Win32:
set(MATLAB_LIBRARIES
  ${MATLAB_MEX_LIBRARY}
  ${MATLAB_MX_LIBRARY}
  ${MATLAB_ENG_LIBRARY}
)

if(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)
  set(MATLAB_FOUND 1)
endif()

mark_as_advanced(
  MATLAB_LIBRARIES
  MATLAB_MEX_LIBRARY
  MATLAB_MX_LIBRARY
  MATLAB_ENG_LIBRARY
  MATLAB_INCLUDE_DIR
  MATLAB_FOUND
  MATLAB_LIB_DIR
)

