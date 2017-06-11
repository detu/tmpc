if("master" STREQUAL "")
  message(FATAL_ERROR "Tag for git checkout should not be empty.")
endif()

set(run 0)

if("/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-stamp/qpDUNES-gitinfo.txt" IS_NEWER_THAN "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-stamp/qpDUNES-gitclone-lastrun.txt")
  set(run 1)
endif()

if(NOT run)
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-stamp/qpDUNES-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES'")
endif()

# try the clone 3 times incase there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git" clone --origin "origin" "https://github.com/mkotlyar/qpDUNES.git" "qpDUNES"
    WORKING_DIRECTORY "/home/molivari/Work/software/projects/tmpc/build/extern/src"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/mkotlyar/qpDUNES.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git" checkout master
  WORKING_DIRECTORY "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'master'")
endif()

execute_process(
  COMMAND "/usr/bin/git" submodule init 
  WORKING_DIRECTORY "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to init submodules in: '/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES'")
endif()

execute_process(
  COMMAND "/usr/bin/git" submodule update --recursive 
  WORKING_DIRECTORY "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-stamp/qpDUNES-gitinfo.txt"
    "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-stamp/qpDUNES-gitclone-lastrun.txt"
  WORKING_DIRECTORY "/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/molivari/Work/software/projects/tmpc/build/extern/src/qpDUNES-stamp/qpDUNES-gitclone-lastrun.txt'")
endif()

