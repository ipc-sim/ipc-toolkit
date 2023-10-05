# Broadphase GPU (https://github.com/dbelgrod/broadphase-gpu)
# License: MIT
if(TARGET STQ::CPU)
    return()
endif()

message(STATUS "Third-party: creating target 'STQ::CPU'")

option(STQ_WITH_CPU  "Enable CPU Implementation"   ON)
option(STQ_WITH_CUDA "Enable CUDA Implementation" OFF) # get this through GPU_CCD

include(CPM)
CPMAddPackage("gh:dbelgrod/broadphase-gpu#48bfb9c4907ae07319574f4726cc451eb98055a3")