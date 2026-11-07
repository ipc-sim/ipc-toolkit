# Tracy (https://github.com/wolfpld/tracy)
# License: BSD-3-Clause
if(TARGET Tracy::TracyClient)
    return()
endif()

message(STATUS "Third-party: creating target 'Tracy::TracyClient'")

include(CPM)
CPMAddPackage(
    URI "gh:wolfpld/tracy@0.13.1"
    OPTIONS
      "TRACY_ENABLE ${IPC_TOOLKIT_WITH_TRACY}"
      "TRACY_ON_DEMAND ON"
)

set_target_properties(TracyClient PROPERTIES FOLDER "ThirdParty")