# Tight Inclusion (https://github.com/Continuous-Collision-Detection/Tight-Inclusion)
# License: MIT
if(TARGET tight_inclusion::tight_inclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'tight_inclusion::tight_inclusion'")

include(CPM)
CPMAddPackage("gh:Continuous-Collision-Detection/Tight-Inclusion@1.0.6")

# Folder name for IDE
set_target_properties(tight_inclusion PROPERTIES FOLDER "ThirdParty")