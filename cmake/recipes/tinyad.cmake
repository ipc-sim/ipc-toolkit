# finite-diff (https://github.com/patr-schm/tinyad)
# License: MIT
if(TARGET TinyAD)
    return()
endif()

message(STATUS "Third-party: creating target 'TinyAD::TinyAD'")

include(CPM)
CPMAddPackage("gh:patr-schm/tinyad#4b48d1a1a588874556a692a3abbdecd0db4c23e1")

# Folder name for IDE
set_target_properties(TinyAD PROPERTIES FOLDER "ThirdParty")