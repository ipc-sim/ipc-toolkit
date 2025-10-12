# TinyAD (https://github.com/patr-schm/tinyad)
# License: MIT
if(TARGET TinyAD::TinyAD)
    return()
endif()

message(STATUS "Third-party: creating target 'TinyAD::TinyAD'")

include(Eigen) # TinyAD will find Eigen if Eigen::Eigen3 is not defined

include(CPM)
CPMAddPackage("gh:patr-schm/tinyad#4b48d1a1a588874556a692a3abbdecd0db4c23e1")

add_library(TinyAD::TinyAD ALIAS TinyAD)

# Folder name for IDE
set_target_properties(TinyAD PROPERTIES FOLDER "ThirdParty")