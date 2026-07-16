# TinyAD (https://github.com/patr-schm/tinyad)
# License: MIT
if(TARGET TinyAD::TinyAD)
    return()
endif()

message(STATUS "Third-party: creating target 'TinyAD::TinyAD'")

include(CPM)
CPMAddPackage(
    URI "gh:zfergus/tinyad#5bddce2c83e85d61b30e99359671037ee4a29349"
    OPTIONS "TINYAD_PARALLEL_BACKEND onetbb"
)

# Folder name for IDE
set_target_properties(TinyAD PROPERTIES FOLDER "ThirdParty")