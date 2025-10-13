# TinyAD (https://github.com/patr-schm/tinyad)
# License: MIT
if(TARGET TinyAD::TinyAD)
    return()
endif()

message(STATUS "Third-party: creating target 'TinyAD::TinyAD'")

include(CPM)
CPMAddPackage(
    URI "gh:zfergus/tinyad#957ede5bcf4b8d924f27230083508c664a9f7a7d"
    OPTIONS "TINYAD_PARALLEL_BACKEND onetbb"
)

# Folder name for IDE
set_target_properties(TinyAD PROPERTIES FOLDER "ThirdParty")