if(TARGET tsl::robin_map)
    return()
endif()

message(STATUS "Third-party: creating target 'tsl::robin_map'")

include(CPM)
CPMAddPackage("gh:Tessil/robin-map@1.2.1")