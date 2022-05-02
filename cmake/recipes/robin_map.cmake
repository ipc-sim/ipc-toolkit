if(TARGET tsl::robin_map)
    return()
endif()

message(STATUS "Third-party: creating target 'tsl::robin_map'")

include(FetchContent)
FetchContent_Declare(
    robin-map
    GIT_REPOSITORY https://github.com/Tessil/robin-map.git
    GIT_TAG v1.0.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(robin-map)
