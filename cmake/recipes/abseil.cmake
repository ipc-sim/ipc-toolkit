if(TARGET absl::flat_hash_map)
    return()
endif()

message(STATUS "Third-party: creating target 'absl::flat_hash_map'")

include(FetchContent)
FetchContent_Declare(
    abseil
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
    GIT_TAG 20220623.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(abseil)
