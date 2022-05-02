if(TARGET absl::flat_hash_map)
    return()
endif()

message(STATUS "Third-party: creating target 'absl::flat_hash_map'")

include(FetchContent)
FetchContent_Declare(
    abseil
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
    GIT_TAG fad73ab077b846ff83e921ea1f5a6f831ea07ec2
)
FetchContent_MakeAvailable(abseil)
