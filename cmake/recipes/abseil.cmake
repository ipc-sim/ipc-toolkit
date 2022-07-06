if(TARGET absl::flat_hash_map)
    return()
endif()

message(STATUS "Third-party: creating target 'absl::flat_hash_map'")

option(ABSL_PROPAGATE_CXX_STD "Use CMake C++ standard meta features (e.g. cxx_std_11) that propagate to targets that link to Abseil" ON)
option(ABSL_BUILD_TESTING "If ON, Abseil will build all of Abseil's own tests." OFF)

include(FetchContent)
FetchContent_Declare(
    abseil
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
    GIT_TAG 20220623.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(abseil)
