# Abseil (https://github.com/abseil/abseil-cpp)
# License: Apache 2.0
if(TARGET absl::flat_hash_map)
    return()
endif()

message(STATUS "Third-party: creating target 'absl::flat_hash_map'")

option(ABSL_PROPAGATE_CXX_STD "Use CMake C++ standard meta features (e.g. cxx_std_11) that propagate to targets that link to Abseil" ON)
option(ABSL_USE_SYSTEM_INCLUDES "Silence warnings in Abseil headers by marking them as SYSTEM includes" ON)
option(ABSL_BUILD_TESTING "If ON, Abseil will build all of Abseil's own tests." OFF)

include(CPM)
CPMAddPackage("gh:abseil/abseil-cpp#20230125.3")