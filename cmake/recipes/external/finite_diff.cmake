if(TARGET FiniteDiff::FiniteDiff)
    return()
endif()

message(STATUS "Third-party: creating target 'FiniteDiff::FiniteDiff'")

include(FetchContent)
FetchContent_Declare(
    finite-diff
    GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
    GIT_TAG f35375d2db00618d19ed7e4d2ce505288006f403
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(finite-diff)

add_library(FiniteDiff::FiniteDiff ALIAS FiniteDiff)
