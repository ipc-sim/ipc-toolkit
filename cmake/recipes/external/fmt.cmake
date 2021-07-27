if(TARGET fmt::fmt)
    return()
endif()

message(STATUS "Third-party: creating target 'fmt::fmt'")

include(FetchContent)
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 8.0.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(fmt)
