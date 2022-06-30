if(TARGET tight_inclusion::tight_inclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'tight_inclusion::tight_inclusion'")

option(TIGHT_INCLUSION_WITH_NO_ZERO_TOI "Enable refinement if CCD produces a zero ToI" ON)

if(EXISTS "${IPC_TOOLKIT_TIGHT_INCLUSION_PATH}")
    message(STATUS "Using Tight Inclusion found at: ${IPC_TOOLKIT_TIGHT_INCLUSION_PATH}")
    add_subdirectory("${IPC_TOOLKIT_TIGHT_INCLUSION_PATH}" "${PROJECT_BINARY_DIR}/tight-inclusion")
else()
    include(FetchContent)
    FetchContent_Declare(
        tight_inclusion
        GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
        GIT_TAG v1.0.2
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(tight_inclusion)
endif()
