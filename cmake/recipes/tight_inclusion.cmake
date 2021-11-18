if(TARGET tight_inclusion::tight_inclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'tight_inclusion::tight_inclusion'")

include(FetchContent)
FetchContent_Declare(
    tight_inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG 1dc56b9d9c202a381fe6bd998cee38484938a77c
    GIT_SHALLOW FALSE
)

option(TIGHT_INCLUSION_WITH_NO_ZERO_TOI "Enable refinement if CCD produces a zero ToI" ON)

FetchContent_MakeAvailable(tight_inclusion)

add_library(tight_inclusion::tight_inclusion ALIAS tight_inclusion)
