if(TARGET TightInclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'TightInclusion'")

include(FetchContent)
FetchContent_Declare(
    Tight-Inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG 5f295da8ed79bc3b1165bb1de1b937c44105c3de
    GIT_SHALLOW FALSE
)

option(TIGHT_INCLUSION_WITH_NO_ZERO_TOI "Enable refinement if CCD produces a zero ToI" ON)

FetchContent_MakeAvailable(Tight-Inclusion)

add_library(TightInclusion ALIAS tight_inclusion)
