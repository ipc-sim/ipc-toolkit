if(TARGET TightInclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'TightInclusion'")

include(FetchContent)
FetchContent_Declare(
    Tight-Inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG f08d717f824b16a297d9b44513aa3807180aca79
    GIT_SHALLOW TRUE
)

option(TIGHT_INCLUSION_WITH_NO_ZERO_TOI "Enable refinement if CCD produces a zero ToI" ON)

FetchContent_MakeAvailable(Tight-Inclusion)

add_library(TightInclusion ALIAS tight_inclusion)
