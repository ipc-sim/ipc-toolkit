if(TARGET tight_inclusion::tight_inclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'tight_inclusion::tight_inclusion'")

option(TIGHT_INCLUSION_WITH_NO_ZERO_TOI "Enable refinement if CCD produces a zero ToI" ON)

include(CPM)
CPMAddPackage("gh:Continuous-Collision-Detection/Tight-Inclusion@1.0.3")