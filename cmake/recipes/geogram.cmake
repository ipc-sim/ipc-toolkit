# finite-diff (https://github.com/BrunoLevy/geogram)
# License: BSD 3-Clause License
if(TARGET geogram::geogram)
    return()
endif()

message(STATUS "Third-party: creating target 'geogram::geogram'")

include(CPM)
CPMAddPackage(
    URI "gh:BrunoLevy/geogram@1.10.1"
    OPTIONS
        "GEOGRAM_WITH_GRAPHICS OFF"  
        "GEOGRAM_WITH_LEGACY_NUMERICS OFF"  
        "GEOGRAM_WITH_HLBFGS OFF"  
        "GEOGRAM_WITH_TETGEN OFF"  
        "GEOGRAM_WITH_TRIANGLE OFF"  
        "GEOGRAM_WITH_LUA OFF"  
        "GEOGRAM_LIB_ONLY ON"
)

if(NOT TARGET geogram::geogram AND TARGET geogram)  
    add_library(geogram::geogram ALIAS geogram)  
endif()