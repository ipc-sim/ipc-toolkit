# geogram (https://github.com/BrunoLevy/geogram)
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

# Geogram builds its own sources with -Weverything on clang. Those are
# third-party warnings, so silence them.
if(TARGET geogram AND NOT MSVC)
    target_compile_options(geogram PRIVATE -w)
endif()

# Geogram adds its source dir joined with its relative build dir as a link
# directory, which does not exist when the two are unrelated (e.g., under CPM),
# and the linker warns about it. Nothing is linked from there.
if(TARGET geogram)
    set_property(TARGET geogram PROPERTY LINK_DIRECTORIES "")
endif()
