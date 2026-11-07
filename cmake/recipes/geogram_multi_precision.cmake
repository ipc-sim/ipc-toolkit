# Geogram's MultiPrecision Pluggable Software Module
# (https://github.com/BrunoLevy/geogram)
# License: BSD 3-Clause License
if(TARGET geogram::multi_precision)
    return()
endif()

message(STATUS "Third-party: creating target 'geogram::multi_precision'")

# The PSM is a standalone extract of geogram's exact arithmetic
# (expansion_nt), so we download its release archive instead of building all
# of geogram.
set(GEOGRAM_MULTI_PRECISION_VERSION "1.10.1")

include(CPM)
CPMAddPackage(
    NAME geogram_multi_precision
    URL "https://github.com/BrunoLevy/geogram/releases/download/v${GEOGRAM_MULTI_PRECISION_VERSION}/MultiPrecision_psm_${GEOGRAM_MULTI_PRECISION_VERSION}.zip"
    URL_HASH SHA256=1c4fa0ef6422547f08b11891a0f934b90a2bc91359dcfaeb89a490b3878c0e18
    DOWNLOAD_ONLY YES
)

# The archive also ships an optional CGAL example that carries an LGPL-3.0
# header. We never build it, so we delete it to keep the build tree free of
# non-permissive files.
file(REMOVE_RECURSE "${geogram_multi_precision_SOURCE_DIR}/CGAL")

# Downstream projects (e.g., PolyFEM) may link their own full copy of geogram,
# which defines the same GEO:: symbols as the PSM. Linking both would either
# fail with duplicate symbols (static) or silently mix two geogram versions
# (shared). Instead, we compile a copy of the PSM with its namespace renamed
# from GEO to ipc_geogram, so the two can never collide. We write the copy to
# the build tree so a shared CPM source cache is never modified.
set(GEOGRAM_MULTI_PRECISION_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/geogram_multi_precision")
foreach(file IN ITEMS MultiPrecision_psm.h MultiPrecision_psm.cpp)
    file(READ "${geogram_multi_precision_SOURCE_DIR}/${file}" content)
    # GEO only appears as the namespace name ("namespace GEO", "GEO::", and
    # "using namespace GEO"). Matching it as a whole word leaves macros such
    # as GEO_OS_WINDOWS untouched.
    string(REGEX REPLACE "([^A-Za-z0-9_])GEO([^A-Za-z0-9_])"
        "\\1ipc_geogram\\2" content "${content}")
    # If a future PSM spells the namespace in a way the pattern above misses,
    # the copy would silently keep GEO:: symbols and could clash again, so we
    # fail loudly instead.
    string(REGEX MATCH "(^|[^A-Za-z0-9_])GEO([^A-Za-z0-9_]|$)"
        leftover "${content}")
    if(leftover)
        message(FATAL_ERROR
            "Failed to rename namespace GEO to ipc_geogram in ${file}.")
    endif()
    file(WRITE "${GEOGRAM_MULTI_PRECISION_DIR}/tmp/${file}" "${content}")
    # Copy only when the content changed, so we don't trigger rebuilds.
    configure_file("${GEOGRAM_MULTI_PRECISION_DIR}/tmp/${file}"
        "${GEOGRAM_MULTI_PRECISION_DIR}/${file}" COPYONLY)
endforeach()

add_library(geogram_multi_precision STATIC
    "${GEOGRAM_MULTI_PRECISION_DIR}/MultiPrecision_psm.cpp"
)
add_library(geogram::multi_precision ALIAS geogram_multi_precision)

target_include_directories(geogram_multi_precision SYSTEM PUBLIC
    "$<BUILD_INTERFACE:${GEOGRAM_MULTI_PRECISION_DIR}>"
)

# The PSM header assumes geogram is a shared library unless told otherwise. On
# MSVC that marks every symbol __declspec(dllimport), which breaks linking
# against our static build, so we define GEO_STATIC_LIBS for the library and
# everything that includes its header.
target_compile_definitions(geogram_multi_precision PUBLIC GEO_STATIC_LIBS)

target_compile_features(geogram_multi_precision PUBLIC cxx_std_17)
set_target_properties(geogram_multi_precision PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# The expansion arithmetic relies on error-free transforms, which break if the
# compiler fuses a*b+c into an FMA, so we turn off contraction. We also silence
# the third-party warnings.
if(MSVC)
    target_compile_options(geogram_multi_precision PRIVATE "/fp:precise" "/w")
else()
    target_compile_options(geogram_multi_precision PRIVATE
        "-ffp-contract=off" "-w"
    )
endif()

# Folder name for IDE
set_target_properties(geogram_multi_precision PROPERTIES FOLDER "ThirdParty")