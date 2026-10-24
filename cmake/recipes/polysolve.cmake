# PolySolve (https://github.com/polyfem/polysolve)
# License: MIT
if(TARGET polysolve::polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'polysolve::polysolve'")

# Only the Eigen built-in linear solvers are needed for the demo simulator;
# disable the optional heavyweight backends and their dependency tails.
# (Accelerate also conflicts with Eigen's BLAS prototypes on recent macOS SDKs.)
set(POLYSOLVE_WITH_ACCELERATE OFF CACHE BOOL "Enable Apple Accelerate" FORCE)
set(POLYSOLVE_WITH_AMGCL      OFF CACHE BOOL "Use AMGCL"               FORCE)
set(POLYSOLVE_WITH_CHOLMOD    OFF CACHE BOOL "Enable Cholmod library"  FORCE)
set(POLYSOLVE_WITH_HYPRE      OFF CACHE BOOL "Enable hypre"            FORCE)
set(POLYSOLVE_WITH_MKL        OFF CACHE BOOL "Enable MKL library"      FORCE)
set(POLYSOLVE_WITH_PARDISO    OFF CACHE BOOL "Enable Pardiso library"  FORCE)
set(POLYSOLVE_WITH_SPECTRA    OFF CACHE BOOL "Enable Spectra library"  FORCE)
set(POLYSOLVE_WITH_SPQR       OFF CACHE BOOL "Enable SPQR library"     FORCE)
set(POLYSOLVE_WITH_SUPERLU    OFF CACHE BOOL "Enable SuperLU library"  FORCE)
set(POLYSOLVE_WITH_TESTS      OFF CACHE BOOL "Build unit-tests"        FORCE)
set(POLYSOLVE_WITH_UMFPACK    OFF CACHE BOOL "Enable UmfPack library"  FORCE)

# Pre-register shared dependencies so this project's pins win over polysolve's
# older ones (both recipe sets guard on the existing targets): finite-diff
# 1.0.4 vs 1.0.2 (the tests need fd::finite_jacobian_tensor), nlohmann/json
# 3.12 vs 3.11, and spdlog 1.17.0 vs 1.9.2.
include(finite_diff)
include(json)
include(spdlog)

include(CPM)
CPMAddPackage("gh:polyfem/polysolve#1b01e6cec13c6813ae48846d1f98654bbbf1402b")

# polysolve always lists its cuSolver backend (CuSolverDN.cu) as a source, even
# with POLYSOLVE_WITH_CUSOLVER off. That is harmless until we enable the CUDA
# language for our own kernels: CMake then compiles the file with nvcc, under
# polysolve's directory-wide MSVC flags (/experimental:external, /FS, ...),
# which nvcc reads as extra input files and rejects on Windows. We never use
# that backend, so we keep the file out of the build.
set_source_files_properties(
    "${polysolve_SOURCE_DIR}/src/polysolve/linear/CuSolverDN.cu"
    TARGET_DIRECTORY polysolve_linear
    PROPERTIES HEADER_FILE_ONLY ON)

# polysolve also hands consumers a bare /utf-8 on MSVC, which would reach nvcc
# in our CUDA tests the same way. We limit it to C++ sources; CUDA sources
# already get -Xcompiler=/utf-8 from ipc_toolkit.
get_target_property(_polysolve_options polysolve_linear INTERFACE_COMPILE_OPTIONS)
if(_polysolve_options)
    list(TRANSFORM _polysolve_options
        REPLACE "^/utf-8$" "$<$<COMPILE_LANGUAGE:CXX>:/utf-8>")
    set_target_properties(polysolve_linear PROPERTIES
        INTERFACE_COMPILE_OPTIONS "${_polysolve_options}")
endif()
unset(_polysolve_options)

# Folder name for IDE
set_target_properties(polysolve PROPERTIES FOLDER "ThirdParty")
set_target_properties(polysolve_linear PROPERTIES FOLDER "ThirdParty")
