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
CPMAddPackage("gh:polyfem/polysolve#231f7e8ed6fc36aad41f21d9b3bf0cdf7816f78c")

# Folder name for IDE
set_target_properties(polysolve PROPERTIES FOLDER "ThirdParty")
set_target_properties(polysolve_linear PROPERTIES FOLDER "ThirdParty")
