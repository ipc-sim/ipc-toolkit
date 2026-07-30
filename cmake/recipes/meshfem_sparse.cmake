# MeshFEMSparse (https://github.com/MeshFEM/MeshFEMSparse)
# License: MIT
#
# Block-CSC sparse matrix data structures and fast Hessian assembly routines
# from the MeshFEM project (Mohammadian et al., "MeshFEM: A Block-accelerated
# Solver for Nonlinear Finite Elements", SIGGRAPH 2026).
#
# Downloaded with DOWNLOAD_ONLY and compiled into our own minimal static
# library rather than through upstream's CMake. We only need the matrix data
# structures and assembly routines; upstream's target additionally compiles
# sparse direct solver wrappers (which clash with Eigen 5's BLAS
# declarations), declares -fvisibility=hidden PUBLIC, and fetches its own
# Eigen/TBB/MeshFEMCore.
if(TARGET MeshFEMSparse)
    return()
endif()

message(STATUS "Third-party: creating target 'MeshFEMSparse'")

include(eigen)
include(onetbb)
find_package(Threads REQUIRED)

include(CPM)
CPMAddPackage(
    NAME MeshFEMCore
    URL "https://github.com/MeshFEM/MeshFEMCore/archive/24e81c425e85eee3ed79af000e82ef1a75bbe696.zip"
    URL_HASH SHA256=7d505f57af2a4b2fbc01668d57c1296d646c6aa0aa923c9d70e814b55987aff7
    DOWNLOAD_ONLY YES
)
CPMAddPackage(
    NAME MeshFEMSparse
    URL "https://github.com/MeshFEM/MeshFEMSparse/archive/efe1af87cf6f7d04359628552e84c666b8612f4c.zip"
    URL_HASH SHA256=e23cb168e9e709662083116e73876f97c23edd60413108bd70bde5a5a1684dd5
    DOWNLOAD_ONLY YES
)

add_library(MeshFEMSparse STATIC
    # Matrix data structures and assembly (no Solvers/)
    "${MeshFEMSparse_SOURCE_DIR}/src/lib/MeshFEMSparse/BlockCSCHessian.cc"
    "${MeshFEMSparse_SOURCE_DIR}/src/lib/MeshFEMSparse/BorderedSparseHessian.cc"
    # MeshFEMCore support code (parallelism arenas, benchmark stubs, types)
    "${MeshFEMCore_SOURCE_DIR}/src/lib/MeshFEMCore/GlobalBenchmark.cc"
    "${MeshFEMCore_SOURCE_DIR}/src/lib/MeshFEMCore/Parallelism.cc"
    "${MeshFEMCore_SOURCE_DIR}/src/lib/MeshFEMCore/Types.cc"
)

add_library(MeshFEM::Sparse ALIAS MeshFEMSparse)

target_include_directories(MeshFEMSparse SYSTEM PUBLIC
    "${MeshFEMCore_SOURCE_DIR}/src/lib"
    "${MeshFEMSparse_SOURCE_DIR}/src/lib"
    "${CMAKE_CURRENT_BINARY_DIR}/meshfem/exports"
)

# MeshFEMCore's headers include the CMake-generated <MeshFEM_export.h>. We
# build a static library, so the export macros are empty.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/meshfem/exports/MeshFEM_export.h" [[
#pragma once
#define MESHFEM_EXPORT
#define MESHFEM_NO_EXPORT
#define MESHFEM_DEPRECATED
#define MESHFEM_DEPRECATED_EXPORT
#define MESHFEM_DEPRECATED_NO_EXPORT
]])

target_link_libraries(MeshFEMSparse PUBLIC
    Eigen3::Eigen
    TBB::tbb
    Threads::Threads
)

target_compile_features(MeshFEMSparse PUBLIC cxx_std_17)

# Matches upstream MeshFEMCore's PUBLIC definitions (Parallelism.hh compiles
# its TBB code paths only when MESHFEM_WITH_TBB is defined).
target_compile_definitions(MeshFEMSparse PUBLIC
    MESHFEM_WITH_TBB
    NOMINMAX
    _ENABLE_EXTENDED_ALIGNED_STORAGE
    _USE_MATH_DEFINES
)

# MeshFEM targets Eigen 3.4, but Eigen 5 removed
# Eigen::internal::make_coherent, which MeshFEMCore's
# AutomaticDifferentiation.hh references. Force-include a shim reimplementing
# it (see src/ipc/utils/meshfem_eigen_compat.hpp); TUs outside this target
# that include MeshFEM headers must include the shim first themselves.
set(MESHFEM_EIGEN_COMPAT_HEADER
    "${PROJECT_SOURCE_DIR}/src/ipc/utils/meshfem_eigen_compat.hpp")
if(MSVC)
    target_compile_options(MeshFEMSparse PRIVATE "/FI${MESHFEM_EIGEN_COMPAT_HEADER}")
else()
    target_compile_options(MeshFEMSparse PRIVATE "-include" "${MESHFEM_EIGEN_COMPAT_HEADER}")
endif()

# ipc_toolkit is compiled with EIGEN_DONT_VECTORIZE=1 when SIMD is enabled;
# compiling the same Eigen templates with different vectorization settings is
# an ODR violation with real alignment/layout consequences.
if(IPC_TOOLKIT_WITH_SIMD)
    target_compile_definitions(MeshFEMSparse PRIVATE EIGEN_DONT_VECTORIZE=1)
endif()

# Folder name for IDE
set_target_properties(MeshFEMSparse PROPERTIES FOLDER "ThirdParty")
