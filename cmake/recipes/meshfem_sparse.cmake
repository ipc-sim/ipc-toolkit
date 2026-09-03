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
if(TARGET MeshFEM::Sparse)
    return()
endif()

if(TARGET MeshFEMSparse)
    # We rely on a MeshFEM::Sparse alias of MeshFEMSparse, so define that alias and return.
    add_library(MeshFEM::Sparse ALIAS MeshFEMSparse)
    return()
endif()

message(STATUS "Third-party: creating target 'MeshFEM::Sparse'")

# Required dependencies of MeshFEMSparse:
include(eigen)
include(onetbb)
find_package(Threads REQUIRED)

include(CPM)
CPMAddPackage(
    URI "gh:MeshFEM/MeshFEMCore#de0d68a16528e59fc491cf27b9f61d37ac464365"
    DOWNLOAD_ONLY YES
)
CPMAddPackage(
    URI "gh:MeshFEM/MeshFEMSparse#0a113141e8288fb1369ecfd465bc12dc09cae77a"
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

# ipc_toolkit is compiled with EIGEN_DONT_VECTORIZE=1 when SIMD is enabled;
# compiling the same Eigen templates with different vectorization settings is
# an ODR violation with real alignment/layout consequences. PUBLIC so the
# setting travels with the target: anything including MeshFEMSparse's headers
# must agree with how its own translation units were compiled.
if(IPC_TOOLKIT_WITH_SIMD)
    target_compile_definitions(MeshFEMSparse PUBLIC EIGEN_DONT_VECTORIZE=1)
endif()

# Folder name for IDE
set_target_properties(MeshFEMSparse PROPERTIES FOLDER "ThirdParty")
