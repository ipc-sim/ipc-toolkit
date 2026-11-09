# Pin Eigen's ABI on the Eigen target itself, so every library that links it
# gets the same settings without us patching each dependency one by one.
#
# Eigen picks the alignment, and so the size and layout of fixed-size types and
# the way it allocates memory, per translation unit from the flags that unit was
# compiled with. If two libraries disagree, an Eigen object that crosses their
# boundary is read or freed with the wrong layout. That shows up as crashes that
# only happen on some compilers and CPUs, because the allocator and the layout
# only disagree when the build flags make Eigen's alignment differ.
#
# Setting a definition on the library that happens to use it is not enough,
# because the libraries we depend on (finite-diff, MeshFEMSparse, polysolve, ...)
# are compiled without it. Instead, we attach the definitions to the Eigen target
# as usage requirements. Anything that links Eigen3::Eigen inherits them, and so
# a future dependency is covered automatically.
#
# Usage:
#   ipc_toolkit_pin_eigen_abi(<eigen_target> [DONT_VECTORIZE] [PIN_ALIGNMENT])
#
#   DONT_VECTORIZE  Define EIGEN_DONT_VECTORIZE=1. We use this when SIMD is on,
#                   because the SIMD flags (-march=...) make Eigen's alignment
#                   depend on which translation unit it is compiled in.
#   PIN_ALIGNMENT   Define EIGEN_MAX_STATIC_ALIGN_BYTES=0 and
#                   EIGEN_MAX_ALIGN_BYTES=0. We use this when CUDA is on so that
#                   C++ and CUDA translation units agree on object layout. nvcc
#                   never sees the host compiler's -march flags, so the two sides
#                   disagree whenever anything makes them differ. Two ways that
#                   happens: EIGEN_DONT_VECTORIZE yields align=0 under the host
#                   compiler but align=16 under nvcc, and -march=native yields
#                   align=32 (AVX) or 64 (AVX-512) in C++ against nvcc's 16.
#                   Fixed-max types (e.g. VectorMax<double, 18>) then get a
#                   different sizeof/alignof in each, which smashes the stack
#                   when one is returned by value across that boundary.
function(ipc_toolkit_pin_eigen_abi eigen_target)
  cmake_parse_arguments(ARG "DONT_VECTORIZE;PIN_ALIGNMENT" "" "" ${ARGN})

  set(definitions)
  if(ARG_DONT_VECTORIZE)
    list(APPEND definitions EIGEN_DONT_VECTORIZE=1)
  endif()
  if(ARG_PIN_ALIGNMENT)
    list(APPEND definitions
      EIGEN_MAX_STATIC_ALIGN_BYTES=0
      EIGEN_MAX_ALIGN_BYTES=0)
  endif()
  if(NOT definitions)
    return()
  endif()

  # Compile definitions cannot be added through an ALIAS target, so we resolve
  # it to the real target first. Eigen3::Eigen is a plain target when we create
  # it in cmake/recipes/eigen.cmake or when the parent project provides one.
  get_target_property(real_target ${eigen_target} ALIASED_TARGET)
  if(NOT real_target)
    set(real_target ${eigen_target})
  endif()

  # A parent project may pin the same target itself, so we skip definitions the
  # target already carries instead of adding duplicates.
  get_target_property(existing ${real_target} INTERFACE_COMPILE_DEFINITIONS)
  if(existing)
    list(REMOVE_ITEM definitions ${existing})
  endif()
  if(NOT definitions)
    return()
  endif()

  set_property(TARGET ${real_target} APPEND PROPERTY
    INTERFACE_COMPILE_DEFINITIONS ${definitions})
endfunction()
