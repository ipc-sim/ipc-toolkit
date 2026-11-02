# finite-diff (https://github.com/zfergus/finite-diff)
# License: MIT
if(TARGET finitediff::finitediff)
    return()
endif()

message(STATUS "Third-party: creating target 'finitediff::finitediff'")

include(CPM)
CPMAddPackage("gh:zfergus/finite-diff@1.0.4")

# ipc_toolkit compiles with EIGEN_DONT_VECTORIZE=1 when SIMD is enabled, which
# changes how Eigen aligns and allocates matrices. finite-diff takes matrices
# that our lambdas allocate and then reshapes and frees them inside its own
# library code. If the two sides disagree about alignment, the library either
# trips Eigen's "data is not aligned" assert (Debug) or frees a pointer that
# was never allocated with Eigen's aligned allocator (Release). We make the
# define PUBLIC so both sides always agree, the same way we do for MeshFEMSparse.
if(IPC_TOOLKIT_WITH_SIMD)
    target_compile_definitions(finitediff_finitediff PUBLIC EIGEN_DONT_VECTORIZE=1)
endif()

# Folder name for IDE
set_target_properties(finitediff_finitediff PROPERTIES FOLDER "ThirdParty")