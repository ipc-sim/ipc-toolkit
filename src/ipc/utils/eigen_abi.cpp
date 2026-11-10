#include "eigen_abi.hpp"

#include <ipc/config.hpp>

// CMake pins Eigen's ABI on the Eigen3::Eigen target (see
// cmake/ipc_toolkit/ipc_toolkit_eigen_abi.cmake). If that pin is lost, the
// library and its dependencies compile with different alignment and corrupt
// the heap at run time, so we fail the build instead.
#if defined(IPC_TOOLKIT_WITH_SIMD) && !defined(EIGEN_DONT_VECTORIZE)
#error "IPC_TOOLKIT_WITH_SIMD requires EIGEN_DONT_VECTORIZE on Eigen3::Eigen."
#endif
#if defined(IPC_TOOLKIT_WITH_CUDA)                                             \
    && (EIGEN_MAX_ALIGN_BYTES != 0 || EIGEN_MAX_STATIC_ALIGN_BYTES != 0)
#error "CUDA builds require EIGEN_MAX_*ALIGN_BYTES=0 on Eigen3::Eigen."
#endif

namespace ipc {

// We define this in a .cpp file, instead of in the header, so that it reports
// the flags the library itself was compiled with. It calls this file's own copy
// of the internal-linkage eigen_abi(), not the consumer's.
EigenAbi library_eigen_abi() { return eigen_abi(); }

} // namespace ipc
