#pragma once

#include <ipc/config.hpp>

namespace ipc::cuda {

#ifdef IPC_TOOLKIT_WITH_CUDA

/// @brief Check whether the ipc::cuda implementations can run here.
///
/// A CUDA build still runs on machines without a GPU or an NVIDIA driver; its
/// ipc::cuda calls then throw on first use. Call this first to choose between
/// the GPU and CPU implementations.
///
/// @return True if the CUDA runtime finds at least one usable device.
bool is_available();

#else

/// @brief Check whether the ipc::cuda implementations can run here.
/// @return Always false: this build has no CUDA support.
inline bool is_available() { return false; }

#endif

} // namespace ipc::cuda
