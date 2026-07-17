#pragma once

#include <cstdint>

namespace ipc {

/// @brief Tag identifying a barrier function type.
///
/// Used where virtual dispatch through the Barrier class hierarchy is not
/// possible (e.g. selecting the barrier function inside CUDA kernels).
enum class BarrierType : uint8_t {
    /// @brief ClampedLogBarrier from [Li et al. 2020].
    CLAMPED_LOG,
    /// @brief NormalizedClampedLogBarrier from [Li et al. 2023].
    NORMALIZED_CLAMPED_LOG,
    /// @brief ClampedLogSqBarrier from [Huang et al. 2024].
    CLAMPED_LOG_SQ,
    /// @brief CubicBarrier from [Ando 2024].
    CUBIC,
    /// @brief TwoStageBarrier from [Chen et al. 2025].
    TWO_STAGE,
};

} // namespace ipc
