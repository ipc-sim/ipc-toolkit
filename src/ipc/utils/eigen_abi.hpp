#pragma once

#include <Eigen/Core>

#include <cstddef>

namespace ipc {

/// @brief The parts of Eigen's ABI that depend on how a translation unit was
/// compiled.
///
/// Eigen picks its alignment, and so the size and layout of fixed-size types
/// and the way it allocates memory, from the flags each translation unit is
/// compiled with. If the library and the code calling it disagree, an Eigen
/// object that crosses the boundary is read or freed with the wrong layout. We
/// compare these values to catch that disagreement before it turns into a
/// crash.
struct EigenAbi {
    /// @brief Maximum alignment Eigen uses for heap allocations (EIGEN_MAX_ALIGN_BYTES).
    int max_align_bytes;
    /// @brief Maximum alignment Eigen uses for stack objects (EIGEN_MAX_STATIC_ALIGN_BYTES).
    int max_static_align_bytes;
    /// @brief Whether vectorization is disabled (EIGEN_DONT_VECTORIZE).
    bool dont_vectorize;
    /// @brief Alignment of a fixed-size vectorizable type.
    std::size_t matrix4d_alignment;
    /// @brief Size of a fixed-size vectorizable type that Eigen may pad.
    std::size_t vector_max_size;

    bool operator==(const EigenAbi& other) const
    {
        return max_align_bytes == other.max_align_bytes
            && max_static_align_bytes == other.max_static_align_bytes
            && dont_vectorize == other.dont_vectorize
            && matrix4d_alignment == other.matrix4d_alignment
            && vector_max_size == other.vector_max_size;
    }

    bool operator!=(const EigenAbi& other) const { return !(*this == other); }
};

/// @brief Get the Eigen ABI of the translation unit that calls this function.
///
/// We give this internal linkage (`static`) on purpose. An `inline` function
/// with external linkage would be merged by the linker into a single copy, so
/// the library and its consumer would both call the same one and always
/// report the same values. With `static`, every translation unit keeps its own
/// copy and reports the flags it was compiled with.
static inline EigenAbi eigen_abi()
{
    return EigenAbi {
        EIGEN_MAX_ALIGN_BYTES,
        EIGEN_MAX_STATIC_ALIGN_BYTES,
#ifdef EIGEN_DONT_VECTORIZE
        true,
#else
        false,
#endif
        alignof(Eigen::Matrix4d),
        sizeof(Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 18, 1>),
    };
}

/// @brief Get the Eigen ABI the IPC Toolkit library was compiled with.
EigenAbi library_eigen_abi();

/// @brief Check that the calling translation unit uses the same Eigen ABI as the library.
///
/// If this returns false, mixing Eigen objects between your code and the IPC
/// Toolkit is unsafe. Link against the `ipc::toolkit` CMake target, which
/// carries the right compile definitions, or define the same Eigen macros
/// (e.g. `EIGEN_DONT_VECTORIZE`) everywhere you include Eigen.
static inline bool eigen_abi_matches_library()
{
    return eigen_abi() == library_eigen_abi();
}

} // namespace ipc
