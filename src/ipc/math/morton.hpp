#pragma once

#include <ipc/config.hpp>     // for IPC_TOOLKIT_HOST_DEVICE
#include <ipc/utils/simd.hpp> // for clamp

#include <Eigen/Core>

#include <cstdint> // for uint64_t

#if !defined(__CUDA_ARCH__) && !defined(__GNUC__) && !defined(__clang__)       \
    && defined(WIN32)
#include <intrin.h> // for __lzcnt / __lzcnt64
#endif

namespace ipc {

/// @brief Expands a 32-bit integer into 64 bits by inserting 1 zero after each bit.
/// @param v The 32-bit integer to expand.
/// @return The expanded 64-bit integer.
IPC_TOOLKIT_HOST_DEVICE inline uint64_t expand_bits_1(uint64_t v)
{
    v = (v | (v << 16)) & 0x0000FFFF0000FFFF;
    v = (v | (v << 8)) & 0x00FF00FF00FF00FF;
    v = (v | (v << 4)) & 0x0F0F0F0F0F0F0F0F;
    v = (v | (v << 2)) & 0x3333333333333333;
    v = (v | (v << 1)) & 0x5555555555555555;
    return v;
}

/// @brief Expands a 21-bit integer into 63 bits by inserting 2 zeros after each bit.
/// @param v The 21-bit integer to expand.
/// @return The expanded 63-bit integer.
IPC_TOOLKIT_HOST_DEVICE inline uint64_t expand_bits_2(uint64_t v)
{
    v = (v | v << 32) & 0x1F00000000FFFF;
    v = (v | v << 16) & 0x1F0000FF0000FF;
    v = (v | v << 8) & 0x100F00F00F00F00F;
    v = (v | v << 4) & 0x10C30C30C30C30C3;
    v = (v | v << 2) & 0x1249249249249249;
    return v;
}

/// @brief Calculates a 64-bit Morton code for the given 2D point located within the unit square [0,1].
/// @param x The x-coordinate of the point.
/// @param y The y-coordinate of the point.
/// @return The 64-bit Morton code.
IPC_TOOLKIT_HOST_DEVICE inline uint64_t morton_2D(double x, double y)
{
    constexpr double scale = 1ULL << 32;
    x = ipc::clamp(x * scale, 0.0, scale - 1);
    y = ipc::clamp(y * scale, 0.0, scale - 1);
    uint64_t xx = expand_bits_1(uint64_t(x));
    uint64_t yy = expand_bits_1(uint64_t(y));
    return (xx << 1) | yy;
}

/// @brief Calculates a 63-bit Morton code for the given 3D point located within the unit cube [0,1].
/// @param x The x-coordinate of the point.
/// @param y The y-coordinate of the point.
/// @param z The z-coordinate of the point.
/// @return The 63-bit Morton code.
IPC_TOOLKIT_HOST_DEVICE inline uint64_t morton_3D(double x, double y, double z)
{
    constexpr double scale = 1ULL << 21;
    x = ipc::clamp(x * scale, 0.0, scale - 1);
    y = ipc::clamp(y * scale, 0.0, scale - 1);
    z = ipc::clamp(z * scale, 0.0, scale - 1);
    uint64_t xx = expand_bits_2(uint64_t(x));
    uint64_t yy = expand_bits_2(uint64_t(y));
    uint64_t zz = expand_bits_2(uint64_t(z));
    return (xx << 2) | (yy << 1) | zz;
}

/// @brief Calculates the Morton code of a box from its center.
///
/// The center is normalized into the unit square/cube by the given domain
/// before being encoded. The domain's width is passed as a reciprocal so this
/// multiplies rather than divides, letting the host and device LBVH builds
/// agree bit-for-bit.
///
/// @param center The center of the box.
/// @param domain_min The minimum corner of the normalization domain.
/// @param domain_width_inv The reciprocal of the normalization domain's width.
/// @param dim The dimension of the simulation (2 or 3).
/// @return The Morton code of the normalized center.
IPC_TOOLKIT_HOST_DEVICE inline uint64_t morton_code(
    const Eigen::Array3d& center,
    const Eigen::Array3d& domain_min,
    const Eigen::Array3d& domain_width_inv,
    const int dim)
{
    const double x = (center.x() - domain_min.x()) * domain_width_inv.x();
    const double y = (center.y() - domain_min.y()) * domain_width_inv.y();
    if (dim == 2) {
        return morton_2D(x, y);
    }
    const double z = (center.z() - domain_min.z()) * domain_width_inv.z();
    return morton_3D(x, y, z);
}

/// @brief Counts the leading zero bits of a 32-bit value.
/// @note Undefined for v == 0, matching the underlying intrinsics.
/// @param v The value to count the leading zeros of.
/// @return The number of leading zero bits.
IPC_TOOLKIT_HOST_DEVICE inline int count_leading_zeros(const uint32_t v)
{
#ifdef __CUDA_ARCH__
    return __clz(static_cast<int>(v));
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_clz(v);
#elif defined(WIN32)
    return static_cast<int>(__lzcnt(v));
#else
#error "count_leading_zeros: no leading-zero-count intrinsic for this compiler"
#endif
}

/// @brief Counts the leading zero bits of a 64-bit value.
/// @note Undefined for v == 0, matching the underlying intrinsics.
/// @param v The value to count the leading zeros of.
/// @return The number of leading zero bits.
IPC_TOOLKIT_HOST_DEVICE inline int count_leading_zeros(const uint64_t v)
{
#ifdef __CUDA_ARCH__
    return __clzll(static_cast<long long>(v));
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(v);
#elif defined(WIN32)
    return static_cast<int>(__lzcnt64(v));
#else
#error "count_leading_zeros: no leading-zero-count intrinsic for this compiler"
#endif
}

/// @brief Computes the length of the common leading-bit prefix of two sorted
/// Morton codes.
///
/// This is the delta of Apetrei [2014]: a larger value means the two positions
/// are separated by a finer split, and so have a nearer common ancestor.
/// Duplicate codes fall back to the leading zeros of the positions' XOR, offset
/// by 32 so that any code-level difference always compares as the shorter
/// prefix.
///
/// @note The two positions must differ (i != j). This holds for every delta the
/// LBVH build evaluates, as it only ever compares adjacent positions.
///
/// @param code_i The Morton code at sorted position i.
/// @param i The first sorted position.
/// @param code_j The Morton code at sorted position j.
/// @param j The second sorted position.
/// @return The length of the common leading-bit prefix.
IPC_TOOLKIT_HOST_DEVICE inline int morton_common_prefix(
    const uint64_t code_i, const int i, const uint64_t code_j, const int j)
{
    if (code_i == code_j) {
        return 32 + count_leading_zeros(static_cast<uint32_t>(i ^ j));
    }
    return count_leading_zeros(code_i ^ code_j);
}

} // namespace ipc