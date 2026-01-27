#pragma once

#include <algorithm> // for std::clamp
#include <cstdint>   // for uint64_t

namespace ipc {

/// @brief Expands a 32-bit integer into 64 bits by inserting 1 zero after each bit.
/// @param v The 32-bit integer to expand.
/// @return The expanded 64-bit integer.
inline uint64_t expand_bits_1(uint64_t v)
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
inline uint64_t expand_bits_2(uint64_t v)
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
inline uint64_t morton_2D(double x, double y)
{
    constexpr double scale = 1ULL << 32;
    x = std::clamp(x * scale, 0.0, scale - 1);
    y = std::clamp(y * scale, 0.0, scale - 1);
    uint64_t xx = expand_bits_1(uint64_t(x));
    uint64_t yy = expand_bits_1(uint64_t(y));
    return (xx << 1) | yy;
}

/// @brief Calculates a 63-bit Morton code for the given 3D point located within the unit cube [0,1].
/// @param x The x-coordinate of the point.
/// @param y The y-coordinate of the point.
/// @param z The z-coordinate of the point.
/// @return The 63-bit Morton code.
inline uint64_t morton_3D(double x, double y, double z)
{
    constexpr double scale = 1ULL << 21;
    x = std::clamp(x * scale, 0.0, scale - 1);
    y = std::clamp(y * scale, 0.0, scale - 1);
    z = std::clamp(z * scale, 0.0, scale - 1);
    uint64_t xx = expand_bits_2(uint64_t(x));
    uint64_t yy = expand_bits_2(uint64_t(y));
    uint64_t zz = expand_bits_2(uint64_t(z));
    return (xx << 2) | (yy << 1) | zz;
}

} // namespace ipc