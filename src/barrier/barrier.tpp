#pragma once
#include <ipc/barrier/barrier.hpp>

namespace ipc {

template <typename T> T barrier(const T& d, double dhat)
{
    if (d <= 0.0) {
        return T(std::numeric_limits<double>::infinity());
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = -(d-d̂)²ln(d / d̂)
    const double tmp = (d - dhat);
    return -tmp * tmp * log(d / dhat);
}

template <typename T> T physical_barrier(const T& d, double dhat)
{
    // units(d) = m and units(d̂) = m ⟹ units(b(d)) = m
    // units(κ) = Pa ⟹ units(κ b(d)) = Pa m

    if (d <= 0.0) {
        return T(std::numeric_limits<double>::infinity());
    }
    if (d >= dhat) {
        return T(0);
    }

    // b(d) = -d̂(d/d̂-1)²ln(d / d̂)
    const double d_over_dhat = d / dhat;
    const double d_over_dhat_minus_1 = d_over_dhat - 1;
    return -dhat * d_over_dhat_minus_1 * d_over_dhat_minus_1 * log(d_over_dhat);
}

} // namespace ipc
