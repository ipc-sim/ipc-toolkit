#pragma once
#include <ipc/barrier/barrier.hpp>

namespace ipc {

template <typename T> T barrier(const T& d, double dhat)
{
    if (d <= T(0)) {
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
    if (d <= T(0)) {
        return T(std::numeric_limits<double>::infinity());
    }
    if (d >= dhat) {
        return T(0);
    }

    // b(d) = -d̂(d/d̂-1)²ln(d / d̂)
    const double tmp = (d / dhat - 1);
    return dhat * tmp * tmp * log(d / dhat);
}

} // namespace ipc
