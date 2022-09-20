#pragma once
#include "barrier.hpp"

#include <limits>
#include <cmath>

namespace ipc {

template <typename T> T barrier(const T& d, const double dhat)
{
    if (d <= 0.0) {
        return T(std::numeric_limits<double>::infinity());
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = -(d-d̂)²ln(d / d̂)
    const T d_minus_dhat = (d - dhat);
    return -d_minus_dhat * d_minus_dhat * log(d / dhat);
}

} // namespace ipc
