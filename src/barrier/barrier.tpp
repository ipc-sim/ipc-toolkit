#pragma once
#include <ipc/barrier/barrier.hpp>

namespace ipc {

template <typename T> T barrier(T d, double dhat)
{
    if (d <= T(0)) {
        return T(std::numeric_limits<double>::infinity());
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = -(d-d̂)²ln(d / d̂)
    T dhat_T = T(dhat);
    return -(d - dhat_T) * (d - dhat_T) * log(d / dhat_T);
}

} // namespace ipc
