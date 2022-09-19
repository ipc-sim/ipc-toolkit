#pragma once
#include <ipc/barrier/barrier.hpp>

#include <ipc/config.hpp>

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
    double b = -d_minus_dhat * d_minus_dhat * log(d / dhat);

#ifdef IPC_TOOLKIT_CONVERGENT
    b /= dhat;
#endif

    return b;
}

} // namespace ipc
