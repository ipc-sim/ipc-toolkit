#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
    template <class scalar>
    scalar smooth_heaviside(const scalar &x);

    template <class scalar>
    scalar heaviside(const scalar &x);

    template <class scalar>
    scalar smooth_max(const scalar &x, const scalar &y);

    template <class scalar>
    scalar smooth_min(const scalar &x, const scalar &y);
}