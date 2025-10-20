#pragma once

#include <ipc/config.hpp> // define DERIVATIVES_WITH_AUTODIFF
#include <ipc/utils/eigen_ext.hpp>

#include <TinyAD/ScalarFunction.hh>

namespace ipc {

using ScalarBase = TinyAD::ScalarBase;

template <int dim> using ADGrad = TinyAD::Scalar<dim, double, false>;
template <int dim> using ADHessian = TinyAD::Scalar<dim, double, true>;

template <typename T> struct IsADGrad {
    static constexpr bool value = false; // NOLINT
};
template <int dim> struct IsADGrad<ADGrad<dim>> {
    static constexpr bool value = true; // NOLINT
};

template <typename T> struct IsADHessian {
    static constexpr bool value = false; // NOLINT
};
template <int dim> struct IsADHessian<ADHessian<dim>> {
    static constexpr bool value = true; // NOLINT
};

} // namespace ipc
