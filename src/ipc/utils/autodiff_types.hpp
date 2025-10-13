#pragma once

#include <TinyAD/ScalarFunction.hh>
#include "eigen_ext.hpp"

namespace ipc {
template <int dim>
using ADGrad = TinyAD::Scalar<dim, double, false>;
template <int dim>
using ADHessian = TinyAD::Scalar<dim, double, true>;

using ScalarBase = TinyAD::ScalarBase;

template <
    typename T,
    int nrows,
    int ncols,
    int options,
    int maxdim1,
    int maxdim2>
Eigen::Matrix<double, nrows, ncols, options, maxdim1, maxdim2>
autodiff_to_double(
    const Eigen::Matrix<T, nrows, ncols, options, maxdim1, maxdim2>& A)
{
    Eigen::Matrix<double, nrows, ncols, options, maxdim1, maxdim2> out;
    out.resize(A.rows(), A.cols());
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            out(i, j) = A(i, j).value;
        }
    }

    return out;
}

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

template <class T> class AutoDiffAllocator {
public:
    T operator()(const int i, const double& v) const { return T(i, v); }
};

template <> class AutoDiffAllocator<double> {
public:
    double operator()(const int i, const double& v) const { return v; }
};

} // namespace ipc
