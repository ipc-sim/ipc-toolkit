#pragma once

#include "autodiff.h"
#include "eigen_ext.hpp"

namespace ipc {
template <int dim, int max_dim = dim>
using ADGrad = DScalar1<double, Eigen::Matrix<double, dim, 1, 0, max_dim, 1>>;
template <int dim, int max_dim = dim>
using ADHessian = DScalar2<
    double,
    Eigen::Matrix<double, dim, 1, 0, max_dim, 1>,
    Eigen::Matrix<double, dim, dim, 0, max_dim, max_dim>>;

template <
    typename T,
    int nrows,
    int ncols,
    int options,
    int maxdim1,
    int maxdim2>
Eigen::Matrix<double, nrows, ncols, options, maxdim1, maxdim2> AutoDiffToDouble(
    const Eigen::Matrix<T, nrows, ncols, options, maxdim1, maxdim2>& A)
{
    Eigen::Matrix<double, nrows, ncols, options, maxdim1, maxdim2> out;
    out.resize(A.rows(), A.cols());
    for (int i = 0; i < nrows; i++)
        for (int j = 0; j < ncols; j++)
            out(i, j) = A(i, j).value;

    return out;
}

template <typename T> struct isADGrad {
    static const bool value = false;
};
template <int dim, int max_dim> struct isADGrad<ADGrad<dim, max_dim>> {
    static const bool value = true;
};

template <typename T> struct isADHessian {
    static const bool value = false;
};
template <int dim, int max_dim> struct isADHessian<ADHessian<dim, max_dim>> {
    static const bool value = true;
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
