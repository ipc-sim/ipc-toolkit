#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace JGSL {

//NOTE: squared distance

template <class T, int dim>
void Point_Point_Distance(const Eigen::Matrix<T, dim, 1>& a, 
    const Eigen::Matrix<T, dim, 1>& b, 
    T& dist2)
{
    dist2 = (a - b).squaredNorm();
}

template <class T, int dim>
void Point_Point_Distance_Gradient(const Eigen::Matrix<T, dim, 1>& a, 
    const Eigen::Matrix<T, dim, 1>& b, 
    Eigen::Matrix<T, dim * 2, 1>& grad)
{
    grad.template segment<dim>(0) = 2.0 * (a - b);
    grad.template segment<dim>(dim) = -grad.template segment<dim>(0);
}

template <class T, int dim>
void Point_Point_Distance_Hessian(const Eigen::Matrix<T, dim, 1>& a, 
    const Eigen::Matrix<T, dim, 1>& b, 
    Eigen::Matrix<T, dim * 2, dim * 2>& Hessian)
{
    Hessian.setZero();
    Hessian.diagonal().setConstant(2.0);
    if constexpr (dim == 2) {
        Hessian(0, 2) = Hessian(1, 3) = Hessian(2, 0) = Hessian(3, 1) = -2.0;
    }
    else {
        Hessian(0, 3) = Hessian(1, 4) = Hessian(2, 5) = Hessian(3, 0) = Hessian(4, 1) = Hessian(5, 2) = -2.0;
    }
}

}