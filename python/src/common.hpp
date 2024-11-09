#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

template <typename Derived>
void assert_2D_or_3D_vector(
    const Eigen::MatrixBase<Derived>& v, const std::string& name)
{
    if ((v.size() != 2 && v.size() != 3) || (v.rows() != 1 && v.cols() != 1)) {
        throw pybind11::value_error(
            "Parameter " + name + " has invalid size: expected " + name
            + " to be a 2D or 3D vector but got " + name + ".shape = ["
            + std::to_string(v.rows()) + ", " + std::to_string(v.cols()) + "]");
    }
}

template <typename Derived>
void assert_3D_vector(
    const Eigen::MatrixBase<Derived>& v, const std::string& name)
{
    if (v.size() != 3 || (v.rows() != 1 && v.cols() != 1)) {
        throw pybind11::value_error(
            "Parameter " + name + " has invalid size: expected " + name
            + " to be a 3D vector but got " + name + ".shape = ["
            + std::to_string(v.rows()) + ", " + std::to_string(v.cols()) + "]");
    }
}

template <typename T>
inline void assert_is_sparse_vector(
    const Eigen::SparseMatrix<T>& M, const std::string& name)
{
    if (M.cols() != 1) {
        throw pybind11::value_error(
            "Parameter " + name + " has invalid size: expected " + name
            + " to be a sparse vector but got " + name + ".shape = ["
            + std::to_string(M.rows()) + ", " + std::to_string(M.cols()) + "]");
    }
}

template <typename DerivedV, typename DerivedVCopy>
inline void copy_vector(
    const Eigen::MatrixBase<DerivedV>& v,
    Eigen::MatrixBase<DerivedVCopy>& v_copy)
{
    if (v.cols() == 1) {
        v_copy = v;
    } else {
        v_copy = v.transpose();
    }
}