#pragma once

#include <ipc/config.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace py::literals;

#include <Eigen/Core>
#include <Eigen/Sparse>

/// @brief Check that a parameter is strictly positive.
/// @throws py::value_error (Python ValueError) if value is not positive.
/// @note The C++ API only enforces this with an assert(), which is compiled
///       out under NDEBUG. Validating here means Python users get an
///       exception rather than undefined behavior in a release build.
/// @note Written as !(value > 0) so that NaN is rejected too.
inline void assert_positive(const double value, const std::string& name)
{
    if (!(value > 0)) {
        throw py::value_error(
            "Parameter " + name + " has invalid value: expected " + name
            + " > 0 but got " + name + " = " + std::to_string(value));
    }
}

/// @brief Check that a shared_ptr parameter is not None.
/// @throws py::value_error (Python ValueError) if ptr is null.
template <typename T>
inline void
assert_not_none(const std::shared_ptr<T>& ptr, const std::string& name)
{
    if (ptr == nullptr) {
        throw py::value_error(
            "Parameter " + name + " has invalid value: expected " + name
            + " to be a valid object but got None");
    }
}

template <typename Derived>
void assert_2D_or_3D_vector(
    const Eigen::MatrixBase<Derived>& v, const std::string& name)
{
    if ((v.size() != 2 && v.size() != 3) || (v.rows() != 1 && v.cols() != 1)) {
        throw py::value_error(
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
        throw py::value_error(
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
        throw py::value_error(
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