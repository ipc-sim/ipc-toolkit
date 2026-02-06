#pragma once

#include <ipc/config.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

enum class HeavisideType : uint8_t { ZERO = 0, ONE = 1, VARIANT = 2 };

struct OrientationTypes {

    static HeavisideType
    compute_type(const double val, const double alpha, const double beta);

    int size() const { return m_size; }
    void set_size(const int size);
    const HeavisideType& tangent_type(const int i) const
    {
        return tangent_types[i];
    }
    const HeavisideType& normal_type(const int i) const
    {
        return normal_types[i];
    }
    HeavisideType& tangent_type(const int i) { return tangent_types[i]; }
    HeavisideType& normal_type(const int i) { return normal_types[i]; }

    int m_size = 0;
    std::vector<HeavisideType> tangent_types, normal_types;
};

constexpr double MOLLIFIER_THRESHOLD_EPS = 1e-2;

template <typename T> struct Math {
    Math() = delete;
    Math(const Math&) = delete;
    Math& operator=(const Math&) = delete;

    // NOTE: Define these in the class definition to allow inlining
    static double sign(const double x) { return x >= 0 ? 1.0 : -1.0; }
    static T abs(const T& x) { return x >= 0 ? x : -x; }
    static T sqr(const T& x) { return x * x; }
    static T cubic(const T& x) { return x * x * x; }

    static T cubic_spline(const T& x);
    static double cubic_spline_grad(const double x);
    static double cubic_spline_hess(const double x);

    /// @brief support is [-1, 1]
    /// @tparam T
    /// @param x
    /// @return
    static T quadratic_spline(const T& x);

    static T
    smooth_heaviside(const T& x, const double alpha, const double beta = 0);
    static double smooth_heaviside_grad(
        const double x, const double alpha, const double beta = 0);
    static double smooth_heaviside_hess(
        const double x, const double alpha, const double beta = 0);

    static T mollifier(const T& x);
    static double mollifier_grad(const double x);
    static double mollifier_hess(const double x);

    // support is [0, 1]
    static T inv_barrier(const T& x, const int r);
    static double inv_barrier_grad(const double x, const int r);
    static double inv_barrier_hess(const double x, const int r);

    static T l_ns(const T& x);

    // NOTE: Define this in the class definition to allow inlining
    static T cross2(
        Eigen::ConstRef<Eigen::Vector2<T>> a,
        Eigen::ConstRef<Eigen::Vector2<T>> b)
    {
        return a[0] * b[1] - a[1] * b[0];
    }
};

template <class T, int rows, int cols, int max_rows = rows>
inline Eigen::Matrix<
    T,
    rows,
    cols,
    (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor),
    max_rows,
    cols>
slice_positions(const Eigen::VectorXd& positions, const int offset = 0)
{
    static_assert(cols > 0);
    const int nrows = rows > 0 ? rows : positions.size() / cols;
    Eigen::Matrix<
        T, rows, cols, (max_rows > 1 ? Eigen::ColMajor : Eigen::RowMajor),
        max_rows, cols>
        points;
    points.setZero(nrows, cols);

    for (int i = 0, id = 0; i < nrows; i++) {
        for (int d = 0; d < cols; d++, id++) {
            if constexpr (std::is_same<T, double>::value) {
                points(i, d) = positions(id);
            } else {
                points(i, d) = T(positions(id), id + offset);
            }
        }
    }

    return points;
}

Eigen::Matrix<double, 3, 6> cross_product_gradient(
    Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2);

std::array<Matrix6d, 3> cross_product_hessian(
    Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2);

// assume unit vector d
double opposite_direction_penalty(
    Eigen::ConstRef<Eigen::Vector3d> t,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta);

GradientType<6> opposite_direction_penalty_grad(
    Eigen::ConstRef<Eigen::Vector3d> t,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta);

HessianType<6> opposite_direction_penalty_hess(
    Eigen::ConstRef<Eigen::Vector3d> t,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta);

// assume unit vector d
double negative_orientation_penalty(
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta);

GradientType<9> negative_orientation_penalty_grad(
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta);

HessianType<9> negative_orientation_penalty_hess(
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    Eigen::ConstRef<Eigen::Vector3d> d,
    const double alpha,
    const double beta);

} // namespace ipc

#include "math.tpp"
