#pragma once

#include "eigen_ext.hpp"

#include <ipc/config.hpp>

namespace ipc {
enum class HeavisideType { ZERO = 0, ONE = 1, VARIANT = 2 };

struct OrientationTypes {

    static HeavisideType
    compute_type(const double& val, const double& alpha, const double& beta);

    int size() const { return m_size; }
    void set_size(const int size);
    const HeavisideType& tangent_type(const int& i) const
    {
        return tangent_types[i];
    }
    const HeavisideType& normal_type(const int& i) const
    {
        return normal_types[i];
    }
    HeavisideType& tangent_type(const int& i) { return tangent_types[i]; }
    HeavisideType& normal_type(const int& i) { return normal_types[i]; }

    bool are_tangent_types_all_one() const;
    bool exists_normal_type_one() const;

    int m_size = 0;
    std::vector<HeavisideType> tangent_types, normal_types;
};

constexpr double MOLLIFIER_THRESHOLD_EPS = 1e-2;

template <typename scalar> struct Math {
    static double sign(const double& x);
    static scalar abs(const scalar& x);
    static scalar sqr(const scalar& x);
    static scalar cubic(const scalar& x);

    static scalar cubic_spline(const scalar& x);
    static double cubic_spline_grad(const double& x);
    static double cubic_spline_hess(const double& x);

    /// @brief support is [-1, 1]
    /// @tparam scalar
    /// @param x
    /// @return
    static scalar quadratic_spline(const scalar& x);

    static scalar smooth_heaviside(
        const scalar& x, const double alpha, const double beta = 0);
    static double smooth_heaviside_grad(
        const double& x, const double alpha, const double beta = 0);
    static double smooth_heaviside_hess(
        const double& x, const double alpha, const double beta = 0);

    static scalar mollifier(const scalar& x);
    static double mollifier_grad(const double& x);
    static double mollifier_hess(const double& x);

    // support is [0, 1]
    static scalar inv_barrier(const scalar& x, const int& r);
    static double inv_barrier_grad(const double& x, const int& r);
    static double inv_barrier_hess(const double& x, const int& r);

    static scalar L_ns(const scalar& x);

    static scalar cross2(
        const Eigen::Ref<const Vector2<scalar>>& a,
        const Eigen::Ref<const Vector2<scalar>>& b);
};

// gradient is symmetric
std::tuple<Eigen::Vector3d, Eigen::Matrix3d>
normalize_vector_grad(const Eigen::Ref<const Eigen::Vector3d>& t);
// hessian is symmetric wrt. the three dimensions
std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix3d,
    std::array<Eigen::Matrix<double, 3, 3>, 3>>
normalize_vector_hess(const Eigen::Ref<const Eigen::Vector3d>& t);

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
    assert(cols > 0);
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
                points(i, d) = T(id + offset, positions(id));
            }
        }
    }

    return points;
}

Eigen::Matrix<double, 3, 6> cross_product_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

std::array<Matrix6d, 3> cross_product_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

// assume unit vector d
double opposite_direction_penalty(
    const Eigen::Ref<const Eigen::Vector3d>& t,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta);

GradType<6> opposite_direction_penalty_grad(
    const Eigen::Ref<const Eigen::Vector3d>& t,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta);

HessianType<6> opposite_direction_penalty_hess(
    const Eigen::Ref<const Eigen::Vector3d>& t,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta);

// assume unit vector d
double negative_orientation_penalty(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta);

GradType<9> negative_orientation_penalty_grad(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta);

HessianType<9> negative_orientation_penalty_hess(
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2,
    const Eigen::Ref<const Eigen::Vector3d>& d,
    const double& alpha,
    const double& beta);

} // namespace ipc

#include "math.tpp"
