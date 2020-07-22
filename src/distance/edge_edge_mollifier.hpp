#pragma once

#include <distance/distance_type.hpp>
#include <distance/line_line.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_point.hpp>

namespace ipc {

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
auto edge_edge_cross_squarednorm(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    return (ea1 - ea0).cross(eb1 - eb0).squaredNorm();
}

// Symbolically generated derivatives;
namespace autogen {
    void edge_edge_cross_squarednorm_gradient(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double g[12]);

    void edge_edge_cross_squarednorm_hessian(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double H[144]);
} // namespace autogen

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
void edge_edge_cross_squarednorm_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::VectorXd& grad)
{
    grad.resize(ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::edge_edge_cross_squarednorm_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
void edge_edge_cross_squarednorm_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::MatrixXd& hess)
{
    hess.resize(
        ea0.size() + ea1.size() + eb0.size() + eb1.size(),
        ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::edge_edge_cross_squarednorm_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
}

template <class T> T edge_edge_mollifier(const T& x, double eps_x)
{
    T x_div_eps_x = x / eps_x;
    return (-x_div_eps_x + 2.0) * x_div_eps_x;
}

template <class T> T edge_edge_mollifier_gradient(const T& x, double eps_x)
{
    T one_div_eps_x = 1.0 / eps_x;
    return 2.0 * one_div_eps_x * (-one_div_eps_x * x + 1.0);
}

template <class T> T edge_edge_mollifier_hessian(const T& x, double eps_x)
{
    return -2.0 / (eps_x * eps_x);
}

/// @brief Compute a mollifier for the edge-edge distance.
///
/// This helps smooth the non-smoothness at close to parallel edges.
///
/// @param ea0,ea1 The points of the first edge.
/// @param eb0,eb1 The points of the second edge.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_mollifier(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const double eps_x)
{
    auto ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        return edge_edge_mollifier(ee_cross_norm_sqr, eps_x);
    } else {
        return decltype(ee_cross_norm_sqr)(1.0);
    }
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedGrad>
void edge_edge_mollifier_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const double eps_x,
    Eigen::MatrixBase<DerivedGrad>& grad)
{
    int dim = ea0.size();
    assert(ea1.size() == dim);
    assert(eb0.size() == dim);
    assert(eb1.size() == dim);

    grad.resize(4 * dim);

    auto ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1, grad);
        grad *= edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x);
    } else {
        grad.setZero();
    }
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedHess>
void edge_edge_mollifier_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const double eps_x,
    Eigen::MatrixBase<DerivedHess>& hess)
{
    int dim = ea0.size();
    assert(ea1.size() == dim);
    assert(eb0.size() == dim);
    assert(eb1.size() == dim);

    hess.resize(4 * dim, 4 * dim);

    auto ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1, hess);
        hess *= edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x);

        Eigen::Matrix<typename DerivedHess::Scalar, Eigen::Dynamic, 1> grad(
            4 * dim);
        edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1, grad);

        hess += (edge_edge_mollifier_hessian(ee_cross_norm_sqr, eps_x) * grad)
            * grad.transpose();
    } else {
        hess.setZero();
    }
}

/// @brief Compute the threshold of the mollifier edge-edge distance.
///
/// This values is computed based on the edges at rest length.
///
/// @param ea0,ea1 The points of the first edge at rest.
/// @param eb0,eb1 The points of the second edge at rest.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
double edge_edge_mollifier_threshold(
    const Eigen::MatrixBase<DerivedEA0>& ea0_rest,
    const Eigen::MatrixBase<DerivedEA1>& ea1_rest,
    const Eigen::MatrixBase<DerivedEB0>& eb0_rest,
    const Eigen::MatrixBase<DerivedEB1>& eb1_rest)
{
    return 1.0e-3 * (ea0_rest - ea1_rest).squaredNorm()
        * (eb0_rest - eb1_rest).squaredNorm();
}

} // namespace ipc
