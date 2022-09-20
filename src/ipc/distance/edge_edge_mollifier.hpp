#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the squared norm of the edge-edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The squared norm of the edge-edge cross product.
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
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    return cross(ea1 - ea0, eb1 - eb0).squaredNorm();
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

/// @brief Compute the gradient of the squared norm of the edge cross product.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param[out] grad The gradient of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedGrad>
void edge_edge_cross_squarednorm_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    grad.resize(ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::edge_edge_cross_squarednorm_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
}

/// @brief Compute the hessian of the squared norm of the edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The squared norm of the edge-edge cross product.
/// @param[out] hess The hessian of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename DerivedHess>
void edge_edge_cross_squarednorm_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    hess.resize(
        ea0.size() + ea1.size() + eb0.size() + eb1.size(),
        ea0.size() + ea1.size() + eb0.size() + eb1.size());
    autogen::edge_edge_cross_squarednorm_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
}

/// @brief Mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The mollifier coefficient to premultiply the edge-edge distance.
template <typename T> T edge_edge_mollifier(const T& x, double eps_x)
{
    T x_div_eps_x = x / eps_x;
    return (-x_div_eps_x + 2.0) * x_div_eps_x;
}

/// @brief The gradient of the mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The gradient of the mollifier function for edge-edge distance wrt x.
template <typename T> T edge_edge_mollifier_gradient(const T& x, double eps_x)
{
    T one_div_eps_x = 1.0 / eps_x;
    return 2.0 * one_div_eps_x * (-one_div_eps_x * x + 1.0);
}

/// @brief The hessian of the mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The hessian of the mollifier function for edge-edge distance wrt x.
template <typename T> T edge_edge_mollifier_hessian(const T& x, double eps_x)
{
    return -2.0 / (eps_x * eps_x);
}

/// @brief Compute a mollifier for the edge-edge distance.
///
/// This helps smooth the non-smoothness at close to parallel edges.
///
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param eps_x Mollifier activation threshold.
/// @return The mollifier coefficient to premultiply the edge-edge distance.
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
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    auto ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        return edge_edge_mollifier(ee_cross_norm_sqr, eps_x);
    } else {
        return decltype(ee_cross_norm_sqr)(1.0);
    }
}

/// @brief Compute the gradient of the mollifier for the edge-edge distance.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param eps_x[in] Mollifier activation threshold.
/// @param[out] grad The gradient of the mollifier.
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
    Eigen::PlainObjectBase<DerivedGrad>& grad)
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

/// @brief Compute the hessian of the mollifier for the edge-edge distance.
/// @param[in] ea0 The first vertex of the first edge.
/// @param[in] ea1 The second vertex of the first edge.
/// @param[in] eb0 The first vertex of the second edge.
/// @param[in] eb1 The second vertex of the second edge.
/// @param eps_x[in] Mollifier activation threshold.
/// @param[out] hess The hessian of the mollifier.
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
    Eigen::PlainObjectBase<DerivedHess>& hess)
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
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @return Threshold for edge-edge mollification.
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
    assert(ea0_rest.size() == 3);
    assert(ea1_rest.size() == 3);
    assert(eb0_rest.size() == 3);
    assert(eb1_rest.size() == 3);

    return 1.0e-3 * (ea0_rest - ea1_rest).squaredNorm()
        * (eb0_rest - eb1_rest).squaredNorm();
}

} // namespace ipc
