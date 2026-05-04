#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    // clang-format off
    template <typename T>
    void edge_edge_cross_squarednorm_gradient(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T g[12]);
    template <typename T>
    void edge_edge_cross_squarednorm_hessian(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T H[144]);
    template <typename T>
    void edge_edge_mollifier_threshold_gradient(
        T ea0x, T ea0y, T ea0z, T ea1x, T ea1y, T ea1z, T eb0x, T eb0y, T eb0z, T eb1x, T eb1y, T eb1z, T grad[12], T scale = T(1e-3));
    // clang-format on
} // namespace autogen

/// @brief Compute the squared norm of the edge-edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The squared norm of the edge-edge cross product.
template <typename T>
inline T edge_edge_cross_squarednorm(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    return (ea1 - ea0).cross(eb1 - eb0).squaredNorm();
}

/// @brief Compute the gradient of the squared norm of the edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The gradient of the squared norm of the edge cross product wrt ea0,
///     ea1, eb0, and eb1.
template <typename T>
inline Eigen::Vector<T, 12> edge_edge_cross_squarednorm_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    Eigen::Vector<T, 12> grad;
    autogen::edge_edge_cross_squarednorm_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
    return grad;
}

/// @brief Compute the hessian of the squared norm of the edge cross product.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The hessian of the squared norm of the edge cross product wrt ea0,
///     ea1, eb0, and eb1.
template <typename T>
inline Eigen::Matrix<T, 12, 12> edge_edge_cross_squarednorm_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    Eigen::Matrix<T, 12, 12> hess;
    autogen::edge_edge_cross_squarednorm_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
    return hess;
}

/// @brief Mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The mollifier coefficient to premultiply the edge-edge distance.
template <typename T> inline T edge_edge_mollifier(const T x, const T eps_x)
{
    if (x < eps_x) {
        const T x_div_eps_x = x / eps_x;
        return (-x_div_eps_x + T(2)) * x_div_eps_x;
    } else {
        return T(1);
    }
}

/// @brief The gradient of the mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The gradient of the mollifier function for edge-edge distance wrt x.
template <typename T>
inline T edge_edge_mollifier_gradient(const T x, const T eps_x)
{
    if (x < eps_x) {
        const T one_div_eps_x = T(1) / eps_x;
        return T(2) * one_div_eps_x * fma(-one_div_eps_x, x, T(1));
    } else {
        return T(0);
    }
}

/// @brief The derivative of the mollifier function for edge-edge distance wrt
///     eps_x.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The derivative of the mollifier function for edge-edge distance wrt
///     eps_x.
template <typename T>
inline T edge_edge_mollifier_derivative_wrt_eps_x(const T x, const T eps_x)
{
    return x < eps_x ? (T(2) * x * (-eps_x + x) / (eps_x * eps_x * eps_x))
                     : T(0);
}

/// @brief The hessian of the mollifier function for edge-edge distance.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The hessian of the mollifier function for edge-edge distance wrt x.
template <typename T>
inline T edge_edge_mollifier_hessian(const T x, const T eps_x)
{
    if (x < eps_x) {
        return T(-2) / (eps_x * eps_x);
    } else {
        return T(0);
    }
}

/// @brief The derivative of the gradient of the mollifier function for
///     edge-edge distance wrt eps_x.
/// @param x Squared norm of the edge-edge cross product.
/// @param eps_x Mollifier activation threshold.
/// @return The derivative of the gradient of the mollifier function for
///     edge-edge distance wrt eps_x.
template <typename T>
inline T
edge_edge_mollifier_gradient_derivative_wrt_eps_x(const T x, const T eps_x)
{
    return x < eps_x ? (T(2) * (-eps_x + T(2) * x) / (eps_x * eps_x * eps_x))
                     : T(0);
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
template <typename T>
inline T edge_edge_mollifier(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    const T eps_x)
{
    const T ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        return edge_edge_mollifier(ee_cross_norm_sqr, eps_x);
    } else {
        return T(1);
    }
}

/// @brief Compute the gradient of the mollifier for the edge-edge distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param eps_x Mollifier activation threshold.
/// @return The gradient of the mollifier.
template <typename T>
inline Eigen::Vector<T, 12> edge_edge_mollifier_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    const T eps_x)
{
    const T ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        return edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x)
            * edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1);
    } else {
        return Eigen::Vector<T, 12>::Zero();
    }
}

/// @brief Compute the hessian of the mollifier for the edge-edge distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param eps_x Mollifier activation threshold.
/// @return The hessian of the mollifier.
template <typename T>
inline Eigen::Matrix<T, 12, 12> edge_edge_mollifier_hessian(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    const T eps_x)
{
    const T ee_cross_norm_sqr = edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        const Eigen::Vector<T, 12> grad =
            edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1);

        return (edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x)
                * edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1))
            + ((edge_edge_mollifier_hessian(ee_cross_norm_sqr, eps_x) * grad)
               * grad.transpose());
    } else {
        return Eigen::Matrix<T, 12, 12>::Zero();
    }
}

/// @brief Compute the gradient of the mollifier for the edge-edge distance wrt
///     rest positions.
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The derivative of the mollifier wrt rest positions.
template <typename T>
Eigen::Vector<T, 12> edge_edge_mollifier_gradient_wrt_x(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

/// @brief Compute the jacobian of the edge-edge distance mollifier's gradient
///     wrt rest positions.
/// @note This is not the hessian of the mollifier wrt rest positions, but the
///     jacobian wrt rest positions of the mollifier's gradient wrt positions.
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @return The jacobian of the mollifier's gradient wrt rest positions.
template <typename T>
Eigen::Matrix<T, 12, 12> edge_edge_mollifier_gradient_jacobian_wrt_x(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

/// @brief Compute the threshold of the mollifier edge-edge distance.
///
/// This values is computed based on the edges at rest length.
///
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @return Threshold for edge-edge mollification.
template <typename T>
T edge_edge_mollifier_threshold(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1_rest)
{
    return T(1e-3) * (ea0_rest - ea1_rest).squaredNorm()
        * (eb0_rest - eb1_rest).squaredNorm();
}

/// @brief Compute the gradient of the threshold of the mollifier edge-edge
///     distance.
///
/// This values is computed based on the edges at rest length.
///
/// @param ea0_rest The rest position of the first vertex of the first edge.
/// @param ea1_rest The rest position of the second vertex of the first edge.
/// @param eb0_rest The rest position of the first vertex of the second edge.
/// @param eb1_rest The rest position of the second vertex of the second edge.
/// @return Gradient of the threshold for edge-edge mollification.
template <typename T>
Eigen::Vector<T, 12> edge_edge_mollifier_threshold_gradient(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0_rest,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1_rest)
{
    Eigen::Vector<T, 12> grad;
    autogen::edge_edge_mollifier_threshold_gradient(
        ea0_rest[0], ea0_rest[1], ea0_rest[2], ea1_rest[0], ea1_rest[1],
        ea1_rest[2], eb0_rest[0], eb0_rest[1], eb0_rest[2], eb1_rest[0],
        eb1_rest[1], eb1_rest[2], grad.data(), /*scale=*/T(1e-3));
    return grad;
}

// --- EigenExpression wrappers ---

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_cross_squarednorm(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1) -> typename DerivedEA0::Scalar
{
    using T = typename DerivedEA0::Scalar;
    return edge_edge_cross_squarednorm(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_cross_squarednorm_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Vector<typename DerivedEA0::Scalar, 12>
{
    using T = typename DerivedEA0::Scalar;
    return edge_edge_cross_squarednorm_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_cross_squarednorm_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Matrix<typename DerivedEA0::Scalar, 12, 12>
{
    using T = typename DerivedEA0::Scalar;
    return edge_edge_cross_squarednorm_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

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
    const typename DerivedEA0::Scalar eps_x) -> typename DerivedEA0::Scalar
{
    using T = typename DerivedEA0::Scalar;
    return edge_edge_mollifier(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1), eps_x);
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_mollifier_gradient(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const typename DerivedEA0::Scalar eps_x)
    -> Eigen::Vector<typename DerivedEA0::Scalar, 12>
{
    using T = typename DerivedEA0::Scalar;
    return edge_edge_mollifier_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1), eps_x);
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_mollifier_hessian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const typename DerivedEA0::Scalar eps_x)
    -> Eigen::Matrix<typename DerivedEA0::Scalar, 12, 12>
{
    using T = typename DerivedEA0::Scalar;
    return edge_edge_mollifier_hessian(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1), eps_x);
}

template <
    typename DerivedEA0Rest,
    typename DerivedEA1Rest,
    typename DerivedEB0Rest,
    typename DerivedEB1Rest,
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_mollifier_gradient_wrt_x(
    const Eigen::MatrixBase<DerivedEA0Rest>& ea0_rest,
    const Eigen::MatrixBase<DerivedEA1Rest>& ea1_rest,
    const Eigen::MatrixBase<DerivedEB0Rest>& eb0_rest,
    const Eigen::MatrixBase<DerivedEB1Rest>& eb1_rest,
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Vector<typename DerivedEA0Rest::Scalar, 12>
{
    using T = typename DerivedEA0Rest::Scalar;
    return edge_edge_mollifier_gradient_wrt_x(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0Rest,
    typename DerivedEA1Rest,
    typename DerivedEB0Rest,
    typename DerivedEB1Rest,
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_mollifier_gradient_jacobian_wrt_x(
    const Eigen::MatrixBase<DerivedEA0Rest>& ea0_rest,
    const Eigen::MatrixBase<DerivedEA1Rest>& ea1_rest,
    const Eigen::MatrixBase<DerivedEB0Rest>& eb0_rest,
    const Eigen::MatrixBase<DerivedEB1Rest>& eb1_rest,
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
    -> Eigen::Matrix<typename DerivedEA0Rest::Scalar, 12, 12>
{
    using T = typename DerivedEA0Rest::Scalar;
    return edge_edge_mollifier_gradient_jacobian_wrt_x(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(ea0),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1));
}

template <
    typename DerivedEA0Rest,
    typename DerivedEA1Rest,
    typename DerivedEB0Rest,
    typename DerivedEB1Rest>
inline auto edge_edge_mollifier_threshold(
    const Eigen::MatrixBase<DerivedEA0Rest>& ea0_rest,
    const Eigen::MatrixBase<DerivedEA1Rest>& ea1_rest,
    const Eigen::MatrixBase<DerivedEB0Rest>& eb0_rest,
    const Eigen::MatrixBase<DerivedEB1Rest>& eb1_rest) ->
    typename DerivedEA0Rest::Scalar
{
    using T = typename DerivedEA0Rest::Scalar;
    return edge_edge_mollifier_threshold(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1_rest));
}

template <
    typename DerivedEA0Rest,
    typename DerivedEA1Rest,
    typename DerivedEB0Rest,
    typename DerivedEB1Rest>
inline auto edge_edge_mollifier_threshold_gradient(
    const Eigen::MatrixBase<DerivedEA0Rest>& ea0_rest,
    const Eigen::MatrixBase<DerivedEA1Rest>& ea1_rest,
    const Eigen::MatrixBase<DerivedEB0Rest>& eb0_rest,
    const Eigen::MatrixBase<DerivedEB1Rest>& eb1_rest)
    -> Eigen::Vector<typename DerivedEA0Rest::Scalar, 12>
{
    using T = typename DerivedEA0Rest::Scalar;
    return edge_edge_mollifier_threshold_gradient(
        Eigen::Ref<const Eigen::Vector3<T>>(ea0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(ea1_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb0_rest),
        Eigen::Ref<const Eigen::Vector3<T>>(eb1_rest));
}

} // namespace ipc
