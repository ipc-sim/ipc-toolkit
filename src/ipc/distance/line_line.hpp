#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Geometry>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    // clang-format off
    template <typename T>
    void line_line_distance_gradient(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T g[12]);
    template <typename T>
    void line_line_distance_hessian(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T H[144]);
    // clang-format on
} // namespace autogen

namespace detail {
    /// @brief Compute the distance between a two infinite lines in 3D.
    /// @note The distance is actually squared distance.
    /// @warning If the lines are parallel this function returns a distance of zero.
    /// @param ea0 The first vertex of the edge defining the first line.
    /// @param ea1 The second vertex of the edge defining the first line.
    /// @param eb0 The first vertex of the edge defining the second line.
    /// @param eb1 The second vertex of the edge defining the second line.
    /// @return The distance between the two lines.
    template <typename T>
    inline T line_line_distance(
        Eigen::ConstRef<Eigen::Vector3<T>> ea0,
        Eigen::ConstRef<Eigen::Vector3<T>> ea1,
        Eigen::ConstRef<Eigen::Vector3<T>> eb0,
        Eigen::ConstRef<Eigen::Vector3<T>> eb1)
    {
        const Eigen::Vector3<T> normal = (ea1 - ea0).cross(eb1 - eb0);
        const T line_to_line = (eb0 - ea0).dot(normal);
        return line_to_line * line_to_line / normal.squaredNorm();
    }
} // namespace detail

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_distance(
    const DerivedEA0& ea0,
    const DerivedEA1& ea1,
    const DerivedEB0& eb0,
    const DerivedEB1& eb1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedEA0, DerivedEA1, DerivedEB0, DerivedEB1);
    using T = typename DerivedEA0::Scalar;
    // NOTE: explicit <T>: the detail overload cannot deduce T from an
    // arbitrary Eigen expression.
    return detail::line_line_distance<T>(ea0, ea1, eb0, eb1);
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_distance_gradient(
    const DerivedEA0& ea0,
    const DerivedEA1& ea1,
    const DerivedEB0& eb0,
    const DerivedEB1& eb1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedEA0, DerivedEA1, DerivedEB0, DerivedEB1);
    Eigen::Vector<typename DerivedEA0::Scalar, 12> grad;
    autogen::line_line_distance_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
    return grad;
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto line_line_distance_hessian(
    const DerivedEA0& ea0,
    const DerivedEA1& ea1,
    const DerivedEB0& eb0,
    const DerivedEB1& eb1)
{
    IPC_ASSERT_EIGEN_ARGS(DerivedEA0, DerivedEA1, DerivedEB0, DerivedEB1);
    Eigen::Matrix<typename DerivedEA0::Scalar, 12, 12> hess;
    autogen::line_line_distance_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
    return hess;
}

} // namespace ipc
