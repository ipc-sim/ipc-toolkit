#pragma once
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/distance/distance_type.hpp>

#include <ipc/smooth_contact/primitives/point3.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>

#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
template <typename PrimitiveA, typename PrimitiveB>
struct PrimitiveDistType { };

template <> struct PrimitiveDistType<Point3, Point3> {
    using type = PointPointDistanceType;
};

template <> struct PrimitiveDistType<Edge3, Point3> {
    using type = PointEdgeDistanceType;
};

template <> struct PrimitiveDistType<Face, Point3> {
    using type = PointTriangleDistanceType;
};

template <> struct PrimitiveDistType<Edge3, Edge3> {
    using type = EdgeEdgeDistanceType;
};

template <typename PrimitiveA, typename PrimitiveB, typename T>
class PrimitiveDistanceTemplate {
    static_assert(
        PrimitiveA::dim == PrimitiveB::dim,
        "Primitives must have the same dimension");
    constexpr static int dim = PrimitiveA::dim;
    constexpr static int n_core_dofs =
        PrimitiveA::n_core_points * PrimitiveA::dim
        + PrimitiveB::n_core_points * PrimitiveB::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr);
};

template <typename PrimitiveA, typename PrimitiveB> class PrimitiveDistance {
    static_assert(
        PrimitiveA::dim == PrimitiveB::dim,
        "Primitives must have the same dimension");
    constexpr static int dim = PrimitiveA::dim;
    constexpr static int n_core_dofs =
        PrimitiveA::n_core_points * PrimitiveA::dim
        + PrimitiveB::n_core_points * PrimitiveB::dim;

public:
    static typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type
    compute_distance_type(const Vector<double, n_core_dofs>& x);

    static double compute_distance(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const long& a,
        const long& b,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);

    // points from primitiveA to primitiveB
    static Vector<double, dim> compute_closest_direction(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const long& a,
        const long& b,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);

    static std::tuple<Vector<double, dim>, Eigen::Matrix<double, dim, n_core_dofs>>
    compute_closest_direction_gradient(
        const Vector<double, n_core_dofs>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
    {
        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADGrad<n_core_dofs>;
        const Vector<T, n_core_dofs> X = slice_positions<T, n_core_dofs, 1>(x);
        const Vector<T, dim> d = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(X, dtype);

        Vector<double, dim> out;
        Eigen::Matrix<double, dim, n_core_dofs> J = Eigen::Matrix<double, dim, n_core_dofs>::Zero();
        for (int i = 0; i < dim; i++)
        {
            out(i) = d(i).getValue();
            J.row(i) = d(i).getGradient();
        }
        return std::make_tuple(out, J);
    }

    static std::tuple<Vector<double, dim>, Eigen::Matrix<double, dim, n_core_dofs>,
        std::array<Eigen::Matrix<double, n_core_dofs, n_core_dofs>, dim>>
    compute_closest_direction_hessian(
        const Vector<double, n_core_dofs>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
    {
        DiffScalarBase::setVariableCount(n_core_dofs);
        using T = ADHessian<n_core_dofs>;
        const Vector<T, n_core_dofs> X = slice_positions<T, n_core_dofs, 1>(x);
        const Vector<T, dim> d = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(X, dtype);

        Vector<double, dim> out;
        Eigen::Matrix<double, dim, n_core_dofs> J = Eigen::Matrix<double, dim, n_core_dofs>::Zero();
        std::array<Eigen::Matrix<double, n_core_dofs, n_core_dofs>, dim> H;
        for (int i = 0; i < dim; i++)
        {
            out(i) = d(i).getValue();
            J.row(i) = d(i).getGradient();
            H[i] = d(i).getHessian();
        }
        return std::make_tuple(out, J, H);
    }

    static std::tuple<double, Vector<double, n_core_dofs + 1>>
    compute_mollifier_gradient(
        const Vector<double, n_core_dofs>& x, const double dist_sqr)
    {
        DiffScalarBase::setVariableCount(n_core_dofs + 1);
        using T = ADGrad<n_core_dofs + 1>;
        const Vector<T, n_core_dofs + 1> X = slice_positions<T, n_core_dofs + 1, 1>((Vector<double, n_core_dofs + 1>() << x, dist_sqr).finished());
        const T out = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(X.head(n_core_dofs), X(n_core_dofs));

        return std::make_tuple(out.getValue(), out.getGradient());
    }

    static std::tuple<double, Vector<double, n_core_dofs + 1>, Eigen::Matrix<double, n_core_dofs + 1, n_core_dofs + 1>>
    compute_mollifier_hessian(
        const Vector<double, n_core_dofs>& x, const double dist_sqr)
    {
        DiffScalarBase::setVariableCount(n_core_dofs + 1);
        using T = ADHessian<n_core_dofs + 1>;
        const Vector<T, n_core_dofs + 1> X = slice_positions<T, n_core_dofs + 1, 1>((Vector<double, n_core_dofs + 1>() << x, dist_sqr).finished());
        const T out = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(X.head(n_core_dofs), X(n_core_dofs));

        return std::make_tuple(out.getValue(), out.getGradient(), out.getHessian());
    }

};

} // namespace ipc

#include "primitive_distance.tpp"
