#pragma once
#include <ipc/collision_mesh.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/smooth_contact/primitives/edge2.hpp>
#include <ipc/smooth_contact/primitives/edge3.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>
#include <ipc/smooth_contact/primitives/point2.hpp>
#include <ipc/smooth_contact/primitives/point3.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
template <typename PrimitiveA, typename PrimitiveB>
struct PrimitiveDistType { };

template <> struct PrimitiveDistType<Point2, Point2> {
    using type = PointPointDistanceType;
    static constexpr std::string_view NAME = "PointPoint";
};

template <> struct PrimitiveDistType<Point3, Point3> {
    using type = PointPointDistanceType;
    static constexpr std::string_view NAME = "PointPoint";
};

template <> struct PrimitiveDistType<Edge2, Point2> {
    using type = PointEdgeDistanceType;
    static constexpr std::string_view NAME = "EdgePoint";
};

template <> struct PrimitiveDistType<Edge3, Point3> {
    using type = PointEdgeDistanceType;
    static constexpr std::string_view NAME = "EdgePoint";
};

template <> struct PrimitiveDistType<Face, Point3> {
    using type = PointTriangleDistanceType;
    static constexpr std::string_view NAME = "FacePoint";
};

template <> struct PrimitiveDistType<Edge3, Edge3> {
    using type = EdgeEdgeDistanceType;
    static constexpr std::string_view NAME = "EDGE_EDGE";
};

template <typename PrimitiveA, typename PrimitiveB, typename T>
class PrimitiveDistanceTemplate {
    static_assert(
        PrimitiveA::DIM == PrimitiveB::DIM,
        "Primitives must have the same dimension");
    static constexpr int DIM = PrimitiveA::DIM;
    static constexpr int N_CORE_DOFS =
        PrimitiveA::N_CORE_POINTS * PrimitiveA::DIM
        + PrimitiveB::N_CORE_POINTS * PrimitiveB::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
    static Eigen::Matrix<T, DIM, 2> compute_closest_point_pairs(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr);
};

template <typename PrimitiveA, typename PrimitiveB> class PrimitiveDistance {
    static_assert(
        PrimitiveA::DIM == PrimitiveB::DIM,
        "Primitives must have the same dimension");

public:
    static constexpr int DIM = PrimitiveA::DIM;
    static constexpr int N_CORE_DOFS =
        PrimitiveA::N_CORE_POINTS * PrimitiveA::DIM
        + PrimitiveB::N_CORE_POINTS * PrimitiveB::DIM;

    static typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type
    compute_distance_type(const Vector<double, N_CORE_DOFS>& x);

    static double compute_distance(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const index_t a,
        const index_t b,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);

    static GradType<N_CORE_DOFS> compute_distance_gradient(
        const Vector<double, N_CORE_DOFS>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
    {
        using T = TinyADGrad<N_CORE_DOFS>;
        const Vector<T, N_CORE_DOFS> X = T::make_active(x);
        const T d = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_distance(X, dtype);

        return { d.val, d.grad };
    }

    static HessianType<N_CORE_DOFS> compute_distance_hessian(
        const Vector<double, N_CORE_DOFS>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
    {
        using T = TinyADHessian<N_CORE_DOFS>;
        const Vector<T, N_CORE_DOFS> X = T::make_active(x);
        const T d = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_distance(X, dtype);

        return { d.val, d.grad, d.Hess };
    }

    // points from primitiveA to primitiveB
    static Vector<double, DIM> compute_closest_direction(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const index_t a,
        const index_t b,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);

    static std::
        tuple<Vector<double, DIM>, Eigen::Matrix<double, DIM, N_CORE_DOFS>>
        compute_closest_direction_gradient(
            const Vector<double, N_CORE_DOFS>& x,
            typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
    {
        using T = TinyADGrad<N_CORE_DOFS>;
        const Vector<T, N_CORE_DOFS> X = T::make_active(x);
        const Vector<T, DIM> d = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_closest_direction(X, dtype);

        Vector<double, DIM> out;
        Eigen::Matrix<double, DIM, N_CORE_DOFS> J =
            Eigen::Matrix<double, DIM, N_CORE_DOFS>::Zero();
        for (int i = 0; i < DIM; i++) {
            out(i) = d(i).val;
            J.row(i) = d(i).grad;
        }
        return std::make_tuple(out, J);
    }

    static std::tuple<
        Vector<double, DIM>,
        Eigen::Matrix<double, DIM, N_CORE_DOFS>,
        std::array<Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS>, DIM>>
    compute_closest_direction_hessian(
        const Vector<double, N_CORE_DOFS>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
    {
        using T = TinyADHessian<N_CORE_DOFS>;
        const Vector<T, N_CORE_DOFS> X = T::make_active(x);
        const Vector<T, DIM> d = PrimitiveDistanceTemplate<
            PrimitiveA, PrimitiveB, T>::compute_closest_direction(X, dtype);

        Vector<double, DIM> out;
        Eigen::Matrix<double, DIM, N_CORE_DOFS> J =
            Eigen::Matrix<double, DIM, N_CORE_DOFS>::Zero();
        std::array<Eigen::Matrix<double, N_CORE_DOFS, N_CORE_DOFS>, DIM> H;
        for (int i = 0; i < DIM; i++) {
            out(i) = d(i).val;
            J.row(i) = d(i).grad;
            H[i] = d(i).Hess;
        }
        return std::make_tuple(out, J, H);
    }

    static GradType<N_CORE_DOFS + 1> compute_mollifier_gradient(
        const Vector<double, N_CORE_DOFS>& x, const double dist_sqr)
    {
        using T = TinyADGrad<N_CORE_DOFS + 1>;
        const Vector<T, N_CORE_DOFS + 1> X =
            T::make_active(
                (Vector<double, N_CORE_DOFS + 1>() << x, dist_sqr).finished());
        const T out =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                X.head(N_CORE_DOFS), X(N_CORE_DOFS));

        return std::make_tuple(out.val, out.grad);
    }

    static HessianType<N_CORE_DOFS + 1> compute_mollifier_hessian(
        const Vector<double, N_CORE_DOFS>& x, const double dist_sqr)
    {
        using T = TinyADHessian<N_CORE_DOFS + 1>;
        const Vector<T, N_CORE_DOFS + 1> X =
            T::make_active(
                (Vector<double, N_CORE_DOFS + 1>() << x, dist_sqr).finished());
        const T out =
            PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::mollifier(
                X.head(N_CORE_DOFS), X(N_CORE_DOFS));

        return std::make_tuple(
            out.val, out.grad, out.Hess);
    }
};

} // namespace ipc

#include "primitive_distance.tpp"
