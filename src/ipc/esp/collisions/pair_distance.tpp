#include "pair_distance.hpp"

#include <ipc/distance/distance_type_exact.hpp>
#include <ipc/gcp/distance/edge_edge.hpp>
#include <ipc/gcp/distance/point_face.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
template <typename T> class PairDistance<Edge3P1, Edge3P1, T> {
public:
    static_assert(
        Edge3P1::DIM == Edge3P1::DIM,
        "Primitives must have the same dimension");
    static constexpr int DIM = Edge3P1::DIM;
    static constexpr int N_DOFS =
        Edge3P1::N_POINTS * Edge3P1::DIM + Edge3P1::N_POINTS * Edge3P1::DIM;
    static PairDistType<Edge3P1, Edge3P1>::type
    compute_distance_type(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X)
    {
        if constexpr (std::is_same_v<T, double>)
            return edge_edge_distance_type_exact(
                X.template head<3>() /* edge 0 */,
                X.template segment<3>(3) /* edge 0 */,
                X.template segment<3>(6) /* edge 1 */,
                X.template segment<3>(9) /* edge 1 */);
        else
            return EdgeEdgeDistanceType::AUTO;
    }
    static T compute_distance(
        Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X,
        PairDistType<Edge3P1, Edge3P1>::type dtype)
    {
        return edge_edge_sqr_distance<T>(
            X.template head<3>() /* edge 0 */,
            X.template segment<3>(3) /* edge 0 */,
            X.template segment<3>(6) /* edge 1 */,
            X.template tail<3>() /* edge 1 */, dtype);
    }
};

// template <typename T>
// class PairDistance<Edge3P1, Vertex3, T> {
// public:
//     static_assert(
//       Edge3P1::DIM == Vertex3::DIM,
//       "Primitives must have the same dimension");
//     static constexpr int DIM = Edge3P1::DIM;
//     static constexpr int N_DOFS =
//         Edge3P1::N_POINTS * Edge3P1::DIM
//         + Vertex3::N_POINTS * Vertex3::DIM;
//     static T compute_distance(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X)
//     {
//         return PointEdgeDistance<T, DIM>::point_edge_sqr_distance(
//             X.template tail<3>(),
//             X.template head<3>(),
//             X.template segment<3>(3));
//     }
// };

// template <typename T>
// class PairDistance<Edge3P1, Face3P1, T> {
// public:
//     static_assert(
//       Edge3P1::DIM == Face3P1::DIM,
//       "Primitives must have the same dimension");
//     static constexpr int DIM = Edge3P1::DIM;
//     static constexpr int N_DOFS =
//         Edge3P1::N_POINTS * Edge3P1::DIM
//         + Face3P1::N_POINTS * Face3P1::DIM;
//     static T compute_distance(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X)
//     {
//         Eigen::ConstRef<Eigen::Vector<T, DIM>> e0 = X.template head<3>();
//         Eigen::ConstRef<Eigen::Vector<T, DIM>> e1 = X.template segment<3>(3);
//
//         Eigen::ConstRef<Eigen::Vector<T, DIM>> f0 = X.template segment<3>(6);
//         Eigen::ConstRef<Eigen::Vector<T, DIM>> f1 = X.template segment<3>(9);
//         Eigen::ConstRef<Eigen::Vector<T, DIM>> f2 = X.template
//         segment<3>(12);
//
//         return std::min({
//             point_triangle_sqr_distance<T>(e0, f0, f1, f2),
//             point_triangle_sqr_distance<T>(e1, f0, f1, f2),
//             edge_edge_sqr_distance<T>(e0, e1, f0, f1),
//             edge_edge_sqr_distance<T>(e0, e1, f2, f1),
//             edge_edge_sqr_distance<T>(e0, e1, f0, f2)});
//     }
// };

template <typename T> class PairDistance<Vertex3, Edge3P1, T> {
public:
    static_assert(
        Vertex3::DIM == Edge3P1::DIM,
        "Primitives must have the same dimension");
    static constexpr int DIM = Vertex3::DIM;
    static constexpr int N_DOFS =
        Vertex3::N_POINTS * Vertex3::DIM + Edge3P1::N_POINTS * Edge3P1::DIM;
    static PairDistType<Vertex3, Edge3P1>::type
    compute_distance_type(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X)
    {
        if constexpr (std::is_same_v<T, double>)
            return point_edge_distance_type_exact(
                X.template head<3>(), X.template segment<3>(3),
                X.template segment<3>(6));
        else
            return PointEdgeDistanceType::AUTO;
    }
    static T compute_distance(
        Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X,
        PairDistType<Vertex3, Edge3P1>::type dtype)
    {
        return PointEdgeDistance<T, DIM>::point_edge_sqr_distance(
            X.template head<3>(), X.template segment<3>(3),
            X.template segment<3>(6), dtype);
    }
};

template <typename T> class PairDistance<Vertex3, Vertex3, T> {
public:
    static_assert(
        Vertex3::DIM == Vertex3::DIM,
        "Primitives must have the same dimension");
    static constexpr int DIM = Vertex3::DIM;
    static constexpr int N_DOFS =
        Vertex3::N_POINTS * Vertex3::DIM + Vertex3::N_POINTS * Vertex3::DIM;
    static PairDistType<Vertex3, Vertex3>::type
    compute_distance_type(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X)
    {
        return PointPointDistanceType::AUTO;
    }
    static T compute_distance(
        Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X,
        PairDistType<Vertex3, Vertex3>::type dtype)
    {
        return (X.template head<3>() - X.template tail<3>()).squaredNorm();
    }
};

template <typename T> class PairDistance<Vertex3, Face3P1, T> {
public:
    static_assert(
        Vertex3::DIM == Face3P1::DIM,
        "Primitives must have the same dimension");
    static constexpr int DIM = Vertex3::DIM;
    static constexpr int N_DOFS =
        Vertex3::N_POINTS * Vertex3::DIM + Face3P1::N_POINTS * Face3P1::DIM;
    static PairDistType<Vertex3, Face3P1>::type
    compute_distance_type(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X)
    {
        Eigen::ConstRef<Eigen::Vector<T, DIM>> v = X.template head<3>();

        Eigen::ConstRef<Eigen::Vector<T, DIM>> f0 = X.template segment<3>(3);
        Eigen::ConstRef<Eigen::Vector<T, DIM>> f1 = X.template segment<3>(6);
        Eigen::ConstRef<Eigen::Vector<T, DIM>> f2 = X.template segment<3>(9);
        if constexpr (std::is_same_v<T, double>)
            return point_triangle_distance_type_exact(v, f0, f1, f2);
        else
            return PointTriangleDistanceType::AUTO;
    }
    static T compute_distance(
        Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X,
        PairDistType<Vertex3, Face3P1>::type dtype)
    {
        Eigen::ConstRef<Eigen::Vector<T, DIM>> v = X.template head<3>();

        Eigen::ConstRef<Eigen::Vector<T, DIM>> f0 = X.template segment<3>(3);
        Eigen::ConstRef<Eigen::Vector<T, DIM>> f1 = X.template segment<3>(6);
        Eigen::ConstRef<Eigen::Vector<T, DIM>> f2 = X.template segment<3>(9);

        return point_triangle_sqr_distance<T>(v, f0, f1, f2, dtype);
    }
};
} // namespace ipc
