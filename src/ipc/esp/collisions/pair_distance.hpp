#pragma once
#include "esp_primitives.hpp"

#include <ipc/distance/distance_type.hpp>

namespace ipc {
template <typename PrimitiveA, typename PrimitiveB> struct PairDistType { };

template <> struct PairDistType<Vertex3, Edge3P1> {
    using type = PointEdgeDistanceType;
    static constexpr std::string_view NAME = "PointEdge";
};

template <> struct PairDistType<Edge3P1, Edge3P1> {
    using type = EdgeEdgeDistanceType;
    static constexpr std::string_view NAME = "EdgeEdge";
};

template <> struct PairDistType<Vertex3, Vertex3> {
    using type = PointPointDistanceType;
    static constexpr std::string_view NAME = "PointPoint";
};

template <> struct PairDistType<Vertex3, Face3P1> {
    using type = PointTriangleDistanceType;
    static constexpr std::string_view NAME = "PointFace";
};

template <typename T, typename PrimitiveA, typename PrimitiveB>
class PairDistance {
public:
    static_assert(
        PrimitiveA::DIM == PrimitiveB::DIM,
        "Primitives must have the same dimension");
    static constexpr int DIM = PrimitiveA::DIM;
    static constexpr int N_DOFS = PrimitiveA::N_POINTS * PrimitiveA::DIM
        + PrimitiveB::N_POINTS * PrimitiveB::DIM;
    static typename PairDistType<PrimitiveA, PrimitiveB>::type
    compute_distance_type(Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X);
    static T compute_distance(
        Eigen::ConstRef<Eigen::Vector<T, N_DOFS>> X,
        typename PairDistType<PrimitiveA, PrimitiveB>::type dtype);
};
} // namespace ipc

#include "pair_distance.tpp"
