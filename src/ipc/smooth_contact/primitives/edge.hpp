#pragma once

#include "primitive.hpp"

namespace ipc {

template <int DIM> class Edge : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    using DVector = Eigen::Vector<double, DIM>;
    using XVector = Eigen::Vector<double, DIM == 2 ? 4 : 12>;
    using GradType = Eigen::Vector<double, DIM == 2 ? 6 : 15>;
    static constexpr int HESSIAN_ROWS = DIM == 2 ? 6 : 15;
    static constexpr int HESSIAN_COLS = HESSIAN_ROWS;
    using HessianType = Eigen::Matrix<double, HESSIAN_ROWS, HESSIAN_COLS>;

    // d is a vector from closest point on the edge to the point outside of the
    // edge
    Edge(
        const index_t id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const SmoothContactParameters& params);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    /// @brief ???
    /// @param d Vector from closest point on the edge to the point outside of the edge
    /// @param x ???
    /// @return ???
    double
    potential(Eigen::ConstRef<DVector> d, Eigen::ConstRef<XVector> x) const;

    /// @brief ???
    /// @param d Vector from closest point on the edge to the point outside of the edge
    /// @param x ???
    /// @return ???
    GradType grad(Eigen::ConstRef<DVector> d, Eigen::ConstRef<XVector> x) const;

    /// @brief ???
    /// @param d Vector from closest point on the edge to the point outside of the edge
    /// @param x ???
    /// @return ???
    HessianType
    hessian(Eigen::ConstRef<DVector> d, Eigen::ConstRef<XVector> x) const;
};

using Edge2 = Edge<2>;
using Edge3 = Edge<3>;

} // namespace ipc
