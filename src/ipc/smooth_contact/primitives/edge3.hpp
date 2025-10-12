#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Edge3 : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int DIM = 3;
    // d is a vector from closest point on the edge to the point outside of the
    // edge
    Edge3(
        const index_t id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const SmoothContactParameters& params);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    double potential(
        Eigen::ConstRef<Eigen::Vector3d> d, Eigen::ConstRef<Vector12d> x) const;
    Vector15d grad(
        Eigen::ConstRef<Eigen::Vector3d> d, Eigen::ConstRef<Vector12d> x) const;
    Matrix15d hessian(
        Eigen::ConstRef<Eigen::Vector3d> d, Eigen::ConstRef<Vector12d> x) const;

private:
    OrientationTypes otypes;

    bool has_neighbor_1, has_neighbor_2;
    bool orientable;
};

double smooth_edge3_normal_term(
    Eigen::ConstRef<Vector3d> dn,
    Eigen::ConstRef<Vector3d> e0,
    Eigen::ConstRef<Vector3d> e1,
    Eigen::ConstRef<Vector3d> f0,
    Eigen::ConstRef<Vector3d> f1,
    const double alpha,
    const double beta,
    const OrientationTypes& otypes);

GradType<15> smooth_edge3_normal_term_gradient(
    Eigen::ConstRef<Vector3d> dn,
    Eigen::ConstRef<Vector3d> e0,
    Eigen::ConstRef<Vector3d> e1,
    Eigen::ConstRef<Vector3d> f0,
    Eigen::ConstRef<Vector3d> f1,
    const double alpha,
    const double beta,
    const OrientationTypes& otypes);

HessianType<15> smooth_edge3_normal_term_hessian(
    Eigen::ConstRef<Vector3d> dn,
    Eigen::ConstRef<Vector3d> e0,
    Eigen::ConstRef<Vector3d> e1,
    Eigen::ConstRef<Vector3d> f0,
    Eigen::ConstRef<Vector3d> f1,
    const double alpha,
    const double beta,
    const OrientationTypes& otypes);
} // namespace ipc