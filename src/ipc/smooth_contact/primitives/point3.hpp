#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Point3 : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 1;
    static constexpr int DIM = 3;
    static constexpr int MAX_SIZE = N_VERT_NEIGHBORS_3D * DIM;
    // d is a vector from this point to the other primitive
    Point3(
        const index_t id,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<VectorMax3d> d,
        const SmoothContactParameters& params);

    Point3(
        const index_t id,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices);
    virtual ~Point3() = default;

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    // assume the following functions are only called if active
    double potential(
        const Eigen::Vector<double, DIM>& d,
        const VectorMax<double, MAX_SIZE>& x) const;
    // derivatives including wrt. d (the closest direction) in front
    VectorMax<double, MAX_SIZE + DIM> grad(
        const Eigen::Vector<double, DIM>& d,
        const VectorMax<double, MAX_SIZE>& x) const;
    MatrixMax<double, MAX_SIZE + DIM, MAX_SIZE + DIM> hessian(
        const Eigen::Vector<double, DIM>& d,
        const VectorMax<double, MAX_SIZE>& x) const;

    /// @brief Compute the smooth point term for this vertex.
    /// @tparam scalar The scalar type.
    /// @tparam n_verts The compile-time row count of X, or -1 if dynamic.
    /// @param X Local vertex positions, one per row: this vertex first,
    ///     followed by its one-ring neighbors in counter-clockwise order.
    /// @param direc Direction pointing from this vertex to the other point.
    ///     It is normalized internally, so it need not arrive normalized.
    /// @return The product of the weight, normal, and tangent terms.
    template <typename scalar, int n_verts = -1>
    scalar smooth_point3_term(
        const Eigen::Matrix<scalar, n_verts, 3>& X,
        Eigen::ConstRef<Eigen::RowVector3<scalar>> direc) const;

    GradientType<-1> smooth_point3_term_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::MatrixX3d> X,
        const SmoothContactParameters& params) const;

    HessianType<-1> smooth_point3_term_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::MatrixX3d> X,
        const SmoothContactParameters& params) const;

    GradientType<-1> smooth_point3_term_tangent_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::MatrixX3d> tangents,
        const double alpha,
        const double beta) const;

    HessianType<-1> smooth_point3_term_tangent_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::MatrixX3d> tangents,
        const double alpha,
        const double beta) const;

    GradientType<-1> smooth_point3_term_normal_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::MatrixX3d> tangents,
        const double alpha,
        const double beta) const;

    HessianType<-1> smooth_point3_term_normal_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::MatrixX3d> tangents,
        const double alpha,
        const double beta) const;

private:
    int n_neighbors;
    OrientationTypes otypes;

    std::vector<index_t> local_to_global_vids;
    std::map<index_t, int> global_to_local_vids;

    Eigen::MatrixX3i faces;
    Eigen::MatrixX2i edges;
    bool orientable;

    bool smooth_point3_term_type(
        Eigen::ConstRef<Eigen::MatrixX3d> X,
        Eigen::ConstRef<Eigen::RowVector3d> direc);
};

} // namespace ipc
