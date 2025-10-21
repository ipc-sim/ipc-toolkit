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
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const SmoothContactParameters& params);

    Point3(
        const index_t id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices);
    virtual ~Point3() = default;

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * DIM; }

    // assume the following functions are only called if active
    double potential(
        const Vector<double, DIM>& d,
        const Vector<double, -1, MAX_SIZE>& x) const;
    // derivatives including wrt. d (the closest direction) in front
    Vector<double, -1, MAX_SIZE + DIM> grad(
        const Vector<double, DIM>& d,
        const Vector<double, -1, MAX_SIZE>& x) const;
    MatrixMax<double, MAX_SIZE + DIM, MAX_SIZE + DIM> hessian(
        const Vector<double, DIM>& d,
        const Vector<double, -1, MAX_SIZE>& x) const;

    /// @brief
    /// @tparam scalar
    /// @param direc normalized
    /// @param v
    /// @param direc points from v to the other point
    /// @param neighbors follow counter-clockwise order
    /// @param params
    /// @return
    template <typename scalar, int n_verts = -1>
    scalar smooth_point3_term(
        const Eigen::Matrix<scalar, n_verts, 3>& X,
        Eigen::ConstRef<Eigen::RowVector3<scalar>> direc) const;

    GradType<-1> smooth_point3_term_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> X,
        const SmoothContactParameters& params) const;

    HessianType<-1> smooth_point3_term_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> X,
        const SmoothContactParameters& params) const;

    GradType<-1> smooth_point3_term_tangent_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
        const double alpha,
        const double beta) const;

    HessianType<-1> smooth_point3_term_tangent_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
        const double alpha,
        const double beta) const;

    GradType<-1> smooth_point3_term_normal_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
        const double alpha,
        const double beta) const;

    HessianType<-1> smooth_point3_term_normal_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direc,
        Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
        const double alpha,
        const double beta) const;

private:
    int n_neighbors;
    OrientationTypes otypes;

    std::vector<index_t> local_to_global_vids;
    std::map<index_t, int> global_to_local_vids;

    Eigen::Matrix<int, -1, 3> faces;
    Eigen::Matrix<int, -1, 2> edges;
    bool orientable;

    bool smooth_point3_term_type(
        const Eigen::Matrix<double, -1, 3>& X,
        Eigen::ConstRef<Eigen::RowVector3d> direc);
};

} // namespace ipc
