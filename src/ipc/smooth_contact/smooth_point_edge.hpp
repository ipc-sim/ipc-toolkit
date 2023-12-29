#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include "smooth_collision_constraint.hpp"
#include <ipc/distance/distance_type.hpp>

namespace ipc {

class SmoothEdgeVertexConstraint : public EdgeVertexCandidate,
                                   public SmoothCollisionConstraint {
public:
    using EdgeVertexCandidate::EdgeVertexCandidate;

    SmoothEdgeVertexConstraint(
        const long _edge_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : EdgeVertexCandidate(_edge_id, _vertex_id)
        , SmoothCollisionConstraint(_weight, _weight_gradient)
    {
    }

    template <typename H>
    friend H AbslHashValue(H h, const SmoothEdgeVertexConstraint& ev)
    {
        return AbslHashValue(
            std::move(h), static_cast<const EdgeVertexCandidate&>(ev));
    }

    double compute_potential(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const SmoothParameter param) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const SmoothParameter param) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const SmoothParameter param,
        const bool project_hessian_to_psd) const override;
};

    /// @brief Compute pointwise potential for a point p and a point specified by uv on edge [e0, e1]
    /// @param p One point outside of edge
    /// @param e0 One end point of the edge
    /// @param e1 One end point of the edge
    /// @param uv Barycentric coordinate
    /// @param dhat The effective distance of barrier
    /// @param alpha The effective angle of barrier
    template <typename scalar>
    scalar smooth_point_edge_potential_pointwise(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &uv,
        const double &dhat,
        const double &alpha);

    /// @brief Compute potential for a point p and an edge [e0, e1], integrated over the edge with high order quadrature
    template <typename scalar>
    scalar smooth_point_edge_potential_quadrature(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);

    /// @brief Compute potential for a list of points and an edge [e0, e1], integrated over the edge with high order quadrature
    template <typename scalar>
    Vector<scalar, -1, -1> smooth_point_edge_potentials_quadrature(
        const Eigen::Ref<Eigen::Matrix<scalar, -1, -1, Eigen::RowMajor, -1, 3>>& points,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const int &N);

    /// @brief Compute potential for a point p and an edge [e0, e1], using the smooth closest point
    template <typename scalar>
    scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const double &dhat,
        const double &alpha,
        const double &r);
}