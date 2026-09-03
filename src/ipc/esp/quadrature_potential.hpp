#pragma once
#include "ipc/candidates/edge_edge.hpp"
#include "ipc/collision_mesh.hpp"
#include "ipc/distance/point_point.hpp"
#include "ipc/distance/point_triangle.hpp"
#include "ipc/esp/esp_collisions.hpp"
#include "ipc/gcp/distance/edge_edge.hpp"

#include <array>

namespace ipc {
namespace PointPotentialHelper {
    double evaluate_potential_at_vertex_with_cached_collisions(
        const Eigen::MatrixXd& V,
        const ESPCollisionDict<PointType::VERTEX>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive);

    Eigen::VectorXd
    evaluate_potential_gradient_at_vertex_with_cached_collisions(
        const Eigen::MatrixXd& V,
        const ESPCollisionDict<PointType::VERTEX>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive);

    Eigen::MatrixXd evaluate_potential_hessian_at_vertex_with_cached_collisions(
        const Eigen::MatrixXd& V,
        const ESPCollisionDict<PointType::VERTEX>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        PSDProjectionMethod project_to_psd);

    std::pair<double, double>
    evaluate_potential_at_vertex_with_cached_collisions_nearfar(
        const Eigen::MatrixXd& V,
        const ESPCollisionDict<PointType::VERTEX>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const NearFarBarrier& nf_barrier);

    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
        const Eigen::MatrixXd& V,
        const ESPCollisionDict<PointType::VERTEX>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const NearFarBarrier& nf_barrier);

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
        const Eigen::MatrixXd& V,
        const ESPCollisionDict<PointType::VERTEX>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        PSDProjectionMethod project_to_psd,
        const NearFarBarrier& nf_barrier);

    double evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::EDGE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        EdgeEdgeDistanceType dtype);

    double
    evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::EDGE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        EdgeEdgeDistanceType dtype,
        const NearFarBarrier& nf_barrier);

    /// @brief Compute the gradient of P(q) for a point q
    /// @return The gradient vector with respect to collisions.m_vertex_ids
    /// @param V_extended Extended vertices matrix, with the point q appended to the last row
    /// @param collisions Primitives that are close in distance to point q
    /// @param q Closest point between two edges, together with the derivatives of q with respect to vids
    template <typename ADType>
    std::enable_if_t<
        IsADGrad<ADType>::value || IsADHessian<ADType>::value,
        Eigen::VectorXd>
    evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::EDGE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        Eigen::ConstRef<Eigen::Vector3<ADType>> q);

    /// @brief Compute the near component gradient (NearFarBarrier)
    template <typename ADType>
    std::enable_if_t<
        IsADGrad<ADType>::value || IsADHessian<ADType>::value,
        Eigen::VectorXd>
    evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::EDGE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        Eigen::ConstRef<Eigen::Vector3<ADType>> q,
        const NearFarBarrier& nf_barrier);

    Eigen::MatrixXd
    evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::EDGE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> q);

    Eigen::MatrixXd
    evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::EDGE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        Eigen::ConstRef<Eigen::Vector3<ADHessian<12>>> q,
        const NearFarBarrier& nf_barrier);

    double evaluate_potential_at_face_center_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive);

    std::pair<double, double>
    evaluate_potential_at_face_center_with_cached_collisions_nearfar(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const NearFarBarrier& nf_barrier);

    Eigen::VectorXd
    evaluate_potential_gradient_at_face_center_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive);

    Eigen::MatrixXd
    evaluate_potential_hessian_at_face_center_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        PSDProjectionMethod project_to_psd);

    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    evaluate_potential_gradient_at_face_center_with_cached_collisions_nearfar(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const NearFarBarrier& nf_barrier);

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    evaluate_potential_hessian_at_face_center_with_cached_collisions_nearfar(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        PSDProjectionMethod project_to_psd,
        const NearFarBarrier& nf_barrier);

    /// @brief Gradient of the face-interior potential for an arbitrary
    ///   interior quadrature point q = λ0·v0 + λ1·v1 + λ2·v2.
    /// @param lambda Barycentric coordinates of the interior point.
    ///   The chain-rule factors λk replace the 1/3 used for the centroid.
    Eigen::VectorXd
    evaluate_potential_gradient_at_face_interior_point_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const std::array<double, 3>& lambda);

    /// @brief Hessian of the face-interior potential for an arbitrary
    ///   interior quadrature point q = λ0·v0 + λ1·v1 + λ2·v2.
    Eigen::MatrixXd
    evaluate_potential_hessian_at_face_interior_point_with_cached_collisions(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const std::array<double, 3>& lambda,
        PSDProjectionMethod project_to_psd);

    std::pair<Eigen::VectorXd, Eigen::VectorXd>
    evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const std::array<double, 3>& lambda,
        const NearFarBarrier& nf_barrier);

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    evaluate_potential_hessian_at_face_interior_point_with_cached_collisions_nearfar(
        VertexMatrixView<3> V_extended,
        const ESPCollisionDict<PointType::FACE>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const std::array<double, 3>& lambda,
        PSDProjectionMethod project_to_psd,
        const NearFarBarrier& nf_barrier);

    // ---- 2D edge quadrature point helpers ----

    /// @brief Evaluate P(q) = sum of barrier values for all pairs in the dict.
    /// @param V_extended Vertices extended with the virtual QP as last row.
    /// @param dict Per-QP collision dict for edge quadrature.
    /// @param params Contact parameters.
    double evaluate_potential_at_edge_qp(
        VertexMatrixView<2> V_extended,
        const ESPCollisionDict<PointType::EDGE, 2>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive);

    /// @brief Gradient of P(q) w.r.t. all real vertices, using chain rule
    /// dP/de_k += lambda[k] * dP/dq.
    /// @param lambda Barycentric coords of QP on the edge: q = lambda[0]*e0 + lambda[1]*e1.
    Eigen::VectorXd evaluate_potential_gradient_at_edge_qp(
        VertexMatrixView<2> V_extended,
        const ESPCollisionDict<PointType::EDGE, 2>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const std::array<double, 2>& lambda);

    /// @brief Hessian of P(q) w.r.t. all real vertices.
    /// @param lambda Barycentric coords of QP on the edge: q = lambda[0]*e0 + lambda[1]*e1.
    Eigen::MatrixXd evaluate_potential_hessian_at_edge_qp(
        VertexMatrixView<2> V_extended,
        const ESPCollisionDict<PointType::EDGE, 2>& collisions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const std::array<double, 2>& lambda,
        PSDProjectionMethod project_to_psd);
} // namespace PointPotentialHelper

class PointPotential {
public:
    constexpr static int r = 2;

    PointPotential(
        const CollisionMesh& mesh_,
        const Candidates& candidates_,
        const ESPParameters params_,
        const AdaptiveSupport* adaptive_ = nullptr)
        : mesh(mesh_)
        , candidates(candidates_)
        , params(params_)
        , adaptive(adaptive_)
    {
    }

    std::unique_ptr<ESPCollisionDict<PointType::VERTEX>>
    build_collisions_at_vertex(
        const Eigen::MatrixXd& V,
        index_t vid,
        size_t& num_collision_pairs) const;

    std::unique_ptr<ESPCollisionDict<PointType::EDGE>>
    build_collisions_at_edge_edge_closest_point(
        const Eigen::MatrixXd& V,
        index_t e0,
        index_t e1,
        EdgeEdgeDistanceType dtype,
        size_t& num_collision_pairs) const;

    std::unique_ptr<ESPCollisionDict<PointType::FACE>>
    build_collisions_at_face_center(
        const Eigen::MatrixXd& V,
        index_t fid,
        size_t& num_collision_pairs) const;

    std::unique_ptr<ESPCollisionDict<PointType::FACE>>
    build_collisions_at_face_interior_point(
        const Eigen::MatrixXd& V,
        index_t fid,
        const std::array<double, 3>& lambda,
        size_t& num_collision_pairs) const;

    /// @brief Build a per-QP collision dict for a 2D edge quadrature point.
    /// @param V Vertex positions (2D).
    /// @param ei Source edge index.
    /// @param lambda Barycentric coords of QP: q = lambda[0]*e0 + lambda[1]*e1.
    /// @param dhat Distance threshold for this edge.
    std::unique_ptr<ESPCollisionDict<PointType::EDGE, 2>>
    build_collisions_at_edge_qp(
        const Eigen::MatrixXd& V,
        index_t ei,
        const std::array<double, 2>& lambda,
        double dhat,
        size_t& num_collision_pairs) const;

    const CollisionMesh& mesh;
    const Candidates& candidates;
    const ESPParameters params;
    const AdaptiveSupport* adaptive;
};
} // namespace ipc
