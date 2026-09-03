#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/high_order_contact/arbitrary_point_bvh.hpp>
#include <ipc/high_order_contact/collisions/high_order_collision_dict.hpp>
#include <ipc/high_order_contact/high_order_contact_parameters.hpp>

#include <memory>
#include <tuple>

namespace ipc {

/// @brief Evaluate the high-order contact potential (ESP) at an arbitrary
/// point in space, not restricted to the mesh's own vertices/edges/faces.
///
/// Reuses the same collision-construction machinery as
/// PointPotential::build_collisions_at_vertex / build_collisions_at_edge_qp
/// (quadrature_potential.cpp): nearby primitives found by ArbitraryPointBVH
/// are classified with exact point-primitive distance-type predicates --
/// which of a primitive's sub-features (face/edge/vertex in 3D,
/// edge/vertex in 2D) the query's closest point actually falls on -- and
/// redundant contributions from different primitives resolving to the same
/// feature (e.g. two adjacent faces and their shared edge, or the two 2D
/// edges meeting at a corner) are merged and symbolically cancelled
/// (integer weight accumulation, dropped exactly at zero) before any
/// barrier value is computed. This avoids the catastrophic cancellation a
/// naive "sum every found primitive independently with a fixed sign"
/// approach suffers near shared features, where barrier derivatives blow up
/// as distance -> 0.
///
/// Value/gradient/Hessian are computed by the same HighOrderCollision::
/// operator()/gradient()/hessian() used by the production
/// HighOrderContactPotential, just evaluated for a virtual point
/// (id == V.rows()) instead of a real mesh vertex, via VertexMatrixView.
///
/// A single fixed params.dhat is used everywhere; AdaptiveSupport
/// (per-primitive dhat) is not supported.
///
/// @tparam dim Spatial dimension of the mesh, 2 or 3.
template <int dim> class ArbitraryPointPotential {
    static_assert(dim == 2 || dim == 3, "dim must be 2 or 3");

public:
    /// @brief A query point in space.
    using Point = Eigen::RowVector<double, dim>;
    /// @brief Gradient of the potential w.r.t. a query point.
    using Gradient = Eigen::Vector<double, dim>;
    /// @brief Hessian of the potential w.r.t. a query point.
    using Hessian = Eigen::Matrix<double, dim, dim>;

    /// @throws std::runtime_error if mesh.dim() != dim.
    ArbitraryPointPotential(
        const CollisionMesh& mesh, HighOrderContactParameters params);

    /// @brief Rebuild the underlying broad-phase index. O(n log n). Call
    /// once per vertex configuration, before any operator()/gradient()/
    /// hessian() calls against that configuration.
    void update(Eigen::ConstRef<Eigen::MatrixXd> V);

    /// @brief Evaluate the potential at q.
    double operator()(
        Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const;

    /// @brief Gradient of the potential with respect to q.
    Gradient gradient(
        Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const;

    /// @brief Hessian of the potential with respect to q.
    Hessian
    hessian(Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const;

    /// @brief Value, gradient, and Hessian at q, computed together.
    ///
    /// Equivalent to calling operator()/gradient()/hessian() separately, but
    /// builds the (BVH-queried, exact-predicate-classified,
    /// symbolically-cancelled) collision dict for q only once instead of
    /// three times -- collision construction, not the final per-collision
    /// barrier evaluation, dominates cost (see the profiling that motivated
    /// this), so calling operator()/gradient()/hessian() separately at the
    /// same point does ~3x the necessary work. Prefer this whenever you need
    /// more than one of the three at the same q (e.g. a Newton step).
    std::tuple<double, Gradient, Hessian> evaluate(
        Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const;

private:
    /// @brief Build the (symbolically-cancelled) collision dict for q,
    /// mirroring PointPotential::build_collisions_at_vertex (3D) /
    /// build_collisions_at_edge_qp (2D) with q as a virtual vertex
    /// (id == V.rows()) instead of a real one or an edge quadrature point,
    /// and candidates sourced from point_bvh instead of
    /// Candidates::vv_set/ve_set/vf_set.
    std::unique_ptr<HighOrderCollisionDict<PointType::VERTEX, dim>>
    build_collisions_at_point(
        Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const;

    const CollisionMesh& mesh;
    HighOrderContactParameters params;
    ArbitraryPointBVH point_bvh;
};

extern template class ArbitraryPointPotential<2>;
extern template class ArbitraryPointPotential<3>;

} // namespace ipc
