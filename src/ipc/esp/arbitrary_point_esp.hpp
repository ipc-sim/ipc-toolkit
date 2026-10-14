#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/esp/arbitrary_point_bvh.hpp>
#include <ipc/esp/collisions/esp_collision_dict.hpp>
#include <ipc/esp/esp_parameters.hpp>

#include <memory>
#include <tuple>
#include <vector>

namespace ipc {

/// @brief Evaluate the Extremum Sum Potential (ESP) at an arbitrary
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
/// Every element enters with its inclusion-exclusion weight, chosen so that
/// the weights of the elements containing any point of the mesh sum to one
/// (supplemental S2). A closed mesh gets the alternating signs; boundary edges
/// and vertices get zero (S4), and edges in no face and isolated vertices get
/// one, so open meshes and lower-dimensional pieces are handled too.
///
/// Value/gradient/Hessian are computed by the same ESPCollision::
/// operator()/gradient()/hessian() used by the production
/// ESPPotential, just evaluated for a virtual point
/// (id == V.rows()) instead of a real mesh vertex, via VertexMatrixView.
///
/// @tparam dim Spatial dimension of the mesh, 2 or 3.
template <int dim> class ArbitraryPointESP {
    static_assert(dim == 2 || dim == 3, "dim must be 2 or 3");

public:
    /// @brief A query point in space.
    using Point = Eigen::RowVector<double, dim>;
    /// @brief Gradient of the potential w.r.t. a query point.
    using Gradient = Eigen::Vector<double, dim>;
    /// @brief Hessian of the potential w.r.t. a query point.
    using Hessian = Eigen::Matrix<double, dim, dim>;

    /// @throws std::runtime_error if mesh.dim() != dim.
    ArbitraryPointESP(const CollisionMesh& mesh, ESPParameters params);

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
    const CollisionMesh& m_mesh;
    ESPParameters m_params;
    ArbitraryPointBVH m_point_bvh;
    /// Inclusion-exclusion weight of each edge and each vertex (every face
    /// has weight one), from the mesh connectivity; see the constructor.
    std::vector<int> m_edge_weights, m_vertex_weights;
};

extern template class ArbitraryPointESP<2>;
extern template class ArbitraryPointESP<3>;

} // namespace ipc
