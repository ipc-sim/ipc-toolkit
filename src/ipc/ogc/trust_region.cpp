#include "trust_region.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/tangent/closest_point.hpp>

namespace ipc::ogc {

TrustRegion::TrustRegion(double dhat) : trust_region_inflation_radius(2 * dhat)
{
    assert(dhat > 0);
}

void TrustRegion::warm_start_time_step(
    const CollisionMesh& mesh,
    Eigen::Ref<Eigen::MatrixXd> x,
    Eigen::ConstRef<Eigen::MatrixXd> pred_x,
    NormalCollisions& collisions,
    const double dhat,
    const double min_distance,
    BroadPhase* broad_phase)
{
    // Assign a default broad phase if none is provided.
    std::unique_ptr<BroadPhase> default_broad_phase;
    if (broad_phase == nullptr) {
        default_broad_phase = make_default_broad_phase();
        broad_phase = default_broad_phase.get();
    }

    assert(x.rows() == pred_x.rows() && x.cols() == pred_x.cols());

    // Compute the norm of the translation for each vertex.
    Eigen::VectorXd dx_norm = (pred_x - x).rowwise().norm();
    assert(dx_norm.size() == x.rows());

    // Compute a trust region inflation radius based on the distance
    // between the current positions x and the predicted positions x̂.
    trust_region_inflation_radius = std::max(2 * dhat, dx_norm.maxCoeff());

    // Initialize the trust region around x
    update(mesh, x, collisions, min_distance, broad_phase);
    should_update_trust_region = false;

    int num_updates = 0;
    for (int i = 0; i < x.rows(); i++) {
        const VectorMax3d x_i = x.row(i);
        const VectorMax3d pred_x_i = pred_x.row(i);

        if (dx_norm(i) <= trust_region_radii(i)) {
            x.row(i) = pred_x_i; // Free to move to the predicted position
        } else {
            // Move to the boundary of the trust region
            x.row(i) =
                x_i + (trust_region_radii(i) / dx_norm(i)) * (pred_x_i - x_i);
            ++num_updates;
        }
    }

    // If a critical mass of vertices are restricted by the trust region,
    // we update the trust region to avoid excessive filtering.
    if (num_updates > update_threshold * mesh.num_vertices()) {
        update(mesh, x, collisions, min_distance, broad_phase);
    }
}

void TrustRegion::update(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> x,
    NormalCollisions& collisions,
    const double min_distance,
    BroadPhase* broad_phase)
{
    assert(0 < relaxed_radius_scaling && relaxed_radius_scaling < 1);

    // Save the x for trust region filtering in filter_step.
    trust_region_centers = x;

    Eigen::MatrixXd vertices;
    if (x.rows() == mesh.num_vertices()) {
        vertices = x;
    } else {
        vertices = mesh.vertices(x);
    }

    Candidates trust_region_candidates;
    trust_region_candidates.build(
        mesh, vertices, trust_region_inflation_radius + min_distance,
        broad_phase);

    // Compute the trust region distances for the reduced set of collision
    // vertices.
    Eigen::VectorXd reduced_trust_region_radii =
        trust_region_candidates.compute_per_vertex_safe_distances(
            mesh, vertices, trust_region_inflation_radius, min_distance);

    // Use < half the safe distance to account for double sided contact
    reduced_trust_region_radii *= relaxed_radius_scaling / 2;

    // Copy the reduced trust region distances to the full set of vertices.
    // Set interior vertices to infinity, so they are not restricted.
    if (x.rows() == mesh.num_vertices()) {
        trust_region_radii = reduced_trust_region_radii;
    } else {
        trust_region_radii.setConstant(
            x.rows(), std::numeric_limits<double>::infinity());
        trust_region_radii(mesh.to_full_vertex_id()) =
            reduced_trust_region_radii;
    }

    // Update collisions using the trust region candidates.
    // NOTE: Use the trust_region_inflation_radius to ensure that the
    // collisions include all pairs until this function is called again.
    collisions.build(
        trust_region_candidates, mesh, vertices, trust_region_inflation_radius,
        min_distance);
}

namespace {
    template <typename T> inline bool approx(T a, T b, T tol = 1e-9)
    {
        return std::abs(a - b) <= tol;
    }

    /// @brief Compute the trust-region truncation ratio for a single vertex.
    ///
    /// Returns the largest β ∈ (0, 1] such that ‖xi + β·dxi − ci‖ ≤ ri.
    /// If xi + dxi is already inside the trust region, returns 1.0.
    ///
    /// @param xi  Current position of the vertex.
    /// @param dxi Proposed displacement of the vertex.
    /// @param ci  Center of the trust region.
    /// @param ri  Radius of the trust region.
    /// @return Truncation ratio β in (0, 1].
    inline double compute_trust_region_beta(
        const VectorMax3d& xi,
        const VectorMax3d& dxi,
        const VectorMax3d& ci,
        double ri)
    {
        if ((xi + dxi - ci).norm() <= ri) {
            return 1.0; // Already inside, no truncation needed
        }

        const double a = dxi.squaredNorm();
        if (a == 0) {
            return 1.0; // Zero displacement
        }

        // Solve ‖x + βΔx - c‖² = r² for β:
        // (x + βΔx - c)ᵀ(x + βΔx - c) - r²
        //  = ‖Δx‖² β² + 2(Δxᵀ(x - c)) β + ‖x - c‖² - r²
        //    { a }      {     b     }     {     c     }
        const double b = 2 * dxi.dot(xi - ci);
        const double c = (xi - ci).squaredNorm() - ri * ri;
        assert(c <= 0);                 // Inside the trust region, so c <= 0
        assert(b * b - 4 * a * c >= 0); // Roots must be real

        const double sign_b = (b >= 0) ? 1.0 : -1.0;
        const double sqrt_d = sqrt(b * b - 4 * a * c);

        // β = (-b + √(b² - 4ac)) / 2a, using a numerically stable form
        double beta;
        if (sign_b > 0) {
            beta = -2 * c / (b + sqrt_d);
        } else {
            beta = (-b + sqrt_d) / (2 * a);
        }
        // β should be the positive root
        assert(beta > 0);

        return beta;
    }
} // namespace

void TrustRegion::filter_step(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> x,
    Eigen::Ref<Eigen::MatrixXd> dx)
{
    assert(x.rows() == dx.rows() && x.cols() == dx.cols());

    int num_updates = 0;
    for (int i = 0; i < x.rows(); i++) {
        const VectorMax3d ci = trust_region_centers.row(i);
        const VectorMax3d xi = x.row(i);
        const VectorMax3d dxi = dx.row(i);

        // Check that the current position is within the trust region.
        assert((xi - ci).norm() <= trust_region_radii(i) + 1e-9);

        const double beta =
            compute_trust_region_beta(xi, dxi, ci, trust_region_radii(i));
        if (beta < 1.0) {
            dx.row(i).array() *= beta;
            assert(approx(
                (xi + dx.row(i).transpose() - ci).norm(),
                trust_region_radii(i)));

            num_updates++;
        }
    }

    // If a critical mass of vertices are restricted by the trust region,
    // we update the trust region to avoid excessive filtering.
    if (num_updates > update_threshold * mesh.num_vertices()) {
        logger().trace(
            "{:.1f}% of vertices restricted by trust region. Updating trust region.",
            100.0 * num_updates / mesh.num_vertices());
        should_update_trust_region = true;
    }
}

namespace {

    /// @brief Compute the closest points between two edges given their distance type.
    /// @return Pair (c_e, c_e') of closest points on ea and eb respectively.
    std::pair<VectorMax3d, VectorMax3d> edge_edge_closest_points(
        const VectorMax3d& ea0,
        const VectorMax3d& ea1,
        const VectorMax3d& eb0,
        const VectorMax3d& eb1,
        const ipc::EdgeEdgeDistanceType dtype)
    {
        switch (dtype) {
        case ipc::EdgeEdgeDistanceType::EA0_EB0:
            return { ea0, eb0 };
        case ipc::EdgeEdgeDistanceType::EA0_EB1:
            return { ea0, eb1 };
        case ipc::EdgeEdgeDistanceType::EA1_EB0:
            return { ea1, eb0 };
        case ipc::EdgeEdgeDistanceType::EA1_EB1:
            return { ea1, eb1 };
        case ipc::EdgeEdgeDistanceType::EA_EB0: {
            const double alpha = ipc::point_edge_closest_point(eb0, ea0, ea1);
            return { alpha * (ea1 - ea0) + ea0, eb0 };
        }
        case ipc::EdgeEdgeDistanceType::EA_EB1: {
            const double alpha = ipc::point_edge_closest_point(eb1, ea0, ea1);
            return { alpha * (ea1 - ea0) + ea0, eb1 };
        }
        case ipc::EdgeEdgeDistanceType::EA0_EB: {
            const double beta = ipc::point_edge_closest_point(ea0, eb0, eb1);
            return { ea0, beta * (eb1 - eb0) + eb0 };
        }
        case ipc::EdgeEdgeDistanceType::EA1_EB: {
            const double beta = ipc::point_edge_closest_point(ea1, eb0, eb1);
            return { ea1, beta * (eb1 - eb0) + eb0 };
        }
        case ipc::EdgeEdgeDistanceType::EA_EB: {
            const Eigen::Vector2d params =
                ipc::edge_edge_closest_point(ea0, ea1, eb0, eb1);
            return {
                params[0] * (ea1 - ea0) + ea0,
                params[1] * (eb1 - eb0) + eb0,
            };
        }
        default:
            assert(false && "unexpected EdgeEdgeDistanceType");
            return { ea0, eb0 }; // unreachable
        }
    }
} // anonymous namespace

void TrustRegion::planar_filter_step(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> x,
    Eigen::Ref<Eigen::MatrixXd> dx,
    const NormalCollisions& collisions)
{
    assert(x.rows() == dx.rows() && x.cols() == dx.cols());
    assert(relaxed_radius_scaling > 0 && relaxed_radius_scaling < 1);

    // Collision mesh vertex positions (may differ from x if hidden DOFs exist)
    const bool has_hidden_dofs = x.rows() != mesh.num_vertices();
    Eigen::MatrixXd vertices;
    if (has_hidden_dofs) {
        vertices = mesh.vertices(x);
    } else {
        vertices = x;
    }

    // Map collision-mesh vertex index → full-mesh row index in x/dx
    const Eigen::VectorXi& full_id_map = mesh.to_full_vertex_id();
    auto to_full = [&](const index_t coll_id) -> int {
        return has_hidden_dofs ? full_id_map(coll_id)
                               : static_cast<int>(coll_id);
    };

    // Per-vertex truncation ratios (1 = no truncation)
    Eigen::VectorXd t_v = Eigen::VectorXd::Ones(x.rows());

    // ---- Face-Vertex Collisions ----------------------------------------
    // FaceVertexNormalCollision always has known_dtype() == P_T (interior),
    // guaranteed by the OGC feasible region filter during build().
    for (const auto& fv : collisions.fv_collisions) {
        const auto ids = fv.vertex_ids(mesh.edges(), mesh.faces());
        // ids = [vi, f0i, f1i, f2i] — collision mesh indices

        const Eigen::Vector3d x_v = vertices.row(ids[0]);
        const Eigen::Vector3d x_a = vertices.row(ids[1]);
        const Eigen::Vector3d x_b = vertices.row(ids[2]);
        const Eigen::Vector3d x_c = vertices.row(ids[3]);

        // Closest point on triangle (t0=x_a, t1=x_b, t2=x_c) to vertex x_v.
        // Returns barycentric (u, v); reconstruct as (1-u-v)*x_a + u*x_b +
        // v*x_c.
        const Eigen::Vector2d uv =
            point_triangle_closest_point(x_v, x_a, x_b, x_c);
        const Eigen::Vector3d c_vt =
            x_a + uv[0] * (x_b - x_a) + uv[1] * (x_c - x_a);

        // Normal pointing from triangle toward vertex
        const Eigen::Vector3d diff = x_v - c_vt;
        const double dist = diff.norm();
        if (dist < 1e-10) {
            continue; // Degenerate (primitives coincident) — skip pair
        }
        const Eigen::Vector3d n = diff / dist;

        // Full-mesh IDs for indexing into dx
        const int fid_v = to_full(ids[0]);
        const int fid_a = to_full(ids[1]);
        const int fid_b = to_full(ids[2]);
        const int fid_c = to_full(ids[3]);

        const Eigen::Vector3d dx_v = dx.row(fid_v);
        const Eigen::Vector3d dx_a = dx.row(fid_a);
        const Eigen::Vector3d dx_b = dx.row(fid_b);
        const Eigen::Vector3d dx_c = dx.row(fid_c);

        // Approach speeds (Eq. 9–10 in DAT paper)
        // δ_v: how fast vertex is approaching triangle (negative n-component)
        const double delta_v = std::max(-dx_v.dot(n), 0.0);
        // δ_t: max approach speed of any triangle vertex toward vertex (Eq. 10)
        const double delta_t =
            std::max({ dx_a.dot(n), dx_b.dot(n), dx_c.dot(n), 0.0 });

        // Adaptive λ: allocate gap proportional to approach speeds (Eq. 11)
        const double lambda = (delta_v == 0.0 && delta_t == 0.0)
            ? 0.5
            : delta_t / (delta_t + delta_v);

        // Division plane point: interpolated between c_vt and x_v
        const Eigen::Vector3d p = c_vt + lambda * diff;

        // Truncate all four vertices using this division plane
        t_v[fid_v] =
            std::min(t_v[fid_v], planar_truncation_ratio(x_v, dx_v, n, p));
        t_v[fid_a] =
            std::min(t_v[fid_a], planar_truncation_ratio(x_a, dx_a, n, p));
        t_v[fid_b] =
            std::min(t_v[fid_b], planar_truncation_ratio(x_b, dx_b, n, p));
        t_v[fid_c] =
            std::min(t_v[fid_c], planar_truncation_ratio(x_c, dx_c, n, p));
    }

    // ---- Edge-Edge Collisions ------------------------------------------
    for (const auto& ee : collisions.ee_collisions) {
        const auto ids = ee.vertex_ids(mesh.edges(), mesh.faces());
        // ids = [ea0i, ea1i, eb0i, eb1i] — collision mesh indices

        const Eigen::Vector3d ea0 = vertices.row(ids[0]);
        const Eigen::Vector3d ea1 = vertices.row(ids[1]);
        const Eigen::Vector3d eb0 = vertices.row(ids[2]);
        const Eigen::Vector3d eb1 = vertices.row(ids[3]);

        // Closest points on each edge (dispatch on stored distance type)
        const auto [c_e, c_ep] =
            edge_edge_closest_points(ea0, ea1, eb0, eb1, ee.known_dtype());

        // Normal pointing from ea toward eb
        const Eigen::Vector3d diff = c_ep - c_e;
        const double dist = diff.norm();
        if (dist < 1e-10) {
            continue; // Degenerate — skip pair
        }
        const Eigen::Vector3d n = diff / dist;

        const int fid_ea0 = to_full(ids[0]);
        const int fid_ea1 = to_full(ids[1]);
        const int fid_eb0 = to_full(ids[2]);
        const int fid_eb1 = to_full(ids[3]);

        const Eigen::Vector3d dx_ea0 = dx.row(fid_ea0);
        const Eigen::Vector3d dx_ea1 = dx.row(fid_ea1);
        const Eigen::Vector3d dx_eb0 = dx.row(fid_eb0);
        const Eigen::Vector3d dx_eb1 = dx.row(fid_eb1);

        // δ_e: max approach speed of ea toward eb in +n direction (Eq. 12)
        const double delta_e = std::max({ dx_ea0.dot(n), dx_ea1.dot(n), 0.0 });
        // δ_e': max approach speed of eb toward ea in −n direction (Eq. 12)
        const double delta_ep =
            std::max({ -dx_eb0.dot(n), -dx_eb1.dot(n), 0.0 });

        // Adaptive λ (Eq. 12 in DAT paper)
        const double lambda = (delta_e == 0.0 && delta_ep == 0.0)
            ? 0.5
            : delta_ep / (delta_e + delta_ep);

        // Division plane point: interpolated between c_e and c_ep
        const Eigen::Vector3d p = c_e + lambda * diff;

        // Truncate all four edge vertices
        t_v[fid_ea0] =
            std::min(t_v[fid_ea0], planar_truncation_ratio(ea0, dx_ea0, n, p));
        t_v[fid_ea1] =
            std::min(t_v[fid_ea1], planar_truncation_ratio(ea1, dx_ea1, n, p));
        t_v[fid_eb0] =
            std::min(t_v[fid_eb0], planar_truncation_ratio(eb0, dx_eb0, n, p));
        t_v[fid_eb1] =
            std::min(t_v[fid_eb1], planar_truncation_ratio(eb1, dx_eb1, n, p));
    }

    // ---- Edge-Vertex Collisions ----------------------------------------
    // ids = [vi, e0i, e1i, -1] — collision mesh indices
    for (const auto& ev : collisions.ev_collisions) {
        const auto ids = ev.vertex_ids(mesh.edges(), mesh.faces());

        const VectorMax3d x_v = vertices.row(ids[0]);
        const VectorMax3d x_e0 = vertices.row(ids[1]);
        const VectorMax3d x_e1 = vertices.row(ids[2]);

        const double alpha = point_edge_closest_point(x_v, x_e0, x_e1);
        assert(-1e-8 <= alpha && alpha <= (1 + 1e-8));
        const VectorMax3d c_e = x_e0 + alpha * (x_e1 - x_e0);

        // Normal pointing from edge toward vertex
        const VectorMax3d diff = x_v - c_e;
        const double dist = diff.norm();
        if (dist < 1e-10) {
            continue; // Degenerate — skip pair
        }
        const VectorMax3d n = diff / dist;

        const int fid_v = to_full(ids[0]);
        const int fid_e0 = to_full(ids[1]);
        const int fid_e1 = to_full(ids[2]);

        const VectorMax3d dx_v = dx.row(fid_v);
        const VectorMax3d dx_e0 = dx.row(fid_e0);
        const VectorMax3d dx_e1 = dx.row(fid_e1);

        // δ_v: vertex approaching edge (moving in −n direction)
        const double delta_v = std::max(-dx_v.dot(n), 0.0);
        // δ_e: max approach speed of any edge vertex toward vertex (in +n)
        const double delta_e = std::max({ dx_e0.dot(n), dx_e1.dot(n), 0.0 });

        const double lambda = (delta_v == 0.0 && delta_e == 0.0)
            ? 0.5
            : delta_e / (delta_e + delta_v);

        const VectorMax3d p = c_e + lambda * diff;

        t_v[fid_v] =
            std::min(t_v[fid_v], planar_truncation_ratio(x_v, dx_v, n, p));
        t_v[fid_e0] =
            std::min(t_v[fid_e0], planar_truncation_ratio(x_e0, dx_e0, n, p));
        t_v[fid_e1] =
            std::min(t_v[fid_e1], planar_truncation_ratio(x_e1, dx_e1, n, p));
    }

    // ---- Vertex-Vertex Collisions --------------------------------------
    // ids = [v0i, v1i, -1, -1] — collision mesh indices
    for (const auto& vv : collisions.vv_collisions) {
        const auto ids = vv.vertex_ids(mesh.edges(), mesh.faces());

        const VectorMax3d x_a = vertices.row(ids[0]);
        const VectorMax3d x_b = vertices.row(ids[1]);

        // Normal pointing from v0 toward v1
        const VectorMax3d diff = x_b - x_a;
        const double dist = diff.norm();
        if (dist < 1e-10) {
            continue; // Degenerate — skip pair
        }
        const VectorMax3d n = diff / dist;

        const int fid_a = to_full(ids[0]);
        const int fid_b = to_full(ids[1]);

        const VectorMax3d dx_a = dx.row(fid_a);
        const VectorMax3d dx_b = dx.row(fid_b);

        // δ_a: v0 approaching v1 (moving in +n); δ_b: v1 approaching v0 (-n)
        const double delta_a = std::max(dx_a.dot(n), 0.0);
        const double delta_b = std::max(-dx_b.dot(n), 0.0);

        const double lambda = (delta_a == 0.0 && delta_b == 0.0)
            ? 0.5
            : delta_b / (delta_a + delta_b);

        const VectorMax3d p = x_a + lambda * diff;

        t_v[fid_a] =
            std::min(t_v[fid_a], planar_truncation_ratio(x_a, dx_a, n, p));
        t_v[fid_b] =
            std::min(t_v[fid_b], planar_truncation_ratio(x_b, dx_b, n, p));
    }

    // ---- Isotropic Fallback --------------------------------------------
    // For primitives beyond the query radius, fall back to isotropic
    // trust-region clamping (same as filter_step) to prevent penetration.
    for (int i = 0; i < x.rows(); i++) {
        const VectorMax3d xi = x.row(i);
        const VectorMax3d dxi = dx.row(i);
        const VectorMax3d ci = trust_region_centers.row(i);
        // Use trust_region_inflation_radius instead of trust_region_radii(v).
        // This is the key difference between Planar-DAT and Isotropic-DAT.
        // This allows us to preserve more motion for vertices that are near the
        // boundary of the trust region, rather than restricting them to a
        // smaller trust region radius.
        const double beta = compute_trust_region_beta(
            xi, dxi, ci, trust_region_inflation_radius);
        t_v[i] = std::min(t_v[i], beta);
    }

    // ---- Apply Truncations ---------------------------------------------
    int num_updates = 0;
    for (int v = 0; v < x.rows(); v++) {
        if (t_v[v] < 1.0) {
            dx.row(v) *= t_v[v];
            num_updates++;
        }
    }

    // Mirror filter_step: if a critical mass of vertices are restricted,
    // flag the trust region for re-centering on the next iteration.
    if (num_updates > update_threshold * mesh.num_vertices()) {
        logger().trace(
            "{:.1f}% of vertices restricted by planar trust region. Updating trust region.",
            100.0 * num_updates / mesh.num_vertices());
        should_update_trust_region = true;
    }
}

double TrustRegion::planar_truncation_ratio(
    Eigen::ConstRef<VectorMax3d> x_u,
    Eigen::ConstRef<VectorMax3d> dx_u,
    Eigen::ConstRef<VectorMax3d> n,
    Eigen::ConstRef<VectorMax3d> p) const
{
    const double denom = dx_u.dot(n);
    if (std::abs(denom) < 1e-15) {
        return 1.0; // dx parallel to plane → no crossing
    }
    // Ray-plane intersection: x_u + t_i * dx_u is on plane when
    // (x_u + t_i * dx_u - p) · n = 0  →  t_i = -(x_u - p) · n / denom
    const double t_i = -(x_u - p).dot(n) / denom;
    // Only constrain if the crossing is meaningfully in the future.
    // Edge/triangle vertices lie geometrically on the division plane (provable
    // from the closest-point construction), so t_i is O(ε_machine) rather
    // than exactly 0. Treating any t_i below this threshold as "no crossing"
    // avoids zeroing out their displacements when they are retreating.
    constexpr double T_EPS = 1e-10;
    if (t_i < T_EPS || t_i >= 1.0 / relaxed_radius_scaling) {
        return 1.0;
    }
    return relaxed_radius_scaling * t_i;
}

} // namespace ipc::ogc