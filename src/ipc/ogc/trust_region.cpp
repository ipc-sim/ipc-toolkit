#include "trust_region.hpp"

#include <ipc/tangent/closest_point.hpp>
#include <ipc/distance/distance_type.hpp>

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

        const double alpha = (xi + dxi - ci).norm();
        if (alpha > trust_region_radii(i)) {
            // Solve ‖x + βΔx - c‖² = r² for β:
            // (x + βΔx - c)ᵀ(x + βΔx - c) - r²
            //  = ‖Δx‖² β² + 2(Δxᵀ(x - c)) β + ‖x - c‖² - r²
            //    { a }      {     b     }     {     c     }

            const double a = dxi.squaredNorm();
            assert(a > 0); // Δx should not be zero

            // b can be positive or negative
            const double b = 2 * dxi.dot(xi - ci);

            const double c = (xi - ci).squaredNorm()
                - trust_region_radii(i) * trust_region_radii(i);
            assert(c <= 0); // Inside the trust region, so c <= 0

            // Roots must be real
            assert(b * b - 4 * a * c >= 0);

            const double sign_b = (b >= 0) ? 1.0 : -1.0;
            const double sqrt_d = sqrt(b * b - 4 * a * c);

            // β = (-b + √(b² - 4ac)) / 2a, but we use a more stable form below
            double beta;
            if (sign_b > 0) {
                beta = -2 * c / (b + sqrt_d);
            } else {
                beta = (-b + sqrt_d) / (2 * a);
            }
            // β should be the positive root
            assert(beta > 0);

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
    /// @brief Compute the truncation ratio for a vertex given a division plane.
    ///
    /// Finds the parameter t_i such that x_u + t_i * dx_u lies on the plane
    /// (n, p). If the crossing is in the future and within one step, returns
    /// relaxation_ratio * t_i; otherwise returns 1 (no truncation).
    ///
    /// @param x_u Current position of vertex u.
    /// @param dx_u Proposed displacement of vertex u.
    /// @param n Division plane normal (unit vector).
    /// @param p A point on the division plane.
    /// @param relaxation_ratio Safety margin γ_r ∈ (0,1).
    /// @return Truncation ratio in (0, 1].
    double planar_truncation_ratio(
        const Eigen::Vector3d& x_u,
        const Eigen::Vector3d& dx_u,
        const Eigen::Vector3d& n,
        const Eigen::Vector3d& p,
        const double relaxation_ratio)
    {
        const double denom = dx_u.dot(n);
        if (std::abs(denom) < 1e-15) {
            return 1.0; // dx parallel to plane → no crossing
        }
        // Ray-plane intersection: x_u + t_i * dx_u is on plane when
        // (x_u + t_i * dx_u - p) · n = 0  →  t_i = -(x_u - p) · n / denom
        const double t_i = -(x_u - p).dot(n) / denom;
        // Only constrain if crossing is in the future and within one step
        if (t_i < 0.0 || t_i >= 1.0 / relaxation_ratio) {
            return 1.0;
        }
        return relaxation_ratio * t_i;
    }

    /// @brief Compute the closest points between two edges given their distance type.
    /// @return Pair (c_e, c_e') of closest points on ea and eb respectively.
    std::pair<Eigen::Vector3d, Eigen::Vector3d> edge_edge_closest_points(
        const Eigen::Vector3d& ea0,
        const Eigen::Vector3d& ea1,
        const Eigen::Vector3d& eb0,
        const Eigen::Vector3d& eb1,
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
            return { (1 - alpha) * ea0 + alpha * ea1, eb0 };
        }
        case ipc::EdgeEdgeDistanceType::EA_EB1: {
            const double alpha = ipc::point_edge_closest_point(eb1, ea0, ea1);
            return { (1 - alpha) * ea0 + alpha * ea1, eb1 };
        }
        case ipc::EdgeEdgeDistanceType::EA0_EB: {
            const double beta = ipc::point_edge_closest_point(ea0, eb0, eb1);
            return { ea0, (1 - beta) * eb0 + beta * eb1 };
        }
        case ipc::EdgeEdgeDistanceType::EA1_EB: {
            const double beta = ipc::point_edge_closest_point(ea1, eb0, eb1);
            return { ea1, (1 - beta) * eb0 + beta * eb1 };
        }
        default: { // EA_EB (interior) or AUTO
            const Eigen::Vector2d params =
                ipc::edge_edge_closest_point(ea0, ea1, eb0, eb1);
            return {
                (1 - params[0]) * ea0 + params[0] * ea1,
                (1 - params[1]) * eb0 + params[1] * eb1,
            };
        }
        }
    }
} // anonymous namespace

void TrustRegion::planar_filter_step(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> x,
    Eigen::Ref<Eigen::MatrixXd> dx,
    const NormalCollisions& collisions,
    const double query_radius,
    const double relaxation_ratio) const
{
    assert(x.rows() == dx.rows() && x.cols() == dx.cols());
    assert(relaxation_ratio > 0 && relaxation_ratio < 1);

    // Collision mesh vertex positions (may differ from x if hidden DOFs exist)
    const bool has_hidden_dofs = x.rows() != mesh.num_vertices();
    const Eigen::MatrixXd vertices =
        has_hidden_dofs ? mesh.vertices(x) : x;

    // Map collision-mesh vertex index → full-mesh row index in x/dx
    const Eigen::VectorXi& full_id_map = mesh.to_full_vertex_id();
    auto to_full = [&](const index_t coll_id) -> int {
        return has_hidden_dofs ? full_id_map(coll_id) : (int)coll_id;
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
        // Returns barycentric (u, v); reconstruct as (1-u-v)*x_a + u*x_b + v*x_c.
        const Eigen::Vector2d uv =
            point_triangle_closest_point(x_v, x_a, x_b, x_c);
        const Eigen::Vector3d c_vt =
            (1 - uv[0] - uv[1]) * x_a + uv[0] * x_b + uv[1] * x_c;

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
        // δ_t: how fast any triangle vertex is approaching vertex (positive n-component)
        const double delta_v = std::max(-dx_v.dot(n), 0.0);
        const double delta_t =
            std::max({ dx_a.dot(n), dx_b.dot(n), dx_c.dot(n), 0.0 });

        // Adaptive λ: allocate gap proportional to approach speeds (Eq. 11)
        const double lambda = (delta_v == 0.0 && delta_t == 0.0)
            ? 0.5
            : delta_t / (delta_t + delta_v);

        // Division plane point: interpolated between c_vt and x_v
        const Eigen::Vector3d p = c_vt + lambda * diff;

        // Truncate all four vertices using this division plane
        t_v[fid_v] = std::min(t_v[fid_v],
            planar_truncation_ratio(x_v, dx_v, n, p, relaxation_ratio));
        t_v[fid_a] = std::min(t_v[fid_a],
            planar_truncation_ratio(x_a, dx_a, n, p, relaxation_ratio));
        t_v[fid_b] = std::min(t_v[fid_b],
            planar_truncation_ratio(x_b, dx_b, n, p, relaxation_ratio));
        t_v[fid_c] = std::min(t_v[fid_c],
            planar_truncation_ratio(x_c, dx_c, n, p, relaxation_ratio));
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

        // δ_e: max approach of ea toward eb (in +n direction)
        // δ_e': max approach of eb toward ea (in −n direction)
        const double delta_e =
            std::max({ dx_ea0.dot(n), dx_ea1.dot(n), 0.0 });
        const double delta_ep =
            std::max({ -dx_eb0.dot(n), -dx_eb1.dot(n), 0.0 });

        // Adaptive λ (Eq. 12 in DAT paper)
        const double lambda = (delta_e == 0.0 && delta_ep == 0.0)
            ? 0.5
            : delta_ep / (delta_e + delta_ep);

        // Division plane point: interpolated between c_e and c_ep
        const Eigen::Vector3d p = c_e + lambda * diff;

        // Truncate all four edge vertices
        t_v[fid_ea0] = std::min(t_v[fid_ea0],
            planar_truncation_ratio(ea0, dx_ea0, n, p, relaxation_ratio));
        t_v[fid_ea1] = std::min(t_v[fid_ea1],
            planar_truncation_ratio(ea1, dx_ea1, n, p, relaxation_ratio));
        t_v[fid_eb0] = std::min(t_v[fid_eb0],
            planar_truncation_ratio(eb0, dx_eb0, n, p, relaxation_ratio));
        t_v[fid_eb1] = std::min(t_v[fid_eb1],
            planar_truncation_ratio(eb1, dx_eb1, n, p, relaxation_ratio));
    }

    // ---- Isotropic Fallback --------------------------------------------
    // For primitives beyond the query radius, apply an isotropic cap so
    // that no vertex moves more than 0.5 * γ_r * r_q (Sec. 4 in DAT paper).
    const double iso_cap = 0.5 * relaxation_ratio * query_radius;
    for (int v = 0; v < (int)x.rows(); v++) {
        const double dx_norm = dx.row(v).norm();
        if (dx_norm > 0 && t_v[v] * dx_norm > iso_cap) {
            t_v[v] = std::min(t_v[v], iso_cap / dx_norm);
        }
    }

    // ---- Apply Truncations ---------------------------------------------
    for (int v = 0; v < (int)x.rows(); v++) {
        if (t_v[v] < 1.0) {
            dx.row(v) *= t_v[v];
        }
    }
}

} // namespace ipc::ogc