#include "trust_region.hpp"

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

    candidates.build(
        mesh, vertices, trust_region_inflation_radius + min_distance,
        broad_phase);

    // Compute the trust region distances for the reduced set of collision
    // vertices.
    Eigen::VectorXd reduced_trust_region_radii =
        candidates.compute_per_vertex_safe_distances(
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
        candidates, mesh, vertices, trust_region_inflation_radius,
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

void TrustRegion::planar_filter_step(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> x,
    Eigen::Ref<Eigen::MatrixXd> dx)
{
    assert(x.rows() == dx.rows() && x.cols() == dx.cols());
    assert(relaxed_radius_scaling > 0 && relaxed_radius_scaling < 1);

    const int d = x.cols();

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

    for (size_t i = 0; i < candidates.size(); ++i) {
        const CollisionStencil& c = candidates[i];
        const auto ids = c.vertex_ids(mesh.edges(), mesh.faces());
        const int nv = c.num_vertices();
        const VectorMax12d pos = c.dof(vertices, mesh.edges(), mesh.faces());

        // dv = c_first − c_second; positive coeffs → first primitive
        VectorMax4d coeffs;
        const VectorMax3d dv = c.compute_distance_vector(pos, coeffs);
        const double dist = dv.norm();
        if (dist < 1e-10) {
            continue;
        }

        const VectorMax3d n = dv / dist;

        VectorMax3d c_first = VectorMax3d::Zero(d);
        VectorMax3d c_second = VectorMax3d::Zero(d);
        double delta_first = 0, delta_second = 0;

        for (int j = 0; j < nv; ++j) {
            const int fid = to_full(ids[j]);
            const VectorMax3d xj = pos.segment(j * d, d);
            const VectorMax3d dxj = dx.row(fid);
            if (coeffs[j] > 0) {
                c_first += coeffs[j] * xj;
                delta_first = std::max(delta_first, -dxj.dot(n));
            } else if (coeffs[j] < 0) {
                c_second -= coeffs[j] * xj;
                delta_second = std::max(delta_second, dxj.dot(n));
            }
        }
        delta_first = std::max(delta_first, 0.0);
        delta_second = std::max(delta_second, 0.0);

        const double lambda = (delta_first == 0 && delta_second == 0)
            ? 0.5
            : delta_second / (delta_first + delta_second);

        // Division plane: p = c_second + λ·dv  (λ=0 → plane at second prim,
        // λ=1 → plane at first prim).
        const VectorMax3d p = c_second + lambda * dv;

        for (int j = 0; j < nv; ++j) {
            assert(ids[j] >= 0); // All candidates should be valid vertices
            const int fid = to_full(ids[j]);
            const VectorMax3d xj = pos.segment(j * d, d);
            const VectorMax3d dxj = dx.row(fid);
            t_v[fid] =
                std::min(t_v[fid], planar_truncation_ratio(xj, dxj, n, p));
        }
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