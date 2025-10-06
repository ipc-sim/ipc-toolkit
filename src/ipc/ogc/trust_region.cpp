#include "trust_region.hpp"

namespace ipc::ogc {

void TrustRegion::warm_start_time_step(
    const CollisionMesh& mesh,
    Eigen::Ref<Eigen::MatrixXd> x,
    Eigen::ConstRef<Eigen::MatrixXd> pred_x,
    NormalCollisions& collisions,
    const double dhat,
    const double min_distance,
    const std::shared_ptr<BroadPhase>& broad_phase)
{
    const int N = x.rows();
    assert(x.rows() == pred_x.rows() && x.cols() == pred_x.cols());
    assert(trust_region_radii.size() == N);

    // Compute the norm of the translation for each vertex.
    Eigen::VectorXd dx_norm = (pred_x - x).rowwise().norm();
    assert(dx_norm.size() == N);

    // Compute a trust region inflation radius based on the distance
    // between the current positions x and the predicted positions x̂.
    trust_region_inflation_radius = std::max(2 * dhat, dx_norm.maxCoeff());

    // Initialize the trust region around x
    update(mesh, x, collisions, min_distance, broad_phase);

    int num_updates = 0;
    for (int i = 0; i < N; i++) {
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
    const std::shared_ptr<BroadPhase>& broad_phase)
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

bool TrustRegion::filter_step(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> x,
    Eigen::Ref<Eigen::MatrixXd> dx)
{
    assert(x.rows() == dx.rows() && x.cols() == dx.cols());
    const int N = x.size() / 3;

    int num_updates = 0;
    for (int i = 0; i < N; i++) {
        const VectorMax3d ci = trust_region_centers.row(i);
        const VectorMax3d xi = x.row(i);
        const VectorMax3d dxi = dx.row(i);

        // Check that the current position is within the trust region.
        assert(approx((xi - ci).norm(), trust_region_radii(i)));

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
            assert(approx((xi + dx.row(i) - ci).norm(), trust_region_radii(i)));

            num_updates++;
        }
    }

    // If a critical mass of vertices are restricted by the trust region,
    // we update the trust region to avoid excessive filtering.
    if (num_updates > update_threshold * mesh.num_vertices()) {
        logger().trace(
            "{:.1f}% of vertices restricted by trust region. Updating trust region.",
            100.0 * num_updates / mesh.num_vertices());
        // init_trust_region(x);
        return true;
    }
    return false;
}

} // namespace ipc::ogc