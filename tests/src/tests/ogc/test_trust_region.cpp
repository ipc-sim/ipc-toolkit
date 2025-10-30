#include <catch2/catch_all.hpp>

#include <ipc/ogc/trust_region.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>

#include <Eigen/Dense>

using namespace ipc;

TEST_CASE("TrustRegion warm_start_time_step sets inflation radius", "[ogc]")
{
    constexpr double dhat = 0.02; // offset distance

    // Create a simple 3-vertex mesh (straight line)
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, /**/ 1, 0, 0, /**/ 2, 0, 0, /**/ 3, 0, 0;

    CollisionMesh mesh(V);

    ipc::ogc::TrustRegion tr(dhat);

    Eigen::MatrixXd x = V;
    // Predicted positions: small motion of 0.01 for all vertices
    Eigen::MatrixXd pred_x = x;
    pred_x.col(0).array() += dhat / 2.0;

    ipc::NormalCollisions collisions;

    tr.warm_start_time_step(mesh, x, pred_x, collisions, dhat);

    CHECK(tr.trust_region_inflation_radius == Catch::Approx(2 * dhat));

    // After warm start, the trust region should be initialized and not
    // scheduled for update on next call.
    CHECK(tr.should_update_trust_region == false);

    x = V;
    pred_x = x;
    pred_x.col(0).array() += 2.0;

    tr.warm_start_time_step(mesh, x, pred_x, collisions, dhat);

    CHECK(tr.trust_region_inflation_radius == Catch::Approx(2.0));
    CHECK(tr.should_update_trust_region == false);

    CHECK(x(0, 0) == Catch::Approx(0.45));
    CHECK(x(1, 0) == Catch::Approx(1.45));
    CHECK(x(2, 0) == Catch::Approx(2.45));
    CHECK(x(3, 0) == Catch::Approx(3.45));

    x = V;
    pred_x = x.array() + 2.0;

    tr.warm_start_time_step(mesh, x, pred_x, collisions, dhat);

    CHECK(tr.trust_region_inflation_radius == Catch::Approx(sqrt(12.0)));
    CHECK(tr.should_update_trust_region == false);

    CHECK(x(0, 0) == Catch::Approx(V(0, 0) + 0.9 / sqrt(12.0)));
    CHECK(x(1, 0) == Catch::Approx(V(1, 0) + 0.9 / sqrt(12.0)));
    CHECK(x(2, 0) == Catch::Approx(V(2, 0) + 0.9 / sqrt(12.0)));
    CHECK(x(3, 0) == Catch::Approx(V(3, 0) + 0.9 / sqrt(12.0)));
}

TEST_CASE("TrustRegion update_if_needed toggles update flag", "[ogc]")
{
    Eigen::MatrixXd V(2, 3);
    V << 0, 0, 0, /**/ 1, 0, 0;

    Eigen::MatrixXi E;
    Eigen::MatrixXi F;

    CollisionMesh mesh(V, E, F);

    ipc::ogc::TrustRegion tr(0.02);

    Eigen::MatrixXd x = V;

    ipc::NormalCollisions collisions;

    // Ensure the flag is honored: set true and call update_if_needed
    tr.should_update_trust_region = true;
    tr.update_if_needed(mesh, x, collisions);
    CHECK(tr.should_update_trust_region == false);

    // Calling again should not call update (flag remains false)
    tr.update_if_needed(mesh, x, collisions);
    CHECK(tr.should_update_trust_region == false);
}

TEST_CASE(
    "TrustRegion filter_step scales step to boundary and requests update",
    "[ogc]")
{
    // Single vertex example
    Eigen::MatrixXd V(1, 3);
    V << 0.5, 0.0, 0.0;

    CollisionMesh mesh(V);

    constexpr double r = 1.0;
    ipc::ogc::TrustRegion tr(r);

    // Set trust region center at origin and radius r > ||xi - ci|| so xi is
    // inside
    tr.trust_region_centers = Eigen::MatrixXd::Zero(1, 3);
    tr.trust_region_radii = Eigen::VectorXd::Constant(V.rows(), r);

    Eigen::MatrixXd x = V; // current positions

    // Proposed displacement: outward along x so that alpha > r
    Eigen::MatrixXd dx(1, 3);
    dx << 1.0, 0.0, 0.0; // adds to x gives 1.5 > r

    ipc::NormalCollisions collisions;

    // Initially not scheduled for update
    tr.should_update_trust_region = false;

    // Call filter_step
    tr.filter_step(mesh, x, dx);

    // Compute expected beta using same stable formula as implementation
    const Eigen::RowVector3d ci = tr.trust_region_centers.row(0);
    const Eigen::RowVector3d xi = x.row(0);
    const Eigen::RowVector3d dxi = Eigen::RowVector3d(1.0, 0.0, 0.0);

    // After filtering, dx should be scaled by expected_beta
    CHECK(dx(0, 0) == Catch::Approx(0.5));
    CHECK(dx(0, 1) == Catch::Approx(0.0));
    CHECK(dx(0, 2) == Catch::Approx(0.0));

    // Since we restricted one vertex, should_update_trust_region should be set
    CHECK(tr.should_update_trust_region == true);
}
