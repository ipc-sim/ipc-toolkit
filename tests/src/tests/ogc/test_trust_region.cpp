#include <catch2/catch_all.hpp>

#include <ipc/ogc/trust_region.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>

#include <Eigen/Dense>

using namespace ipc;

// ============================================================================
// Helper: build a minimal mesh with one triangle and one floating vertex above
// its interior, close enough to generate an FV collision with the given dhat.
//
// Layout (y-up):
//   V[0] = (0,  0, 0)   \
//   V[1] = (2,  0, 0)    > triangle in y=0 plane
//   V[2] = (1,  0, 2)   /
//   V[3] = (1, gap, 1)    <- vertex directly above triangle centroid
//
// Returned mesh has one face [0,1,2] and three boundary edges.
static std::pair<CollisionMesh, Eigen::MatrixXd> make_fv_mesh(double gap)
{
    Eigen::MatrixXd V(4, 3);
    // clang-format off
    V << 0, 0,   0,
         2, 0,   0,
         1, 0,   2,
         1, gap, 1;  // vertex above centroid
    // clang-format on

    Eigen::MatrixXi E(3, 2);
    E << 0, 1, /**/ 1, 2, /**/ 0, 2;

    Eigen::MatrixXi F(1, 3);
    F << 0, 1, 2;

    return { CollisionMesh(V, E, F), V };
}

// Helper: build a mesh with two perpendicular edges separated by a gap.
//   Edge A: (0, 0, 0) → (1, 0, 0)     (along x)
//   Edge B: (0.5, gap, -0.5) → (0.5, gap, 0.5)  (along z, above A's midpoint)
static std::pair<CollisionMesh, Eigen::MatrixXd> make_ee_mesh(double gap)
{
    Eigen::MatrixXd V(4, 3);
    // clang-format off
    V << 0.0,   0,    0.0,   // ea0
         1.0,   0,    0.0,   // ea1
         0.5, gap,   -0.5,   // eb0
         0.5, gap,    0.5;   // eb1
    // clang-format on

    Eigen::MatrixXi E(2, 2);
    E << 0, 1, /**/ 2, 3;

    Eigen::MatrixXi F(0, 3); // no faces

    return { CollisionMesh(V, E, F), V };
}

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

// ============================================================================
// Planar-DAT (planar_filter_step) tests
// ============================================================================

TEST_CASE(
    "planar_filter_step: separating vertex is capped by trust region",
    "[ogc][planar_dat]")
{
    // Vertex moving AWAY from the triangle is not restricted by the planar
    // constraint, but the trust-region fallback (same as filter_step) still
    // applies to prevent penetration with undetected geometry beyond the
    // query radius.
    const double gap = 0.5;
    const double dhat = 0.6; // dhat > gap so the FV pair is active
    auto [mesh, V] = make_fv_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);
    // Vertex 3 moves upward (away from triangle) with a large step that
    // exceeds the trust region inflation radius (2 * dhat = 1.2).
    dx.row(3) << 0.0, 2.0, 0.0;

    tr.planar_filter_step(mesh, x, dx, collisions);

    // The planar constraint leaves the step alone (moving away), but the
    // isotropic trust-region fallback clamps the displacement to the
    // trust region inflation radius.
    CHECK(dx(3, 0) == Catch::Approx(0.0));
    CHECK(dx(3, 2) == Catch::Approx(0.0));
    // Displacement was reduced and is still in the correct direction
    CHECK(dx(3, 1) > 0.0);
    CHECK(dx(3, 1) < 2.0);
    // Final position must lie on the inflation-radius sphere boundary
    const double dist_from_center =
        (x.row(3) + dx.row(3) - tr.trust_region_centers.row(3)).norm();
    CHECK(
        dist_from_center
        == Catch::Approx(tr.trust_region_inflation_radius).epsilon(1e-6));
}

TEST_CASE(
    "planar_filter_step: tangential vertex is not truncated by planar DAT",
    "[ogc][planar_dat]")
{
    // Vertex sliding sideways (perpendicular to contact normal) should not
    // be restricted by Planar-DAT, but the trust-region fallback still
    // applies.
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_fv_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);
    // Vertex 3 moves sideways (contact normal is (0,1,0), dx is in xz-plane)
    // with a large step that exceeds the trust region inflation radius
    // (2 * dhat = 1.2).
    dx.row(3) << 2.0, 0.0, 0.0;

    tr.planar_filter_step(mesh, x, dx, collisions);

    // Tangential motion: Planar-DAT does not truncate, but the isotropic
    // trust-region fallback caps the displacement to the inflation radius.
    CHECK(dx(3, 1) == Catch::Approx(0.0));
    CHECK(dx(3, 2) == Catch::Approx(0.0));
    // Displacement was reduced and is still in the correct direction
    CHECK(dx(3, 0) > 0.0);
    CHECK(dx(3, 0) < 2.0);
    // Final position must lie on the inflation-radius sphere boundary
    const double dist_from_center =
        (x.row(3) + dx.row(3) - tr.trust_region_centers.row(3)).norm();
    CHECK(
        dist_from_center
        == Catch::Approx(tr.trust_region_inflation_radius).epsilon(1e-6));
}

TEST_CASE(
    "planar_filter_step: approaching vertex is truncated", "[ogc][planar_dat]")
{
    // Vertex moving directly toward the triangle should be stopped before
    // penetration. With gap=0.5 and dx_v=(0,-0.6,0), the vertex would
    // penetrate if not truncated.
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_fv_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);
    // Vertex 3 approaches the triangle (would overshoot if not truncated)
    dx.row(3) << 0.0, -0.6, 0.0;

    tr.planar_filter_step(mesh, x, dx, collisions);

    // After truncation the vertex must still be above y=0 (no penetration)
    const double final_y = x(3, 1) + dx(3, 1);
    CHECK(final_y > 0.0);
    // The vertex must have been restricted (dx_y must be less in magnitude than
    // 0.6)
    CHECK(std::abs(dx(3, 1)) < 0.6);
    // And it must be moving in the correct direction (still downward)
    CHECK(dx(3, 1) < 0.0);
}

TEST_CASE(
    "planar_filter_step: preserves more displacement than isotropic for "
    "separating pair",
    "[ogc][planar_dat]")
{
    // When the vertex moves AWAY with a step that exceeds trust_region_radii
    // (tight, per-vertex) but stays within trust_region_inflation_radius
    // (loose), planar_filter_step passes the full step while filter_step
    // clips it to trust_region_radii.
    //
    // For gap=0.5, dhat=0.6: inflation_radius = 1.2,
    // trust_region_radii(3) ≈ gap * relaxed_radius_scaling / 2 ≈ 0.225.
    // A separating step of 0.3 sits between those two values.
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_fv_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.empty());

    // Guard: per-vertex radius must be < 0.3 for the test to be meaningful
    REQUIRE(tr.trust_region_radii(3) < 0.3);

    Eigen::MatrixXd x = V;

    // Planar-DAT: separating step (t_i < 0 — no plane crossing) and within
    // inflation radius → full step preserved
    Eigen::MatrixXd dx_planar = Eigen::MatrixXd::Zero(V.rows(), 3);
    dx_planar.row(3) << 0.0, 0.3, 0.0;
    tr.planar_filter_step(mesh, x, dx_planar, collisions);
    CHECK(dx_planar(3, 1) == Catch::Approx(0.3));

    // Isotropic: clips to trust_region_radii(3) < 0.3
    Eigen::MatrixXd dx_iso = Eigen::MatrixXd::Zero(V.rows(), 3);
    dx_iso.row(3) << 0.0, 0.3, 0.0;
    tr.filter_step(mesh, x, dx_iso);
    CHECK(dx_iso(3, 1) < 0.3 - 1e-10);
}

// Helper: two vertices with no edges/faces, separated by gap along y.
//   V[0] = (0, 0,   0)
//   V[1] = (0, gap, 0)
static std::pair<CollisionMesh, Eigen::MatrixXd> make_vv_mesh(double gap)
{
    Eigen::MatrixXd V(2, 3);
    V << 0.0, 0.0, 0.0, /**/ 0.0, gap, 0.0;
    Eigen::MatrixXi E(0, 2);
    Eigen::MatrixXi F(0, 3);
    return { CollisionMesh(V, E, F), V };
}

// Helper: one edge along x with a vertex above its midpoint, gap along y.
//   V[0] = (0,   0, 0)  edge e0
//   V[1] = (1,   0, 0)  edge e1
//   V[2] = (0.5, gap, 0)  vertex above midpoint
static std::pair<CollisionMesh, Eigen::MatrixXd> make_ev_mesh(double gap)
{
    Eigen::MatrixXd V(3, 3);
    // clang-format off
    V << 0.0,   0,   0.0,
         1.0,   0,   0.0,
         0.5, gap,   0.0;
    // clang-format on
    Eigen::MatrixXi E(1, 2);
    E << 0, 1;
    Eigen::MatrixXi F(0, 3);
    return { CollisionMesh(V, E, F), V };
}

// Helper: 2-D edge (x-axis) with a vertex above midpoint, using 2D positions.
//   V[0] = (0,   0)  edge e0
//   V[1] = (1,   0)  edge e1
//   V[2] = (0.5, gap)  vertex above midpoint
static std::pair<CollisionMesh, Eigen::MatrixXd> make_ev_mesh_2d(double gap)
{
    Eigen::MatrixXd V(3, 2);
    // clang-format off
    V << 0.0,   0,
         1.0,   0,
         0.5, gap;
    // clang-format on
    Eigen::MatrixXi E(1, 2);
    E << 0, 1;
    Eigen::MatrixXi F(0, 3);
    return { CollisionMesh(V, E, F), V };
}

TEST_CASE(
    "planar_filter_step: edge-edge approaching pair is truncated",
    "[ogc][planar_dat]")
{
    // Two perpendicular edges approaching each other should both be truncated.
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_ee_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);
    // Edge A moves up, Edge B moves down — both approach each other
    dx.row(0) << 0.0, 0.4, 0.0;
    dx.row(1) << 0.0, 0.4, 0.0;
    dx.row(2) << 0.0, -0.4, 0.0;
    dx.row(3) << 0.0, -0.4, 0.0;

    tr.planar_filter_step(mesh, x, dx, collisions);

    // After truncation, no penetration: gap between edges must remain > 0
    const double ea_y = x(0, 1) + dx(0, 1); // both ea vertices move same
    const double eb_y = x(2, 1) + dx(2, 1); // both eb vertices move same
    CHECK(eb_y - ea_y > 0.0);               // eb must stay above ea
    // Both edges must have been restricted
    CHECK(std::abs(dx(0, 1)) < 0.4);
    CHECK(std::abs(dx(2, 1)) < 0.4);
}

TEST_CASE(
    "planar_filter_step: vertex-vertex approaching pair is truncated",
    "[ogc][planar_dat]")
{
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_vv_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.vv_collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);
    // v0 approaches v1 (upward), v1 approaches v0 (downward)
    dx.row(0) << 0.0, 0.4, 0.0;
    dx.row(1) << 0.0, -0.4, 0.0;

    tr.planar_filter_step(mesh, x, dx, collisions);

    // After truncation gap must remain positive
    const double final_gap = (x(1, 1) + dx(1, 1)) - (x(0, 1) + dx(0, 1));
    CHECK(final_gap > 0.0);
    // Both vertices must have been restricted
    CHECK(std::abs(dx(0, 1)) < 0.4);
    CHECK(std::abs(dx(1, 1)) < 0.4);
}

TEST_CASE(
    "planar_filter_step: edge-vertex approaching pair is truncated",
    "[ogc][planar_dat]")
{
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_ev_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.ev_collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);
    // Vertex approaches edge (downward), edge stationary
    dx.row(2) << 0.0, -0.6, 0.0;

    tr.planar_filter_step(mesh, x, dx, collisions);

    // After truncation vertex must remain above edge (y > 0)
    const double final_y = x(2, 1) + dx(2, 1);
    CHECK(final_y > 0.0);
    CHECK(std::abs(dx(2, 1)) < 0.6);
    CHECK(dx(2, 1) < 0.0); // still moving downward

    // Tangential motion of the vertex is not truncated by the planar
    // constraint (only the isotropic fallback applies if needed)
    Eigen::MatrixXd dx_tan = Eigen::MatrixXd::Zero(V.rows(), 3);
    dx_tan.row(2) << 0.3, 0.0, 0.0; // tangential (within inflation radius)
    tr.planar_filter_step(mesh, x, dx_tan, collisions);
    CHECK(dx_tan(2, 0) == Catch::Approx(0.3)); // not truncated
}

TEST_CASE(
    "planar_filter_step: approaching vertex, retreating triangle — only "
    "vertex is truncated",
    "[ogc][planar_dat]")
{
    // This tests the image scenario: P approaches triangle ABC while
    // A, B, C move away from P.
    //
    // With delta_t = 0 (triangle not closing), λ = 0 and the division
    // plane is placed at the triangle surface (c_vt). The planar constraint
    // truncates only P; A, B, C have t_i < 0 (moving away from the plane)
    // so their displacements are left untouched.
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_fv_mesh(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.fv_collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 3);

    // Triangle (A=0, B=1, C=2) moves away from P (downward, −y).
    // P (vertex 3) approaches the triangle (would overshoot gap=0.5).
    const double retreat = 0.3;
    const double approach = 0.6;
    dx.row(0) << 0.0, -retreat, 0.0;
    dx.row(1) << 0.0, -retreat, 0.0;
    dx.row(2) << 0.0, -retreat, 0.0;
    dx.row(3) << 0.0, -approach, 0.0; // P towards triangle

    tr.planar_filter_step(mesh, x, dx, collisions);

    // A, B, C move away from the division plane (t_i < 0) → no truncation.
    CHECK(dx(0, 1) == Catch::Approx(-retreat));
    CHECK(dx(1, 1) == Catch::Approx(-retreat));
    CHECK(dx(2, 1) == Catch::Approx(-retreat));

    // P must be truncated: less displacement in magnitude, still downward.
    CHECK(dx(3, 1) > -approach + 1e-10);
    CHECK(dx(3, 1) < 0.0);

    // P must not penetrate the triangle (final P.y > final triangle.y).
    const double final_tri_y = x(0, 1) + dx(0, 1); // all triangle verts same
    const double final_p_y = x(3, 1) + dx(3, 1);
    CHECK(final_p_y > final_tri_y);
}

TEST_CASE(
    "planar_filter_step: works on 2D mesh (edge-vertex)",
    "[ogc][planar_dat]")
{
    const double gap = 0.5;
    const double dhat = 0.6;
    auto [mesh, V] = make_ev_mesh_2d(gap);

    ipc::ogc::TrustRegion tr(dhat);
    ipc::NormalCollisions collisions;
    collisions.set_collision_set_type(
        ipc::NormalCollisions::CollisionSetType::OGC);
    tr.update(mesh, V, collisions);

    REQUIRE(!collisions.ev_collisions.empty());

    Eigen::MatrixXd x = V;
    Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(V.rows(), 2);
    // Vertex approaches edge (downward)
    dx.row(2) << 0.0, -0.6;

    // Must not crash and must truncate the approaching vertex
    tr.planar_filter_step(mesh, x, dx, collisions);

    const double final_y = x(2, 1) + dx(2, 1);
    CHECK(final_y > 0.0);
    CHECK(std::abs(dx(2, 1)) < 0.6);
}
