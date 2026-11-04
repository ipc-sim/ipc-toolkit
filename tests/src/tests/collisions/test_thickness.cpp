#include <tests/utils.hpp>

#include <ipc/ipc.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <igl/edges.h>

using namespace ipc;

namespace {

/// Two triangles facing each other, the upper one gap above the lower and
/// over its interior.
void two_triangles(
    const double gap,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    V.resize(6, 3);
    V << 0, 0, 0, //
        2, 0, 0,  //
        0, 2, 0,  //
        0.2, 0.2, gap, 1.2, 0.2, gap, 0.2, 1.2, gap;
    F.resize(2, 3);
    F << 0, 1, 2, 3, 4, 5;
    igl::edges(F, E);
}

/// Every collision keeps its own two primitives' thicknesses apart.
void check_dmin(const NormalCollisions& collisions, const CollisionMesh& mesh)
{
    for (const auto& c : collisions.vv_collisions) {
        CHECK(c.dmin == Catch::Approx(mesh.pair_thickness(c)));
    }
    for (const auto& c : collisions.ev_collisions) {
        CHECK(c.dmin == Catch::Approx(mesh.pair_thickness(c)));
    }
    for (const auto& c : collisions.ee_collisions) {
        CHECK(c.dmin == Catch::Approx(mesh.pair_thickness(c)));
    }
    for (const auto& c : collisions.fv_collisions) {
        CHECK(c.dmin == Catch::Approx(mesh.pair_thickness(c)));
    }
    for (const auto& c : collisions.pv_collisions) {
        CHECK(c.dmin == Catch::Approx(mesh.pair_thickness(c)));
    }
}

} // namespace

TEST_CASE(
    "A uniform thickness matches a uniform dmin", "[collisions][thickness]")
{
    const auto broad_phase = GENERATE(tests::BroadPhaseGenerator::create());
#ifndef NDEBUG
    if (broad_phase->name() == "BruteForce") {
        SKIP("Brute force is too slow in debug mode");
    }
#endif
    const NormalCollisions::CollisionSetType collision_set_type = GENERATE(
        NormalCollisions::CollisionSetType::IPC,
        NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX,
        NormalCollisions::CollisionSetType::OGC);
    CAPTURE(broad_phase->name(), collision_set_type);

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("two-cubes-close.ply", V, E, F));

    const double dhat = 1e-1;
    const double thickness = 1e-2;

    // Every pair of primitives is 2 * thickness thick, so contacts rest
    // exactly where a dmin of 2 * thickness puts them.
    const CollisionMesh plain_mesh(V, E, F);
    NormalCollisions expected;
    expected.set_collision_set_type(collision_set_type);
    expected.build(plain_mesh, V, dhat, 2 * thickness, broad_phase.get());
    REQUIRE(!expected.empty());

    CollisionMesh mesh(V, E, F);
    mesh.set_vertex_thickness(
        Eigen::VectorXd::Constant(mesh.num_vertices(), thickness));
    NormalCollisions collisions;
    collisions.set_collision_set_type(collision_set_type);
    collisions.build(mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());

    CHECK(collisions.vv_collisions.size() == expected.vv_collisions.size());
    CHECK(collisions.ev_collisions.size() == expected.ev_collisions.size());
    CHECK(collisions.ee_collisions.size() == expected.ee_collisions.size());
    CHECK(collisions.fv_collisions.size() == expected.fv_collisions.size());
    check_dmin(collisions, mesh);

    const BarrierPotential barrier(dhat, /*stiffness=*/1.0);
    // The improved max approximation's sum of positive and negative terms
    // already differs between two builds of the same set (~1e-5 in energy and
    // ~1e-3 in gradient, as its merge order is not deterministic), so it is
    // compared more loosely.
    const double tol = collision_set_type
            == NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX
        ? 1e-2
        : 1e-8;
    CHECK(
        barrier(collisions, mesh, V)
        == Catch::Approx(barrier(expected, plain_mesh, V)).epsilon(tol));
    CHECK(barrier.gradient(collisions, mesh, V)
              .isApprox(barrier.gradient(expected, plain_mesh, V), tol));
}

TEST_CASE(
    "Contacts rest at the sum of their thicknesses", "[collisions][thickness]")
{
    const auto broad_phase = GENERATE(tests::BroadPhaseGenerator::create());
    CAPTURE(broad_phase->name());

    const double dhat = 0.25;
    const double lower = 0.2, upper = 0.1; // face thicknesses

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;

    SECTION("Within dhat of the thicknesses")
    {
        two_triangles(/*gap=*/0.5, V, E, F);

        const CollisionMesh plain_mesh(V, E, F);
        NormalCollisions plain;
        plain.build(
            plain_mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());
        CHECK(plain.empty()); // 0.5 > dhat

        CollisionMesh mesh(V, E, F);
        mesh.set_face_thickness(Eigen::Vector2d(lower, upper));
        NormalCollisions collisions;
        collisions.build(
            mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());
        REQUIRE(!collisions.fv_collisions.empty());
        check_dmin(collisions, mesh);
        for (const auto& fv : collisions.fv_collisions) {
            CHECK(fv.dmin == Catch::Approx(lower + upper));
        }

        const BarrierPotential barrier(dhat, /*stiffness=*/1.0);
        CHECK(barrier(collisions, mesh, V) > 0);

        // A dmin (deprecated) overrides the thickness.
        collisions.build(mesh, V, dhat, /*dmin=*/0.0, broad_phase.get());
        CHECK(collisions.empty());
        collisions.build(mesh, V, dhat, /*dmin=*/0.4, broad_phase.get());
        REQUIRE(!collisions.empty());
        for (size_t i = 0; i < collisions.size(); i++) {
            CHECK(collisions[i].dmin == 0.4);
        }
    }

    SECTION("Beyond dhat of the thicknesses")
    {
        two_triangles(/*gap=*/lower + upper + dhat + 1e-3, V, E, F);

        CollisionMesh mesh(V, E, F);
        mesh.set_face_thickness(Eigen::Vector2d(lower, upper));
        NormalCollisions collisions;
        collisions.build(
            mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());
        CHECK(collisions.empty());
    }

    SECTION("Each primitive type keeps its own thickness")
    {
        two_triangles(/*gap=*/0.5, V, E, F);

        CollisionMesh mesh(V, E, F);
        Eigen::VectorXd vertex(6), edge(E.rows()), face(2);
        vertex << 0.12, 0.12, 0.12, 0.1, 0.1, 0.1;
        for (int i = 0; i < E.rows(); i++) {
            edge[i] = E(i, 0) < 3 ? 0.15 : 0.1;
        }
        face << 0.2, 0.1;
        mesh.set_thickness(vertex, edge, face);

        NormalCollisions collisions;
        collisions.build(
            mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());
        REQUIRE(!collisions.fv_collisions.empty());
        check_dmin(collisions, mesh);
    }
}

TEST_CASE("Codim. vertices with a thickness", "[collisions][thickness][codim]")
{
    const auto broad_phase = GENERATE(tests::BroadPhaseGenerator::create());
    CAPTURE(broad_phase->name());

    Eigen::MatrixXd V(2, 3);
    V << 0, 0, 0, 1.2, 0, 0;
    CollisionMesh mesh(V);
    REQUIRE(mesh.num_codim_vertices() == 2);

    const double dhat = 0.3;
    NormalCollisions collisions;
    collisions.build(mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());
    CHECK(collisions.empty());

    mesh.set_vertex_thickness(Eigen::Vector2d(0.5, 0.5));
    collisions.build(mesh, V, dhat, /*dmin=*/std::nullopt, broad_phase.get());
    REQUIRE(collisions.vv_collisions.size() == 1);
    CHECK(collisions.vv_collisions[0].dmin == Catch::Approx(1.0));
}

TEST_CASE("Plane-vertex with a thickness", "[collisions][thickness][plane]")
{
    Eigen::MatrixXd V(1, 3);
    V << 0, 0.25, 0;
    CollisionMesh mesh(V);
    mesh.planes.emplace_back(Eigen::Vector3d::UnitY(), 0);

    const double dhat = 0.15;
    NormalCollisions collisions;
    collisions.build(mesh, V, dhat);
    CHECK(collisions.empty());

    mesh.set_vertex_thickness(Eigen::VectorXd::Constant(1, 0.2));
    collisions.build(mesh, V, dhat);
    REQUIRE(collisions.pv_collisions.size() == 1);
    CHECK(collisions.pv_collisions[0].dmin == Catch::Approx(0.2));

    V(0, 1) = 0.4;
    collisions.build(mesh, V, dhat);
    CHECK(collisions.empty());
}

TEST_CASE("CCD stops at the thickness", "[ccd][thickness]")
{
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(NDEBUG)
    SKIP("'CCD stops at the thickness' test is skipped in debug mode");
#endif

    const auto broad_phase = GENERATE(tests::BroadPhaseGenerator::create());
    CAPTURE(broad_phase->name());

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    two_triangles(/*gap=*/1, V0, E, F);
    V1 = V0;
    V1.bottomRows(3).col(2).setConstant(-1); // the upper one falls through

    CollisionMesh mesh(V0, E, F);
    const double thickness = 0.1; // per face, so a pair rests 0.2 apart

    const double plain_step = compute_collision_free_stepsize(
        mesh, V0, V1, /*min_distance=*/std::nullopt, broad_phase.get());
    CHECK(plain_step <= 0.5);
    CHECK(plain_step > 0.3);

    mesh.set_face_thickness(Eigen::Vector2d::Constant(thickness));
    const double step = compute_collision_free_stepsize(
        mesh, V0, V1, /*min_distance=*/std::nullopt, broad_phase.get());
    // Contact at a gap of 0.2, so z = 0.2 at t = 0.4.
    CHECK(step <= 0.4);
    CHECK(step > 0.3);
    CHECK(step < plain_step);

    // A min_distance (deprecated) overrides the thickness.
    CHECK(
        compute_collision_free_stepsize(
            mesh, V0, V1, /*min_distance=*/0.0, broad_phase.get())
        == Catch::Approx(plain_step));

    Eigen::MatrixXd V_close = V0;
    V_close.bottomRows(3).col(2).setConstant(0.15);
    CHECK(!is_step_collision_free(
        mesh, V0, V_close, /*min_distance=*/std::nullopt, broad_phase.get()));
    CHECK(is_step_collision_free(
        mesh, V0, V_close, /*min_distance=*/0.0, broad_phase.get()));
    mesh.clear_thickness();
    CHECK(is_step_collision_free(
        mesh, V0, V_close, /*min_distance=*/std::nullopt, broad_phase.get()));
}

TEST_CASE(
    "Per-vertex safe distances with a thickness", "[candidates][thickness]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    two_triangles(/*gap=*/0.5, V, E, F);

    CollisionMesh mesh(V, E, F);
    mesh.set_face_thickness(Eigen::Vector2d(0.2, 0.1));

    const double inflation_radius = 0.5;
    Candidates candidates;
    candidates.build(mesh, V, inflation_radius);
    REQUIRE(!candidates.empty());

    const Eigen::VectorXd safe =
        candidates.compute_per_vertex_safe_distances(mesh, V, inflation_radius);
    // Vertex 3 is 0.5 above the lower face's interior and the pair is 0.3
    // thick.
    CHECK(safe[3] == Catch::Approx(0.2));
    CHECK(safe.minCoeff() >= 0);
}

TEST_CASE("Friction keeps each collision's thickness", "[friction][thickness]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    two_triangles(/*gap=*/0.5, V, E, F);

    CollisionMesh mesh(V, E, F);
    mesh.set_face_thickness(Eigen::Vector2d(0.2, 0.1));

    const double dhat = 0.25;
    const BarrierPotential barrier(dhat, /*stiffness=*/1.0);
    NormalCollisions collisions;
    collisions.build(mesh, V, dhat);
    REQUIRE(!collisions.empty());

    TangentialCollisions tangential_collisions;
    tangential_collisions.build(mesh, V, collisions, barrier, /*mu=*/0.5);
    // Mollified edge-edge pairs have no friction, so there may be fewer.
    // Every pair here is 0.2 + 0.1 thick.
    REQUIRE(!tangential_collisions.empty());
    for (size_t i = 0; i < tangential_collisions.size(); i++) {
        CHECK(tangential_collisions[i].dmin == Catch::Approx(0.3));
        CHECK(tangential_collisions[i].normal_force_magnitude > 0);
    }

    // Without a dmin, the friction force uses each collision's own, so it
    // matches passing the pairs' (uniform) thickness explicitly.
    const FrictionPotential friction(/*eps_v=*/1e-3);
    const Eigen::MatrixXd U = Eigen::MatrixXd::Zero(V.rows(), V.cols());
    Eigen::MatrixXd velocities = U;
    velocities.col(0).setConstant(1e-4);
    const Eigen::VectorXd force =
        friction.force(tangential_collisions, mesh, V, U, velocities, barrier);
    CHECK(force.squaredNorm() > 0);
    CHECK(force.isApprox(friction.force(
        tangential_collisions, mesh, V, U, velocities, barrier,
        /*dmin=*/0.3)));
    CHECK(!force.isApprox(friction.force(
        tangential_collisions, mesh, V, U, velocities, barrier,
        /*dmin=*/0.0)));
}

TEST_CASE(
    "Semi-implicit stiffness keeps each stencil's thickness",
    "[barrier][thickness]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    two_triangles(/*gap=*/0.5, V, E, F);

    CollisionMesh mesh(V, E, F);
    mesh.set_face_thickness(Eigen::Vector2d(0.2, 0.1));

    NormalCollisions collisions;
    collisions.build(mesh, V, /*dhat=*/0.25);
    REQUIRE(!collisions.empty());
    Candidates candidates;
    candidates.build(mesh, V, /*inflation_radius=*/0.125);
    REQUIRE(!candidates.empty());

    const Eigen::VectorXd masses = Eigen::VectorXd::Ones(V.rows());
    Eigen::SparseMatrix<double> hess(V.size(), V.size());
    hess.setIdentity();

    // Every pair here is 0.2 + 0.1 thick, so the default matches passing it.
    CHECK(semi_implicit_stiffness(mesh, V, collisions, masses, hess)
              .isApprox(semi_implicit_stiffness(
                  mesh, V, collisions, masses, hess, /*dmin=*/0.3)));
    CHECK(semi_implicit_stiffness(mesh, V, candidates, masses, hess)
              .isApprox(semi_implicit_stiffness(
                  mesh, V, candidates, masses, hess, /*dmin=*/0.3)));
    CHECK(!semi_implicit_stiffness(mesh, V, candidates, masses, hess)
               .isApprox(semi_implicit_stiffness(
                   mesh, V, candidates, masses, hess, /*dmin=*/0.0)));
}
