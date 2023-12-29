#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <ipc/ipc.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>
#include <ipc/ccd/point_static_plane.hpp>

using namespace ipc;

#ifdef NDEBUG
TEST_CASE("Repeated CCD", "[ccd][repeat]")
#else
TEST_CASE("Repeated CCD", "[ccd][repeat][.]")
#endif
{
    constexpr double FIRST_TOL = 1e-6, SECOND_TOL = 1e-7;
    constexpr long FIRST_MAX_ITER = 1'000'000, SECOND_MAX_ITER = 1'000'000;
    constexpr double MIN_DISTANCE = 0.0;

    // BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();
    BroadPhaseMethod broadphase_method = BroadPhaseMethod::HASH_GRID;
    double inflation_radius = 0;

    bool recompute_candidates = GENERATE(false, true);

    std::string t0_filename, t1_filename;
    SECTION("tooth")
    {
        t0_filename = "private/ccd-failure/repeated_toi_tooth_0.obj";
        t1_filename = "private/ccd-failure/repeated_toi_tooth_1.obj";
    }
    SECTION("hip")
    {
        t0_filename = "private/ccd-failure/repeated_toi_hip_0.obj";
        t1_filename = "private/ccd-failure/repeated_toi_hip_1.obj";
    }
    SECTION("hip small repeated toi")
    {
        t0_filename = "private/ccd-failure/small_repeated_toi_hip_0.obj";
        t1_filename = "private/ccd-failure/small_repeated_toi_hip_1.obj";
    }
    SECTION("hip inf-repeat 0")
    {
        t0_filename = "private/ccd-failure/inf_repeated_toi_hip_0.obj";
        t1_filename = "private/ccd-failure/inf_repeated_toi_hip_1.obj";
    }
    SECTION("hip inf-repeat 1")
    {
        t0_filename = "private/ccd-failure/s0121.obj";
        t1_filename = "private/ccd-failure/s1121.obj";
    }
    SECTION("hip inf-repeat 2")
    {
        t0_filename = "private/ccd-failure/s0110.obj";
        t1_filename = "private/ccd-failure/s1110.obj";
    }
    SECTION("cloth-ball")
    {
        t0_filename = "cloth_ball92.ply";
        t1_filename = "cloth_ball93.ply";
    }

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    bool success = tests::load_mesh(t0_filename, V0, E, F)
        && tests::load_mesh(t1_filename, V1, E, F);
    if (!success) {
        return;
    }
    // REQUIRE(success);

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    Candidates candidates;
    candidates.build(mesh, V0, V1, inflation_radius, broadphase_method);

    bool has_collisions = !candidates.is_step_collision_free(
        mesh, V0, V1, MIN_DISTANCE, FIRST_TOL, FIRST_MAX_ITER);

    double stepsize = candidates.compute_collision_free_stepsize(
        mesh, V0, V1, MIN_DISTANCE, FIRST_TOL, FIRST_MAX_ITER);

    if (!has_collisions) {
        CHECK(stepsize == 1.0);
        return;
    }

    double collision_free_step_size = stepsize;
    bool has_collisions_repeated;
    double stepsize_repeated;
    do {
        Eigen::MatrixXd Vt = (V1 - V0) * collision_free_step_size + V0;
        // CHECK(!has_intersections(Vt, E, F));

        if (recompute_candidates) {
            candidates.build(mesh, V0, Vt, inflation_radius, broadphase_method);
        }

        has_collisions_repeated = !candidates.is_step_collision_free(
            mesh, V0, Vt, MIN_DISTANCE, SECOND_TOL, SECOND_MAX_ITER);

        stepsize_repeated = candidates.compute_collision_free_stepsize(
            mesh, V0, Vt, MIN_DISTANCE, SECOND_TOL, SECOND_MAX_ITER);

        CAPTURE(
            t0_filename, t1_filename, broadphase_method, recompute_candidates,
            has_collisions, collision_free_step_size, has_collisions_repeated,
            stepsize_repeated);
        CHECK(!has_collisions_repeated);
        CHECK(stepsize_repeated == 1.0);

        collision_free_step_size *= stepsize_repeated;
    } while (has_collisions_repeated && stepsize_repeated != 1.0);
}

TEST_CASE("Point-Plane CCD", "[ccd][point-plane]")
{
    Eigen::Vector3d p_t0 = Eigen::Vector3d::Random();

    Eigen::Vector3d origin = Eigen::Vector3d::Random();
    origin.x() = origin.z() = 0;
    Eigen::Vector3d normal(0, 1, 0);

    Eigen::Vector3d p_t1 = p_t0;
    double t = 0;
    SECTION("Known t")
    {
        t = GENERATE(1e-16, 1e-8, 1e-6, 1e-4, 0.1, 0.5, 0.8, 0.88, 0.9, 1.0);
    }
    SECTION("Random t") { t = GENERATE(take(100, random(0.0, 1.0))); }
    if (t == 0) {
        return;
    }
    p_t1.y() = (origin.y() + (t - 1) * p_t0.y()) / t;

    double toi;
    bool is_colliding =
        point_static_plane_ccd(p_t0, p_t1, origin, normal, toi, 0.99);

    CAPTURE(p_t0.x(), p_t0.y(), p_t0.z(), p_t1.x(), p_t1.y(), p_t1.z());
    CHECK(is_colliding);
    CHECK(toi <= t);
}

TEST_CASE("Squash Tet", "[ccd]")
{
    const double dhat = 1e-3;

    Eigen::MatrixXd rest_vertices(4, 3);
    rest_vertices << 0.0, 0.0, 0.0, //
        1.0, 0.0, 0.0,              //
        0.0, 1.0, 0.0,              //
        0.0, 0.0, 1.0;

    Eigen::MatrixXd deformed_vertices = rest_vertices;
    deformed_vertices(3, 0) += 0.1;
    deformed_vertices(3, 1) -= 0.1;
    deformed_vertices(3, 2) = -0.5 * dhat;

    Eigen::MatrixXi edges(6, 2);
    edges << 0, 1, //
        0, 2,      //
        0, 3,      //
        1, 2,      //
        1, 3,      //
        2, 3;
    Eigen::MatrixXi faces(4, 3);
    faces << 0, 2, 1, //
        0, 1, 3,      //
        0, 3, 2,      //
        1, 2, 3;

    ipc::CollisionMesh mesh =
        ipc::CollisionMesh::build_from_full_mesh(rest_vertices, edges, faces);

    // BENCHMARK("compute toi")
    // {
    logger().debug(
        "toi={}",
        ipc::compute_collision_free_stepsize(
            mesh, rest_vertices, deformed_vertices));
    // };
}

TEST_CASE("Degenerate tolerance", "[ccd]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("height-field.ply", V, E, F));

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V, E, F);

    const double t0 = ipc::compute_collision_free_stepsize(mesh, V, V);

    CHECK(t0 == 1.0);
}