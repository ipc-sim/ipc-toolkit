#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/esp/esp_collisions.hpp>
#include <ipc/esp/esp_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <igl/edges.h>
#include <igl/PI.h>

#include <algorithm>

using namespace ipc;

namespace {

// Two tetrahedra whose equilateral bases form a star: a dynamic one with its
// base at z = gap (vertices 0-3) over an obstacle one with its base at z = 0,
// rotated by 60° (vertices 4-7). Each base's tips lie near the other base's
// edges, and the base edges cross. Same-body primitives are >= 1 apart.
CollisionMesh star_mesh(const std::vector<bool>& obstacle_vertex)
{
    constexpr double gap = 0.05;
    Eigen::MatrixXd V(8, 3);
    for (int i = 0; i < 3; i++) {
        const double a = 2 * igl::PI * i / 3;
        V.row(i) << std::cos(a), std::sin(a), gap;
        V.row(4 + i) << std::cos(a + igl::PI / 3), std::sin(a + igl::PI / 3), 0;
    }
    V.row(3) << 0, 0, gap + 1;
    V.row(7) << 0, 0, -1;
    Eigen::MatrixXi F(8, 3);
    F << 0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3, //
        4, 5, 6, 5, 4, 7, 6, 5, 7, 4, 6, 7;
    Eigen::MatrixXi E;
    igl::edges(F, E);

    return CollisionMesh(
        std::vector<bool>(V.rows(), true), std::vector<bool>(V.rows(), false),
        obstacle_vertex, V, E, F);
}

/// The three vertices and the centroid of each face.
FaceQuadRule vertices_and_centroid_rule()
{
    return {
        { { { 1., 0., 0. } }, 0.25 },
        { { { 0., 1., 0. } }, 0.25 },
        { { { 0., 0., 1. } }, 0.25 },
        { { { 1. / 3, 1. / 3, 1. / 3 } }, 0.25 },
    };
}

const std::vector<bool> LOWER_OBSTACLE = { false, false, false, false,
                                           true,  true,  true,  true };

constexpr double DHAT = 0.6;

// An obstacle whose tip, a vertex shared by several obstacle edges, points at a
// deformable surface from height h. In 3D, the tip is the apex of a tetrahedron
// above the interior edge (0, 1) of a two-triangle sheet. In 2D, it is the apex
// of a coarse Λ above the apex (vertex 0) of a fine V. Vertex 0 is deformable
// in both, and the tip is the last vertex before the obstacle's other vertices.
Eigen::MatrixXd apex_vertices(const int dim, const double h)
{
    Eigen::MatrixXd V;
    if (dim == 3) {
        V.resize(8, 3);
        V << -1, 0, 0, 1, 0, 0, 0, -2, 0, 0, 2, 0,               // sheet
            0, 0, h, 0.3, 0, 3, -0.15, 0.26, 3, -0.15, -0.26, 3; // tetrahedron
    } else {
        V.resize(6, 2);
        V << 0, 0, -0.1, -0.1, 0.1, -0.1, // V
            0, h, -1, 1 + h, 1, 1 + h;    // Λ
    }
    return V;
}

index_t apex_obstacle_tip(const int dim) { return dim == 3 ? 4 : 3; }

CollisionMesh apex_mesh(const int dim)
{
    const Eigen::MatrixXd V = apex_vertices(dim, 0.05);
    std::vector<bool> obstacle(V.rows(), false);
    std::fill(obstacle.begin() + apex_obstacle_tip(dim), obstacle.end(), true);

    Eigen::MatrixXi E, F;
    if (dim == 3) {
        F.resize(6, 3);
        F << 0, 2, 1, 0, 1, 3, //
            4, 5, 6, 4, 6, 7, 4, 7, 5, 5, 7, 6;
        igl::edges(F, E);
    } else {
        E.resize(4, 2);
        E << 1, 0, 0, 2, 4, 3, 3, 5;
    }

    return CollisionMesh(
        std::vector<bool>(V.rows(), true), std::vector<bool>(V.rows(), false),
        obstacle, V, E, F);
}

} // namespace

TEST_CASE("Collision mesh obstacle flags", "[collision_mesh][obstacle]")
{
    SECTION("Default constructor has no obstacles")
    {
        const CollisionMesh mesh = star_mesh(std::vector<bool>(8, false));
        const CollisionMesh plain(
            mesh.rest_positions(), mesh.edges(), mesh.faces());
        for (int i = 0; i < plain.num_vertices(); i++) {
            CHECK(!plain.is_obstacle_vertex(i));
        }
    }

    SECTION("Flags follow the vertices")
    {
        const CollisionMesh mesh = star_mesh(LOWER_OBSTACLE);
        for (int i = 0; i < mesh.num_vertices(); i++) {
            CHECK(mesh.is_obstacle_vertex(i) == (i >= 4));
        }
        for (int e = 0; e < mesh.num_edges(); e++) {
            CHECK(mesh.is_obstacle_edge(e) == (mesh.edges()(e, 0) >= 4));
        }
        for (int f = 0; f < mesh.num_faces(); f++) {
            CHECK(mesh.is_obstacle_face(f) == (f >= 4));
        }
    }

    SECTION("Flags are mapped to the included vertices")
    {
        // Full mesh with an extra, excluded obstacle vertex at index 0.
        const CollisionMesh star = star_mesh(LOWER_OBSTACLE);
        Eigen::MatrixXd V(9, 3);
        V.row(0).setZero();
        V.bottomRows(8) = star.rest_positions();
        const Eigen::MatrixXi E = (star.edges().array() + 1).matrix();
        const Eigen::MatrixXi F = (star.faces().array() + 1).matrix();

        std::vector<bool> include(9, true), obstacle(9, true);
        include[0] = false;
        for (int i = 1; i <= 4; i++) {
            obstacle[i] = false;
        }

        const CollisionMesh mesh(
            include, std::vector<bool>(9, false), obstacle, V, E, F);
        REQUIRE(mesh.num_vertices() == 8);
        for (int i = 0; i < mesh.num_vertices(); i++) {
            CHECK(mesh.is_obstacle_vertex(i) == (i >= 4));
        }
    }

    SECTION("Edges and faces touching a dynamic vertex are dynamic")
    {
        // Vertex 0 of the dynamic body is flagged as an obstacle, like the end
        // of a rope anchored to a static object.
        std::vector<bool> obstacle = LOWER_OBSTACLE;
        obstacle[0] = true;
        const CollisionMesh mesh = star_mesh(obstacle);
        for (int e = 0; e < mesh.num_edges(); e++) {
            CHECK(mesh.is_obstacle_edge(e) == (mesh.edges()(e, 0) >= 4));
        }
        for (int f = 0; f < mesh.num_faces(); f++) {
            CHECK(mesh.is_obstacle_face(f) == (f >= 4));
        }
    }

    SECTION("Mask must cover the full mesh")
    {
        const CollisionMesh star = star_mesh(LOWER_OBSTACLE);
        const Eigen::MatrixXd& V = star.rest_positions();
        const std::vector<bool> all(V.rows(), true);
        for (const size_t size : { V.rows() - 1, V.rows() + 1 }) {
            CHECK_THROWS_AS(
                CollisionMesh(
                    all, all, std::vector<bool>(size, false), V, star.edges(),
                    star.faces()),
                std::invalid_argument);
        }
    }
}

TEST_CASE("Normal collisions skip obstacles", "[collisions][obstacle]")
{
    CollisionMesh mesh = star_mesh(LOWER_OBSTACLE);
    const Eigen::MatrixXd& V = mesh.rest_positions();
    // A ground plane just below the obstacle's tip (vertex 7).
    mesh.planes.emplace_back(
        Eigen::Vector3d::UnitZ(), Eigen::Vector3d(0, 0, -1.1));

    NormalCollisions all;
    all.build(mesh, V, DHAT);

    NormalCollisions skipped;
    skipped.set_skip_obstacles(true);
    CHECK(skipped.skip_obstacles());
    skipped.build(mesh, V, DHAT);

    // Pairs whose vertex is an obstacle are dropped.
    const auto count_obstacle_vertices = [&](const NormalCollisions& c) {
        int n = 0;
        for (const auto& ev : c.ev_collisions) {
            n += mesh.is_obstacle_vertex(ev.vertex_id);
        }
        for (const auto& fv : c.fv_collisions) {
            n += mesh.is_obstacle_vertex(fv.vertex_id);
        }
        for (const auto& pv : c.pv_collisions) {
            n += mesh.is_obstacle_vertex(pv.vertex_id);
        }
        return n;
    };
    REQUIRE(count_obstacle_vertices(all) > 0);
    REQUIRE(!all.pv_collisions.empty());
    CHECK(count_obstacle_vertices(skipped) == 0);
    CHECK(skipped.size() < all.size());

    // Edge-edge pairs with one obstacle edge are kept at half weight.
    REQUIRE(!all.ee_collisions.empty());
    REQUIRE(skipped.ee_collisions.size() == all.ee_collisions.size());
    for (const auto& ee : skipped.ee_collisions) {
        REQUIRE(
            mesh.is_obstacle_edge(ee.edge0_id)
            != mesh.is_obstacle_edge(ee.edge1_id));
        const auto it = std::find_if(
            all.ee_collisions.begin(), all.ee_collisions.end(),
            [&](const auto& other) { return other == ee; });
        REQUIRE(it != all.ee_collisions.end());
        CHECK(ee.weight == Catch::Approx(0.5 * it->weight));
    }

    // With every vertex an obstacle, nothing is left.
    const CollisionMesh all_obstacle = star_mesh(std::vector<bool>(8, true));
    NormalCollisions none;
    none.set_skip_obstacles(true);
    none.build(all_obstacle, V, DHAT);
    CHECK(none.empty());
}

TEST_CASE(
    "Skipping obstacles keeps the barrier repulsive", "[collisions][obstacle]")
{
    using CollisionSetType = NormalCollisions::CollisionSetType;

    const int dim = GENERATE(2, 3);
    const CollisionSetType type =
        GENERATE(CollisionSetType::IPC, CollisionSetType::IMPROVED_MAX_APPROX);
    const bool use_area_weighting = GENERATE(false, true);
    CAPTURE(dim, type, use_area_weighting);

    const CollisionMesh mesh = apex_mesh(dim);
    const index_t tip = apex_obstacle_tip(dim);
    constexpr double dhat = 0.1;
    const BarrierPotential barrier(dhat, 1.0);

    // As the tip approaches the deformable surface, the energy must grow and
    // push the two apart. The improved max approximation's corrections used to
    // outweigh the one-sided contact here, making the barrier attractive.
    double prev_energy = 0;
    for (const double h : { 0.05, 0.01, 1e-3, 1e-5 }) {
        CAPTURE(h);
        const Eigen::MatrixXd V = apex_vertices(dim, h);

        NormalCollisions collisions;
        collisions.set_use_area_weighting(use_area_weighting);
        collisions.set_collision_set_type(type);
        collisions.set_skip_obstacles(true);
        collisions.build(mesh, V, dhat);

        const double energy = barrier(collisions, mesh, V);
        CHECK(energy > prev_energy);
        prev_energy = energy;

        // The last coordinate points from the deformable surface to the tip.
        const Eigen::VectorXd grad = barrier.gradient(collisions, mesh, V);
        CHECK(grad[dim * tip + dim - 1] < 0);
        CHECK(grad[dim - 1] > 0);
    }
}

TEST_CASE(
    "Obstacle size does not change edge-edge contact stiffness",
    "[collisions][obstacle]")
{
    using CollisionSetType = NormalCollisions::CollisionSetType;

    // A short deformable edge crossing under an obstacle edge of varying
    // length.
    const double obstacle_length = GENERATE(0.2, 2.0, 20.0);
    const CollisionSetType type =
        GENERATE(CollisionSetType::IPC, CollisionSetType::IMPROVED_MAX_APPROX);
    CAPTURE(obstacle_length, type);

    constexpr double dhat = 0.05;
    const auto potential_value = [&](const double length, const bool skip) {
        Eigen::MatrixXd V(4, 3);
        V << -0.1, 0, 0, 0.1, 0, 0, //
            0, -length / 2, 0.01, 0, length / 2, 0.01;
        Eigen::MatrixXi E(2, 2);
        E << 0, 1, 2, 3;
        const CollisionMesh mesh(
            std::vector<bool>(4, true), std::vector<bool>(4, false),
            { false, false, true, true }, V, E, Eigen::MatrixXi());

        NormalCollisions collisions;
        collisions.set_use_area_weighting(true);
        collisions.set_collision_set_type(type);
        collisions.set_skip_obstacles(skip);
        collisions.build(mesh, V, dhat);
        return BarrierPotential(dhat, 1.0)(collisions, mesh, V);
    };

    // Two equal edges contribute equally, so integrating only over the
    // deformable one gives half the two-sided energy, whatever the obstacle's
    // size.
    const double two_sided = potential_value(0.2, false);
    REQUIRE(two_sided > 0);
    CHECK(
        potential_value(obstacle_length, true)
        == Catch::Approx(two_sided / 2).epsilon(1e-12));
}

TEST_CASE("ESP integration types", "[esp_potential][obstacle]")
{
    using IntegrationType = ESPParameters::IntegrationType;

    // Vertex integration or a face quadrature rule.
    const bool face_rule = GENERATE(false, true);
    CAPTURE(face_rule);

    const auto potential_value = [&](const std::vector<bool>& obstacle_vertex,
                                     const IntegrationType type) {
        const CollisionMesh mesh = star_mesh(obstacle_vertex);
        ESPParameters params(DHAT, 1.0, 1, true, type);
        if (face_rule) {
            params.set_quad_rule(vertices_and_centroid_rule());
        }
        ESPCollisions collisions;
        collisions.build(mesh, mesh.rest_positions(), params);
        return ESPPotential(params)(collisions, mesh, mesh.rest_positions());
    };

    const std::vector<bool> none(8, false), all(8, true);

    // Without obstacles, every integration type is the same.
    const double expected = potential_value(none, IntegrationType::BRUTE_FORCE);
    REQUIRE(expected > 0);
    CHECK(
        potential_value(none, IntegrationType::NORMAL)
        == Catch::Approx(expected).epsilon(1e-12));
    CHECK(
        potential_value(none, IntegrationType::NO_OBST)
        == Catch::Approx(expected).epsilon(1e-12));

    // BRUTE_FORCE ignores the obstacle flags.
    for (const auto& flags : { LOWER_OBSTACLE, all }) {
        CHECK(
            potential_value(flags, IntegrationType::BRUTE_FORCE)
            == Catch::Approx(expected).epsilon(1e-12));
    }

    // The two bodies are mirror images, so each surface integrates half of the
    // contact; NO_OBST integrates the dynamic one only.
    CHECK(
        potential_value(LOWER_OBSTACLE, IntegrationType::NO_OBST)
        == Catch::Approx(expected / 2).epsilon(1e-12));

    // NORMAL keeps the obstacle's sources but drops their obstacle-obstacle
    // pairs: the base edges cross within dhat of the obstacle's own vertices.
    const double normal =
        potential_value(LOWER_OBSTACLE, IntegrationType::NORMAL);
    CHECK(normal > expected / 2);
    CHECK(normal != Catch::Approx(expected).epsilon(1e-12));

    // Obstacle-only contact is integrated only by BRUTE_FORCE.
    CHECK(potential_value(all, IntegrationType::NORMAL) == 0);
    CHECK(potential_value(all, IntegrationType::NO_OBST) == 0);
}

TEST_CASE("ESP integration types 2D", "[esp_potential][obstacle]")
{
    using IntegrationType = ESPParameters::IntegrationType;

    // Two squares 0.2 apart, mirror images of each other about x = 0.
    Eigen::MatrixXd V(8, 2);
    V << -1, 1, -1, 0, -.1, 0, -.1, 1, .1, 1, .1, 0, 1, 0, 1, 1;
    Eigen::MatrixXi E(8, 2);
    E << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4;

    const auto potential_value = [&](const std::vector<bool>& obstacle_vertex,
                                     const IntegrationType type) {
        const CollisionMesh mesh(
            std::vector<bool>(V.rows(), true),
            std::vector<bool>(V.rows(), false), obstacle_vertex, V, E,
            Eigen::MatrixXi());
        const ESPParameters params(0.6, 1.0, 2, true, type);
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        return ESPPotential(params)(collisions, mesh, V);
    };

    const std::vector<bool> none(8, false), all(8, true);
    const std::vector<bool> left = { true,  true,  true,  true,
                                     false, false, false, false };

    const double expected = potential_value(none, IntegrationType::BRUTE_FORCE);
    REQUIRE(expected > 0);
    for (const auto type :
         { IntegrationType::NORMAL, IntegrationType::NO_OBST }) {
        CHECK(
            potential_value(none, type)
            == Catch::Approx(expected).epsilon(1e-12));
    }
    for (const auto& flags : { left, all }) {
        CHECK(
            potential_value(flags, IntegrationType::BRUTE_FORCE)
            == Catch::Approx(expected).epsilon(1e-12));
    }

    // NO_OBST integrates the dynamic square only.
    CHECK(
        potential_value(left, IntegrationType::NO_OBST)
        == Catch::Approx(expected / 2).epsilon(1e-12));
    // The obstacle-obstacle terms NORMAL drops sum to zero here.
    CHECK(
        potential_value(left, IntegrationType::NORMAL)
        == Catch::Approx(expected).epsilon(1e-12));

    CHECK(potential_value(all, IntegrationType::NORMAL) == 0);
    CHECK(potential_value(all, IntegrationType::NO_OBST) == 0);
}

TEST_CASE(
    "ESP obstacle vertices next to dynamic vertices",
    "[esp_potential][obstacle]")
{
    using IntegrationType = ESPParameters::IntegrationType;

    // Two unit cubes 0.001 apart with aligned vertices; the first is an
    // obstacle, so its vertices have dynamic vertex-vertex candidates.
    Eigen::MatrixXd V_cube;
    Eigen::MatrixXi E_cube, F_cube;
    REQUIRE(tests::load_mesh("cube.ply", V_cube, E_cube, F_cube));
    const int n = V_cube.rows();
    Eigen::MatrixXd V(2 * n, 3);
    V << V_cube, V_cube.rowwise() + Eigen::RowVector3d(1.001, 0, 0);
    Eigen::MatrixXi F(2 * F_cube.rows(), 3);
    F << F_cube, (F_cube.array() + n).matrix();
    Eigen::MatrixXi E;
    igl::edges(F, E);
    std::vector<bool> obstacle(2 * n, false);
    std::fill_n(obstacle.begin(), n, true);
    const CollisionMesh mesh(
        std::vector<bool>(2 * n, true), std::vector<bool>(2 * n, false),
        obstacle, V, E, F);

    const bool face_rule = GENERATE(false, true);
    CAPTURE(face_rule);
    const auto potential_value = [&](const IntegrationType type) {
        ESPParameters params(0.5, 1.0, 0, true, type);
        if (face_rule) {
            params.set_quad_rule(vertices_and_centroid_rule());
        }
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        return ESPPotential(params)(collisions, mesh, V);
    };

    const double expected = potential_value(IntegrationType::BRUTE_FORCE);
    REQUIRE(expected > 0);
    // No obstacle-obstacle pair is within dhat, so NORMAL drops nothing.
    CHECK(
        potential_value(IntegrationType::NORMAL)
        == Catch::Approx(expected).epsilon(1e-12));
    CHECK(
        potential_value(IntegrationType::NO_OBST)
        == Catch::Approx(expected / 2).epsilon(1e-12));
}

TEST_CASE(
    "ESP integrates edges and faces with mixed obstacle flags",
    "[esp_potential][obstacle]")
{
    using IntegrationType = ESPParameters::IntegrationType;

    const int dim = GENERATE(2, 3);
    CAPTURE(dim);

    // In 3D, the star with one dynamic vertex flagged as an obstacle. In 2D, a
    // rope hanging from the top right corner of a static block, 0.05 from its
    // side; the rope's first edge has one obstacle and one dynamic vertex.
    CollisionMesh mesh;
    double dhat;
    if (dim == 3) {
        std::vector<bool> obstacle = LOWER_OBSTACLE;
        obstacle[0] = true;
        mesh = star_mesh(obstacle);
        dhat = DHAT;
    } else {
        Eigen::MatrixXd V(6, 2);
        V << -0.5, 0, -0.5, -1, 0.5, -1, 0.5, 0, // block
            0.55, -0.3, 0.55, -0.7;              // rope
        Eigen::MatrixXi E(6, 2);
        E << 0, 1, 1, 2, 2, 3, 3, 0, 3, 4, 4, 5;
        mesh = CollisionMesh(
            std::vector<bool>(6, true), std::vector<bool>(6, false),
            { true, true, true, true, false, false }, V, E, Eigen::MatrixXi());
        dhat = 0.1;
    }
    const Eigen::MatrixXd& V = mesh.rest_positions();

    const auto potential_value = [&](const CollisionMesh& m,
                                     const IntegrationType type) {
        const ESPParameters params(dhat, 1.0, 1, true, type);
        ESPCollisions collisions;
        collisions.build(m, V, params);
        return ESPPotential(params)(collisions, m, V);
    };

    // BRUTE_FORCE ignores the flags, so it matches the unflagged mesh.
    const CollisionMesh unflagged(V, mesh.edges(), mesh.faces());
    const double expected =
        potential_value(unflagged, IntegrationType::BRUTE_FORCE);
    REQUIRE(expected > 0);
    CHECK(
        potential_value(mesh, IntegrationType::BRUTE_FORCE)
        == Catch::Approx(expected).epsilon(1e-12));

    // The mixed primitives move, so the filtered types integrate them too.
    for (const auto type :
         { IntegrationType::NORMAL, IntegrationType::NO_OBST }) {
        CAPTURE(type);
        double value = 0;
        REQUIRE_NOTHROW(value = potential_value(mesh, type));
        CHECK(value > 0);
    }
}

TEST_CASE(
    "ESP NORMAL integration sees dynamic edges near obstacle corners in 2D",
    "[esp_potential][obstacle]")
{
    using IntegrationType = ESPParameters::IntegrationType;

    // A dynamic box resting on a static obstacle whose top lies under the
    // box's bottom edge but has no vertex near it: a spike pointing at the
    // edge, or a block narrower than the box. Only the obstacle's corner
    // vertices are close to the box, so the box's edge reaches the obstacle's
    // edges through those corners.
    const bool spike = GENERATE(true, false);
    const double gap = GENERATE(0.05, 0.01, 1e-3, 1e-6);
    const int quad_order = GENERATE(1, 2);
    CAPTURE(spike, gap, quad_order);

    Eigen::MatrixXd V(spike ? 7 : 8, 2);
    V.topRows(4) << -1, gap, 1, gap, 1, gap + 1, -1, gap + 1;
    Eigen::MatrixXi E(spike ? 7 : 8, 2);
    E.topRows(4) << 0, 1, 1, 2, 2, 3, 3, 0;
    if (spike) {
        V.bottomRows(3) << 0, 0, -0.5, -1, 0.5, -1;
        E.bottomRows(3) << 4, 5, 5, 6, 6, 4;
    } else {
        V.bottomRows(4) << -0.5, 0, -0.5, -1, 0.5, -1, 0.5, 0;
        E.bottomRows(4) << 4, 5, 5, 6, 6, 7, 7, 4;
    }
    std::vector<bool> obstacle(V.rows(), true);
    std::fill_n(obstacle.begin(), 4, false);
    const CollisionMesh mesh(
        std::vector<bool>(V.rows(), true), std::vector<bool>(V.rows(), false),
        obstacle, V, E, Eigen::MatrixXi());

    const auto potential_value = [&](const IntegrationType type) {
        const ESPParameters params(0.1, 0.2, quad_order, true, type);
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        return ESPPotential(params)(collisions, mesh, V);
    };

    // No two obstacle primitives are within dhat of each other, so NORMAL has
    // nothing to drop and must match BRUTE_FORCE.
    const double expected = potential_value(IntegrationType::BRUTE_FORCE);
    REQUIRE(expected > 0);
    CHECK(
        potential_value(IntegrationType::NORMAL)
        == Catch::Approx(expected).epsilon(1e-12));
}
