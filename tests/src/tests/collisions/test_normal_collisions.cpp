#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <igl/edges.h>

using namespace ipc;

TEST_CASE("Codim. vertex-vertex collisions", "[collisions][codim]")
{
    constexpr double thickness = 0.4;
    constexpr double min_distance = 2 * thickness;

    Eigen::MatrixXd vertices(8, 3);
    vertices << 0, 0, 0, //
        0, 0, 1,         //
        0, 1, 0,         //
        0, 1, 1,         //
        1, 0, 0,         //
        1, 0, 1,         //
        1, 1, 0,         //
        1, 1, 1;
    vertices.rowwise() -= vertices.colwise().mean();

    CollisionMesh mesh(vertices);
    mesh.init_area_jacobians();

    CHECK(mesh.num_vertices() == 8);
    CHECK(mesh.num_codim_vertices() == 8);
    CHECK(mesh.num_edges() == 0);
    CHECK(mesh.num_faces() == 0);

    auto broad_phase = GENERATE(tests::BroadPhaseGenerator::create());
    CAPTURE(broad_phase->name());

    SECTION("Candidates")
    {
        Eigen::MatrixXd V1 = vertices;
        V1.col(1) *= 0.5;

        Candidates candidates;
        candidates.build(mesh, vertices, V1, thickness, broad_phase);

        CHECK(!candidates.empty());
        CHECK(candidates.vv_candidates.size() == candidates.size());

        CHECK(!candidates.is_step_collision_free(
            mesh, vertices, V1, min_distance));

        // Account for conservative rescaling
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
        constexpr double conservative_min_dist = 0.2 * (1 - min_distance);
#else
        constexpr double conservative_min_dist = 1e-4;
#endif
        constexpr double expected_toi =
            (1 - (min_distance + conservative_min_dist)) / 2.0 / 0.25;
        CHECK(
            candidates.compute_collision_free_stepsize(
                mesh, vertices, V1, min_distance)
            == Catch::Approx(expected_toi));
    }

    SECTION("Collisions")
    {
        const bool use_area_weighting = GENERATE(false, true);
        const NormalCollisions::CollisionSetType collision_set_type = GENERATE(
            NormalCollisions::CollisionSetType::IPC,
            NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX,
            NormalCollisions::CollisionSetType::OGC);
        const bool use_physical_barrier = GENERATE(false, true);
        const bool enable_shape_derivatives = GENERATE(false, true);
        const double dhat = 0.25;

        NormalCollisions collisions;
        collisions.set_use_area_weighting(use_area_weighting);
        collisions.set_collision_set_type(collision_set_type);
        collisions.set_enable_shape_derivatives(enable_shape_derivatives);

        collisions.build(mesh, vertices, dhat, min_distance, broad_phase);

        CHECK(collisions.size() == 12);
        CHECK(collisions.vv_collisions.size() == 12);

        BarrierPotential barrier_potential(dhat, use_physical_barrier);

        CHECK(barrier_potential(collisions, mesh, vertices) > 0.0);
        const Eigen::VectorXd grad =
            barrier_potential.gradient(collisions, mesh, vertices);
        for (int i = 0; i < vertices.rows(); i++) {
            const Eigen::Vector3d f = -grad.segment<3>(3 * i);
            CHECK(f.normalized().isApprox(
                vertices.row(i).normalized().transpose()));
        }
    }
}

TEST_CASE("Codim. edge-vertex collisions", "[collisions][codim]")
{
    constexpr double thickness = 1e-3;
    constexpr double min_distance = 2 * thickness;

    Eigen::MatrixXd vertices(8, 3);
    vertices << 0, 0, 0, //
        1, 0, 0,         //
        0, 0, -1,        //
        -1, 0, 0,        //
        0, 0, 1,         //
        0, 1, 0,         //
        0, 2, 0,         //
        0, 3, 0;
    Eigen::MatrixXi edges(4, 2);
    edges << 0, 1, //
        0, 2,      //
        0, 3,      //
        0, 4;

    CollisionMesh mesh(vertices, edges);
    mesh.init_area_jacobians();

    CHECK(mesh.num_vertices() == 8);
    CHECK(mesh.num_codim_vertices() == 3);
    CHECK(mesh.num_codim_edges() == 4);
    CHECK(mesh.num_edges() == 4);
    CHECK(mesh.num_faces() == 0);

    auto broad_phase = GENERATE(tests::BroadPhaseGenerator::create());
    CAPTURE(broad_phase->name());

    SECTION("Candidates")
    {
        Eigen::MatrixXd V1 = vertices;
        V1.bottomRows(3).col(1).array() -= 4; // Translate the codim vertices

        Candidates candidates;
        candidates.build(mesh, vertices, V1, thickness, broad_phase);

        CHECK(candidates.size() == 15);
        CHECK(candidates.vv_candidates.size() == 3);
        CHECK(candidates.ev_candidates.size() == 12);
        CHECK(candidates.ee_candidates.size() == 0);
        CHECK(candidates.fv_candidates.size() == 0);

        CHECK(!candidates.is_step_collision_free(
            mesh, vertices, V1, min_distance));

        // Account for conservative rescaling
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
        constexpr double conservative_min_dist = 0.2 * (1 - min_distance);
#else
        constexpr double conservative_min_dist = 1e-4;
#endif
        constexpr double expected_toi =
            (1 - (min_distance + conservative_min_dist)) / 4;
        CHECK(
            candidates.compute_collision_free_stepsize(
                mesh, vertices, V1, min_distance)
            == Catch::Approx(expected_toi));
    }

    SECTION("Collisions")
    {
        const bool use_area_weighting = GENERATE(false, true);
        const NormalCollisions::CollisionSetType collision_set_type = GENERATE(
            NormalCollisions::CollisionSetType::IPC,
            NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX,
            NormalCollisions::CollisionSetType::OGC);
        const bool use_physical_barrier = GENERATE(false, true);
        const bool enable_shape_derivatives = GENERATE(false, true);

        NormalCollisions collisions;
        collisions.set_use_area_weighting(use_area_weighting);
        collisions.set_collision_set_type(collision_set_type);
        collisions.set_enable_shape_derivatives(enable_shape_derivatives);

        const double dhat = 0.25;
        collisions.build(
            mesh, vertices, dhat, /*min_distance=*/0.8, broad_phase);

        const int expected_num_collisions = 6
            + int(collision_set_type
                  == NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX);
        const int expected_num_vv_collisions = 2
            + int(collision_set_type
                  == NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX);

        CHECK(collisions.size() == expected_num_collisions);
        CHECK(collisions.vv_collisions.size() == expected_num_vv_collisions);
        CHECK(collisions.ev_collisions.size() == 4);
        CHECK(collisions.ee_collisions.size() == 0);
        CHECK(collisions.fv_collisions.size() == 0);

        CHECK(
            BarrierPotential(dhat, use_physical_barrier)(
                collisions, mesh, vertices)
            > 0.0);
    }
}

TEST_CASE("Vertex-Vertex NormalCollision", "[collision][vertex-vertex]")
{
    CHECK(
        VertexVertexNormalCollision(0, 1) == VertexVertexNormalCollision(0, 1));
    CHECK(
        VertexVertexNormalCollision(0, 1) == VertexVertexNormalCollision(1, 0));
    CHECK(
        VertexVertexNormalCollision(0, 1) != VertexVertexNormalCollision(0, 2));
    CHECK(
        VertexVertexNormalCollision(0, 1) != VertexVertexNormalCollision(2, 0));
    CHECK(
        VertexVertexNormalCollision(0, 1) < VertexVertexNormalCollision(0, 2));
    CHECK(
        VertexVertexNormalCollision(0, 1) < VertexVertexNormalCollision(2, 0));

    CHECK(
        VertexVertexNormalCollision(VertexVertexCandidate(0, 1))
        == VertexVertexNormalCollision(0, 1));
}

TEST_CASE("Edge-Vertex NormalCollision", "[collision][edge-vertex]")
{
    CHECK(EdgeVertexNormalCollision(0, 1) == EdgeVertexNormalCollision(0, 1));
    CHECK(EdgeVertexNormalCollision(0, 1) != EdgeVertexNormalCollision(1, 0));
    CHECK(EdgeVertexNormalCollision(0, 1) != EdgeVertexNormalCollision(0, 2));
    CHECK(EdgeVertexNormalCollision(0, 1) != EdgeVertexNormalCollision(2, 0));
    CHECK(EdgeVertexNormalCollision(0, 1) < EdgeVertexNormalCollision(0, 2));
    CHECK(!(EdgeVertexNormalCollision(1, 1) < EdgeVertexNormalCollision(0, 2)));
    CHECK(EdgeVertexNormalCollision(0, 1) < EdgeVertexNormalCollision(2, 0));

    CHECK(
        EdgeVertexNormalCollision(EdgeVertexCandidate(0, 1))
        == EdgeVertexNormalCollision(0, 1));
}

TEST_CASE("Edge-Edge NormalCollision", "[collision][edge-edge]")
{
    CHECK(
        EdgeEdgeNormalCollision(0, 1, 0.0)
        == EdgeEdgeNormalCollision(0, 1, 1.0));
    CHECK(
        EdgeEdgeNormalCollision(0, 1, 0.0)
        == EdgeEdgeNormalCollision(1, 0, 1.0));
    CHECK(
        EdgeEdgeNormalCollision(0, 1, 0.0)
        != EdgeEdgeNormalCollision(0, 2, 1.0));
    CHECK(
        EdgeEdgeNormalCollision(0, 1, 0.0)
        != EdgeEdgeNormalCollision(2, 0, 1.0));
    CHECK(
        EdgeEdgeNormalCollision(0, 1, 0.0)
        < EdgeEdgeNormalCollision(0, 2, 1.0));
    CHECK(
        EdgeEdgeNormalCollision(0, 1, 0.0)
        < EdgeEdgeNormalCollision(2, 0, 1.0));

    CHECK(
        EdgeEdgeNormalCollision(EdgeEdgeCandidate(0, 1), 0.0)
        == EdgeEdgeNormalCollision(0, 1, 0.0));
}

TEST_CASE("Face-Vertex NormalCollision", "[collision][face-vertex]")
{
    CHECK(FaceVertexNormalCollision(0, 1) == FaceVertexNormalCollision(0, 1));
    CHECK(FaceVertexNormalCollision(0, 1) != FaceVertexNormalCollision(1, 0));
    CHECK(FaceVertexNormalCollision(0, 1) != FaceVertexNormalCollision(0, 2));
    CHECK(FaceVertexNormalCollision(0, 1) != FaceVertexNormalCollision(2, 0));
    CHECK(FaceVertexNormalCollision(0, 1) < FaceVertexNormalCollision(0, 2));
    CHECK(!(FaceVertexNormalCollision(1, 1) < FaceVertexNormalCollision(0, 2)));
    CHECK(FaceVertexNormalCollision(0, 1) < FaceVertexNormalCollision(2, 0));

    CHECK(
        FaceVertexNormalCollision(FaceVertexCandidate(0, 1))
        == FaceVertexNormalCollision(0, 1));
}

TEST_CASE("Plane-Vertex NormalCollision", "[collision][plane-vertex]")
{
    Eigen::MatrixXi edges, faces;
    const Eigen::Vector3d n(0, 1, 0), o(0, 0, 0);
    const PlaneVertexNormalCollision c(o, n, 0);
    CHECK(c.num_vertices() == 1);
    CHECK(
        c.vertex_ids(edges, faces)
        == std::array<index_t, 4> { { 0, -1, -1, -1 } });
    CHECK(c.plane_origin == o);
    CHECK(c.plane_normal == n);
    CHECK(c.vertex_id == 0);

    CHECK(c.compute_distance(Eigen::Vector3d(0, -2, 0)) == 4.0);
    CHECK(c.compute_distance(Eigen::Vector3d(0, 2, 0)) == 4.0);
    CHECK(
        c.compute_distance_gradient(Eigen::Vector3d(0, 2, 0))
        == Eigen::Vector3d(0, 4, 0));
    CHECK(
        c.compute_distance_hessian(Eigen::Vector3d(0, 2, 0))
        == 2 * n * n.transpose());
}

TEST_CASE("NormalCollisions::is_*", "[collisions]")
{
    NormalCollisions collisions;
    collisions.vv_collisions.emplace_back(0, 1);
    collisions.ev_collisions.emplace_back(0, 1);
    collisions.ee_collisions.emplace_back(0, 1, 0.0);
    collisions.fv_collisions.emplace_back(0, 1);
    collisions.pv_collisions.emplace_back(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 0), 0);

    for (int i = 0; i < collisions.size(); i++) {
        CHECK(collisions.is_vertex_vertex(i) == (i == 0));
        CHECK(collisions.is_edge_vertex(i) == (i == 1));
        CHECK(collisions.is_edge_edge(i) == (i == 2));
        CHECK(collisions.is_face_vertex(i) == (i == 3));
        CHECK(collisions.is_plane_vertex(i) == (i == 4));
    }
}

TEST_CASE("NormalCollisions::to_string", "[collisions]")
{
    const double dhat = 1e-1;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    const bool success =
        tests::load_mesh("two-cubes-close.ply", vertices, edges, faces);
    REQUIRE(success);

    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions collisions;
    collisions.set_use_area_weighting(true);
    collisions.set_collision_set_type(
        NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX);
    collisions.build(mesh, vertices, dhat);

    std::sort(collisions.vv_collisions.begin(), collisions.vv_collisions.end());
    collisions.vv_collisions.erase(
        collisions.vv_collisions.begin() + 1,
        collisions.vv_collisions.end() - 1);

    std::sort(collisions.ev_collisions.begin(), collisions.ev_collisions.end());
    collisions.ev_collisions.erase(
        collisions.ev_collisions.begin() + 1,
        collisions.ev_collisions.end() - 1);

    std::sort(collisions.ee_collisions.begin(), collisions.ee_collisions.end());
    collisions.ee_collisions.erase(
        collisions.ee_collisions.begin() + 1,
        collisions.ee_collisions.end() - 1);

    std::sort(collisions.fv_collisions.begin(), collisions.fv_collisions.end());
    collisions.fv_collisions.erase(
        collisions.fv_collisions.begin() + 1,
        collisions.fv_collisions.end() - 1);

    std::string s = collisions.to_string(mesh, vertices);

    CHECK(s == R"ipc_Qu8mg5v7(
vv: 4 140, w: -0.0128182, d: 0.00568928
vv: 119 243, w: -0.0263598, d: 0.00727215
ev: 14=(18, 37) 187, w: 0.0140559, d: 0.00500364
ev: 719=(219, 243) 119, w: 0.0122856, d: 0.00727213
ee: 14=(18, 37) 456=(144, 187), w: -0.00297852, dtype: 0, d: 0.00538125
ee: 346=(67, 119) 718=(218, 243), w: -0.00328539, dtype: 5, d: 0.00717647
fv: 17=(37, 18, 63) 187, w: 0.0137957, d: 0.00500269
fv: 471=(140, 220, 241) 55, w: 0.0160469, d: 0.006392)ipc_Qu8mg5v7");
}