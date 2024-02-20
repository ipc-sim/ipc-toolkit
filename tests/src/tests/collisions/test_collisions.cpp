#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/collisions/collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>

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

    CollisionMesh mesh(vertices, Eigen::MatrixXi(), Eigen::MatrixXi());
    mesh.init_area_jacobians();

    CHECK(mesh.num_vertices() == 8);
    CHECK(mesh.num_codim_vertices() == 8);
    CHECK(mesh.num_edges() == 0);
    CHECK(mesh.num_faces() == 0);

    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();
    CAPTURE(method);

    SECTION("Candidates")
    {
        Eigen::MatrixXd V1 = vertices;
        V1.col(1) *= 0.5;

        Candidates candidates;
        candidates.build(mesh, vertices, V1, thickness, method);

        CHECK(candidates.size() > 0);
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
        const bool use_convergent_formulation = GENERATE(false, true);
        const bool are_shape_derivatives_enabled = GENERATE(false, true);
        const double dhat = 0.25;

        Collisions collisions;
        collisions.set_use_convergent_formulation(use_convergent_formulation);
        collisions.set_are_shape_derivatives_enabled(
            are_shape_derivatives_enabled);

        collisions.build(mesh, vertices, dhat, min_distance, method);

        CHECK(collisions.size() == 12);
        CHECK(collisions.vv_collisions.size() == 12);

        BarrierPotential barrier_potential(dhat);

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

    CollisionMesh mesh(vertices, edges, Eigen::MatrixXi());
    mesh.init_area_jacobians();

    CHECK(mesh.num_vertices() == 8);
    CHECK(mesh.num_codim_vertices() == 3);
    CHECK(mesh.num_codim_edges() == 4);
    CHECK(mesh.num_edges() == 4);
    CHECK(mesh.num_faces() == 0);

    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();
    CAPTURE(method);

    SECTION("Candidates")
    {
        Eigen::MatrixXd V1 = vertices;
        V1.bottomRows(3).col(1).array() -= 4; // Translate the codim vertices

        Candidates candidates;
        candidates.build(mesh, vertices, V1, thickness, method);

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
        const bool use_convergent_formulation = GENERATE(false, true);
        const bool are_shape_derivatives_enabled = GENERATE(false, true);

        Collisions collisions;
        collisions.set_use_convergent_formulation(use_convergent_formulation);
        collisions.set_are_shape_derivatives_enabled(
            are_shape_derivatives_enabled);

        const double dhat = 0.25;
        collisions.build(mesh, vertices, dhat, /*min_distance=*/0.8, method);

        const int expected_num_collisions = 6 + int(use_convergent_formulation);
        const int expected_num_vv_collisions =
            2 + int(use_convergent_formulation);

        CHECK(collisions.size() == expected_num_collisions);
        CHECK(collisions.vv_collisions.size() == expected_num_vv_collisions);
        CHECK(collisions.ev_collisions.size() == 4);
        CHECK(collisions.ee_collisions.size() == 0);
        CHECK(collisions.fv_collisions.size() == 0);

        CHECK(BarrierPotential(dhat)(collisions, mesh, vertices) > 0.0);
    }
}

TEST_CASE("Vertex-Vertex Collision", "[collision][vertex-vertex]")
{
    CHECK(VertexVertexCollision(0, 1) == VertexVertexCollision(0, 1));
    CHECK(VertexVertexCollision(0, 1) == VertexVertexCollision(1, 0));
    CHECK(VertexVertexCollision(0, 1) != VertexVertexCollision(0, 2));
    CHECK(VertexVertexCollision(0, 1) != VertexVertexCollision(2, 0));
    CHECK(VertexVertexCollision(0, 1) < VertexVertexCollision(0, 2));
    CHECK(VertexVertexCollision(0, 1) < VertexVertexCollision(2, 0));

    CHECK(
        VertexVertexCollision(VertexVertexCandidate(0, 1))
        == VertexVertexCollision(0, 1));
}

TEST_CASE("Edge-Vertex Collision", "[collision][edge-vertex]")
{
    CHECK(EdgeVertexCollision(0, 1) == EdgeVertexCollision(0, 1));
    CHECK(EdgeVertexCollision(0, 1) != EdgeVertexCollision(1, 0));
    CHECK(EdgeVertexCollision(0, 1) != EdgeVertexCollision(0, 2));
    CHECK(EdgeVertexCollision(0, 1) != EdgeVertexCollision(2, 0));
    CHECK(EdgeVertexCollision(0, 1) < EdgeVertexCollision(0, 2));
    CHECK(!(EdgeVertexCollision(1, 1) < EdgeVertexCollision(0, 2)));
    CHECK(EdgeVertexCollision(0, 1) < EdgeVertexCollision(2, 0));

    CHECK(
        EdgeVertexCollision(EdgeVertexCandidate(0, 1))
        == EdgeVertexCollision(0, 1));
}

TEST_CASE("Edge-Edge Collision", "[collision][edge-edge]")
{
    CHECK(EdgeEdgeCollision(0, 1, 0.0) == EdgeEdgeCollision(0, 1, 1.0));
    CHECK(EdgeEdgeCollision(0, 1, 0.0) == EdgeEdgeCollision(1, 0, 1.0));
    CHECK(EdgeEdgeCollision(0, 1, 0.0) != EdgeEdgeCollision(0, 2, 1.0));
    CHECK(EdgeEdgeCollision(0, 1, 0.0) != EdgeEdgeCollision(2, 0, 1.0));
    CHECK(EdgeEdgeCollision(0, 1, 0.0) < EdgeEdgeCollision(0, 2, 1.0));
    CHECK(EdgeEdgeCollision(0, 1, 0.0) < EdgeEdgeCollision(2, 0, 1.0));

    CHECK(
        EdgeEdgeCollision(EdgeEdgeCandidate(0, 1), 0.0)
        == EdgeEdgeCollision(0, 1, 0.0));
}

TEST_CASE("Face-Vertex Collision", "[collision][face-vertex]")
{
    CHECK(FaceVertexCollision(0, 1) == FaceVertexCollision(0, 1));
    CHECK(FaceVertexCollision(0, 1) != FaceVertexCollision(1, 0));
    CHECK(FaceVertexCollision(0, 1) != FaceVertexCollision(0, 2));
    CHECK(FaceVertexCollision(0, 1) != FaceVertexCollision(2, 0));
    CHECK(FaceVertexCollision(0, 1) < FaceVertexCollision(0, 2));
    CHECK(!(FaceVertexCollision(1, 1) < FaceVertexCollision(0, 2)));
    CHECK(FaceVertexCollision(0, 1) < FaceVertexCollision(2, 0));

    CHECK(
        FaceVertexCollision(FaceVertexCandidate(0, 1))
        == FaceVertexCollision(0, 1));
}

TEST_CASE("Plane-Vertex Collision", "[collision][plane-vertex]")
{
    Eigen::MatrixXi edges, faces;
    const Eigen::Vector3d n(0, 1, 0), o(0, 0, 0);
    const PlaneVertexCollision c(o, n, 0);
    CHECK(c.num_vertices() == 1);
    CHECK(
        c.vertex_ids(edges, faces)
        == std::array<long, 4> { { 0, -1, -1, -1 } });
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

TEST_CASE("Collisions::is_*", "[collisions]")
{
    Collisions collisions;
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