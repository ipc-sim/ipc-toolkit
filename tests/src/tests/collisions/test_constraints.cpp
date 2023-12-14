#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/potentials/barrier_potential.hpp>

using namespace ipc;

TEST_CASE("Codim. vertex-vertex constraints", "[constraints][codim]")
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

    // These methods do not support vertex-vertex candidates
    if (method == ipc::BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE
        || method == ipc::BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU) {
        return;
    }

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
#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
        constexpr double conservative_min_dist = 1e-4;
#else
        constexpr double conservative_min_dist = 0.2 * (1 - min_distance);
#endif
        constexpr double expected_toi =
            (1 - (min_distance + conservative_min_dist)) / 2.0 / 0.25;
        CHECK(
            candidates.compute_collision_free_stepsize(
                mesh, vertices, V1, min_distance)
            == Catch::Approx(expected_toi));
    }

    SECTION("Constraints")
    {
        const bool use_convergent_formulation = GENERATE(false, true);
        const bool are_shape_derivatives_enabled = GENERATE(false, true);
        const double dhat = 0.25;

        CollisionConstraints constraints;
        constraints.set_use_convergent_formulation(use_convergent_formulation);
        constraints.set_are_shape_derivatives_enabled(
            are_shape_derivatives_enabled);

        constraints.build(mesh, vertices, dhat, min_distance, method);

        CHECK(constraints.size() == 12);
        CHECK(constraints.vv_constraints.size() == 12);

        BarrierPotential barrier_potential(dhat);

        CHECK(barrier_potential(mesh, vertices, constraints) > 0.0);
        const Eigen::VectorXd grad =
            barrier_potential.gradient(mesh, vertices, constraints);
        for (int i = 0; i < vertices.rows(); i++) {
            const Eigen::Vector3d f = -grad.segment<3>(3 * i);
            CHECK(f.normalized().isApprox(
                vertices.row(i).normalized().transpose()));
        }
    }
}

TEST_CASE("Codim. edge-vertex constraints", "[constraints][codim]")
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

    // These methods do not support vertex-vertex candidates
    if (method == ipc::BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE
        || method == ipc::BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU) {
        return;
    }

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
#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
        constexpr double conservative_min_dist = 1e-4;
#else
        constexpr double conservative_min_dist = 0.2 * (1 - min_distance);
#endif
        constexpr double expected_toi =
            (1 - (min_distance + conservative_min_dist)) / 4;
        CHECK(
            candidates.compute_collision_free_stepsize(
                mesh, vertices, V1, min_distance)
            == Catch::Approx(expected_toi));
    }

    SECTION("Constraints")
    {
        const bool use_convergent_formulation = GENERATE(false, true);
        const bool are_shape_derivatives_enabled = GENERATE(false, true);

        CollisionConstraints constraints;
        constraints.set_use_convergent_formulation(use_convergent_formulation);
        constraints.set_are_shape_derivatives_enabled(
            are_shape_derivatives_enabled);

        const double dhat = 0.25;
        constraints.build(mesh, vertices, dhat, /*min_distance=*/0.8, method);

        const int expected_num_constraints =
            6 + int(use_convergent_formulation);
        const int expected_num_vv_constraints =
            2 + int(use_convergent_formulation);

        CHECK(constraints.size() == expected_num_constraints);
        CHECK(constraints.vv_constraints.size() == expected_num_vv_constraints);
        CHECK(constraints.ev_constraints.size() == 4);
        CHECK(constraints.ee_constraints.size() == 0);
        CHECK(constraints.fv_constraints.size() == 0);

        CHECK(BarrierPotential(dhat)(mesh, vertices, constraints) > 0.0);
    }
}

TEST_CASE("Vertex-Vertex Constraint", "[constraint][vertex-vertex]")
{
    CHECK(VertexVertexConstraint(0, 1) == VertexVertexConstraint(0, 1));
    CHECK(VertexVertexConstraint(0, 1) == VertexVertexConstraint(1, 0));
    CHECK(VertexVertexConstraint(0, 1) != VertexVertexConstraint(0, 2));
    CHECK(VertexVertexConstraint(0, 1) != VertexVertexConstraint(2, 0));
    CHECK(VertexVertexConstraint(0, 1) < VertexVertexConstraint(0, 2));
    CHECK(VertexVertexConstraint(0, 1) < VertexVertexConstraint(2, 0));

    CHECK(
        VertexVertexConstraint(VertexVertexCandidate(0, 1))
        == VertexVertexConstraint(0, 1));
}

TEST_CASE("Edge-Vertex Constraint", "[constraint][edge-vertex]")
{
    CHECK(EdgeVertexConstraint(0, 1) == EdgeVertexConstraint(0, 1));
    CHECK(EdgeVertexConstraint(0, 1) != EdgeVertexConstraint(1, 0));
    CHECK(EdgeVertexConstraint(0, 1) != EdgeVertexConstraint(0, 2));
    CHECK(EdgeVertexConstraint(0, 1) != EdgeVertexConstraint(2, 0));
    CHECK(EdgeVertexConstraint(0, 1) < EdgeVertexConstraint(0, 2));
    CHECK(!(EdgeVertexConstraint(1, 1) < EdgeVertexConstraint(0, 2)));
    CHECK(EdgeVertexConstraint(0, 1) < EdgeVertexConstraint(2, 0));

    CHECK(
        EdgeVertexConstraint(EdgeVertexCandidate(0, 1))
        == EdgeVertexConstraint(0, 1));
}

TEST_CASE("Edge-Edge Constraint", "[constraint][edge-edge]")
{
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) == EdgeEdgeConstraint(0, 1, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) == EdgeEdgeConstraint(1, 0, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) != EdgeEdgeConstraint(0, 2, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) != EdgeEdgeConstraint(2, 0, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) < EdgeEdgeConstraint(0, 2, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) < EdgeEdgeConstraint(2, 0, 1.0));

    CHECK(
        EdgeEdgeConstraint(EdgeEdgeCandidate(0, 1), 0.0)
        == EdgeEdgeConstraint(0, 1, 0.0));
}

TEST_CASE("Face-Vertex Constraint", "[constraint][face-vertex]")
{
    CHECK(FaceVertexConstraint(0, 1) == FaceVertexConstraint(0, 1));
    CHECK(FaceVertexConstraint(0, 1) != FaceVertexConstraint(1, 0));
    CHECK(FaceVertexConstraint(0, 1) != FaceVertexConstraint(0, 2));
    CHECK(FaceVertexConstraint(0, 1) != FaceVertexConstraint(2, 0));
    CHECK(FaceVertexConstraint(0, 1) < FaceVertexConstraint(0, 2));
    CHECK(!(FaceVertexConstraint(1, 1) < FaceVertexConstraint(0, 2)));
    CHECK(FaceVertexConstraint(0, 1) < FaceVertexConstraint(2, 0));

    CHECK(
        FaceVertexConstraint(FaceVertexCandidate(0, 1))
        == FaceVertexConstraint(0, 1));
}

TEST_CASE("Plane-Vertex Constraint", "[constraint][plane-vertex]")
{
    Eigen::MatrixXi edges, faces;
    const Eigen::Vector3d n(0, 1, 0), o(0, 0, 0);
    const PlaneVertexConstraint c(o, n, 0);
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

TEST_CASE("is_*", "[constraints]")
{
    CollisionConstraints constraints;
    constraints.vv_constraints.emplace_back(0, 1);
    constraints.ev_constraints.emplace_back(0, 1);
    constraints.ee_constraints.emplace_back(0, 1, 0.0);
    constraints.fv_constraints.emplace_back(0, 1);
    constraints.pv_constraints.emplace_back(
        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 0), 0);

    for (int i = 0; i < constraints.size(); i++) {
        CHECK(constraints.is_vertex_vertex(i) == (i == 0));
        CHECK(constraints.is_edge_vertex(i) == (i == 1));
        CHECK(constraints.is_edge_edge(i) == (i == 2));
        CHECK(constraints.is_face_vertex(i) == (i == 3));
        CHECK(constraints.is_plane_vertex(i) == (i == 4));
    }
}