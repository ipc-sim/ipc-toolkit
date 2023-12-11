#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE("Constraint Shape Derivative", "[constraint][shape_derivative]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.obj", V, E, F);

    const bool use_convergent_formulation = GENERATE(false);
    const double dhat = 1e-1;

    // Stack cube on top of itself
    E.conservativeResize(E.rows() * 2, E.cols());
    E.bottomRows(E.rows() / 2) = E.topRows(E.rows() / 2).array() + V.rows();

    F.conservativeResize(F.rows() * 2, F.cols());
    F.bottomRows(F.rows() / 2) = F.topRows(F.rows() / 2).array() + V.rows();

    V.conservativeResize(V.rows() * 2, V.cols());
    V.bottomRows(V.rows() / 2) = V.topRows(V.rows() / 2);
    V.bottomRows(V.rows() / 2).col(1).array() += 1 + 0.1 * dhat;

    // Rest positions
    Eigen::MatrixXd X = V;
    X.bottomRows(V.rows() / 2).col(1).array() += 1.0;

    // Displacements
    const Eigen::MatrixXd U = V - X;

    const int ndof = V.size();

    // ------------------------------------------------------------------------

    CollisionMesh mesh(X, E, F);
    mesh.init_area_jacobians();

    Candidates candidates;
    candidates.build(mesh, V, dhat);

    CollisionConstraints constraints;
    constraints.set_use_convergent_formulation(use_convergent_formulation);
    constraints.set_are_shape_derivatives_enabled(true);
    constraints.build(candidates, mesh, V, dhat);

    REQUIRE(constraints.ee_constraints.size() > 0);

    for (int i = 0; i < constraints.size(); i++) {
        if (use_convergent_formulation)
            break;

        std::vector<Eigen::Triplet<double>> triplets;
        constraints[i].compute_shape_derivative(X, V, E, F, dhat, triplets);
        Eigen::SparseMatrix<double> JF_wrt_X_sparse(ndof, ndof);
        JF_wrt_X_sparse.setFromTriplets(triplets.begin(), triplets.end());
        const Eigen::MatrixXd JF_wrt_X = Eigen::MatrixXd(JF_wrt_X_sparse);

        auto F_X = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            // TODO: Recompute weight based on x
            assert(use_convergent_formulation == false);
            // Recompute eps_x based on x
            double prev_eps_x;
            if (constraints.is_edge_edge(i)) {
                EdgeEdgeConstraint& c =
                    dynamic_cast<EdgeEdgeConstraint&>(constraints[i]);
                prev_eps_x = c.eps_x;
                c.eps_x = edge_edge_mollifier_threshold(
                    x.segment<3>(3 * E(c.edge0_id, 0)),
                    x.segment<3>(3 * E(c.edge0_id, 1)),
                    x.segment<3>(3 * E(c.edge1_id, 0)),
                    x.segment<3>(3 * E(c.edge1_id, 1)));
            }

            Eigen::VectorXd grad = Eigen::VectorXd::Zero(ndof);
            local_gradient_to_global_gradient(
                constraints[i].compute_potential_gradient(
                    fd::unflatten(x, X.cols()) + U, E, F, dhat),
                constraints[i].vertex_ids(E, F), V.cols(), grad);

            // Restore eps_x
            if (constraints.is_edge_edge(i)) {
                dynamic_cast<EdgeEdgeConstraint&>(constraints[i]).eps_x =
                    prev_eps_x;
            }

            return grad;
        };

        Eigen::MatrixXd fd_JF_wrt_X;
        fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);
        CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
        if (!fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X)) {
            tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
        }
    }

    // ------------------------------------------------------------------------

    const Eigen::MatrixXd JF_wrt_X =
        constraints.compute_shape_derivative(mesh, V, dhat);

    Eigen::MatrixXd sum = Eigen::MatrixXd::Zero(ndof, ndof);
    for (int i = 0; i < constraints.size(); i++) {
        std::vector<Eigen::Triplet<double>> triplets;
        constraints[i].compute_shape_derivative(X, V, E, F, dhat, triplets);
        Eigen::SparseMatrix<double> JF_wrt_X_sparse(ndof, ndof);
        JF_wrt_X_sparse.setFromTriplets(triplets.begin(), triplets.end());
        sum += Eigen::MatrixXd(JF_wrt_X_sparse);
    }
    CHECK(fd::compare_jacobian(JF_wrt_X, sum));

    auto F_X = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_X = fd::unflatten(x, X.cols());
        const Eigen::MatrixXd fd_V = fd_X + U;

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());

        // WARNING: This breaks the tests because EE distances are C0 when edges
        // are parallel
        // CollisionConstraints fd_constraints;
        // fd_constraints.set_use_convergent_formulation(
        //     constraints.use_convergent_formulation());
        // fd_constraints.build(fd_mesh, fd_V, dhat);

        return constraints.compute_potential_gradient(fd_mesh, fd_V, dhat);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);
    CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
    if (!fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X)) {
        tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }
}

TEST_CASE("Codim. Vertex-Vertex Constraints", "[constraints][codim]")
{
    constexpr double thickness = 0.4;
    constexpr double min_distance = 2 * thickness;

    Eigen::MatrixXd V(8, 3);
    V << 0, 0, 0, //
        0, 0, 1,  //
        0, 1, 0,  //
        0, 1, 1,  //
        1, 0, 0,  //
        1, 0, 1,  //
        1, 1, 0,  //
        1, 1, 1;
    V.rowwise() -= V.colwise().mean();

    CollisionMesh mesh(V, Eigen::MatrixXi(), Eigen::MatrixXi());
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
        Eigen::MatrixXd V1 = V;
        V1.col(1) *= 0.5;

        Candidates candidates;
        candidates.build(mesh, V, V1, thickness, method);

        CHECK(candidates.size() > 0);
        CHECK(candidates.vv_candidates.size() == candidates.size());

        CHECK(!candidates.is_step_collision_free(mesh, V, V1, min_distance));

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
                mesh, V, V1, min_distance)
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

        constraints.build(mesh, V, 0.25, min_distance, method);

        CHECK(constraints.size() == 12);
        CHECK(constraints.vv_constraints.size() == 12);

        CHECK(constraints.compute_potential(mesh, V, 0.25) > 0.0);
        const Eigen::VectorXd grad =
            constraints.compute_potential_gradient(mesh, V, 0.25);
        for (int i = 0; i < V.rows(); i++) {
            const Eigen::Vector3d f = -grad.segment<3>(3 * i);
            CHECK(f.normalized().isApprox(V.row(i).normalized().transpose()));
        }
    }
}

TEST_CASE("Codim. Edge-Vertex Constraints", "[constraints][codim]")
{
    constexpr double thickness = 1e-3;
    constexpr double min_distance = 2 * thickness;

    Eigen::MatrixXd V(8, 3);
    V << 0, 0, 0, //
        1, 0, 0,  //
        0, 0, -1, //
        -1, 0, 0, //
        0, 0, 1,  //
        0, 1, 0,  //
        0, 2, 0,  //
        0, 3, 0;
    Eigen::MatrixXi E(4, 2);
    E << 0, 1, //
        0, 2,  //
        0, 3,  //
        0, 4;

    CollisionMesh mesh(V, E, Eigen::MatrixXi());
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
        Eigen::MatrixXd V1 = V;
        V1.bottomRows(3).col(1).array() -= 4; // Translate the codim vertices

        Candidates candidates;
        candidates.build(mesh, V, V1, thickness, method);

        CHECK(candidates.size() == 15);
        CHECK(candidates.vv_candidates.size() == 3);
        CHECK(candidates.ev_candidates.size() == 12);
        CHECK(candidates.ee_candidates.size() == 0);
        CHECK(candidates.fv_candidates.size() == 0);

        CHECK(!candidates.is_step_collision_free(mesh, V, V1, min_distance));

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
                mesh, V, V1, min_distance)
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
        constraints.build(mesh, V, dhat, /*min_distance=*/0.8, method);

        const int expected_num_constraints =
            6 + int(use_convergent_formulation);
        const int expected_num_vv_constraints =
            2 + int(use_convergent_formulation);

        CHECK(constraints.size() == expected_num_constraints);
        CHECK(constraints.vv_constraints.size() == expected_num_vv_constraints);
        CHECK(constraints.ev_constraints.size() == 4);
        CHECK(constraints.ee_constraints.size() == 0);
        CHECK(constraints.fv_constraints.size() == 0);

        CHECK(constraints.compute_potential(mesh, V, dhat) > 0.0);
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
    Eigen::MatrixXi E, F;
    const Eigen::Vector3d n(0, 1, 0), o(0, 0, 0);
    const PlaneVertexConstraint c(o, n, 0);
    CHECK(c.num_vertices() == 1);
    CHECK(c.vertex_ids(E, F) == std::array<long, 4> { { 0, -1, -1, -1 } });
    CHECK(c.plane_origin == o);
    CHECK(c.plane_normal == n);
    CHECK(c.vertex_id == 0);

    CHECK(c.compute_distance(Eigen::RowVector3d(0, -2, 0), E, F) == 4.0);
    CHECK(c.compute_distance(Eigen::RowVector3d(0, 2, 0), E, F) == 4.0);
    CHECK(
        c.compute_distance_gradient(Eigen::RowVector3d(0, 2, 0), E, F)
        == Eigen::Vector3d(0, 4, 0));
    CHECK(
        c.compute_distance_hessian(Eigen::RowVector3d(0, 2, 0), E, F)
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