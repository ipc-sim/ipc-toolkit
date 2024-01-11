#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>

using namespace ipc;

TEST_CASE(
    "Barrier potential full gradient and hessian",
    "[potential][barrier_potential][gradient][hessian]")
{
    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();
    const bool use_convergent_formulation = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name = "";
    bool all_vertices_on_surface = true;
    SECTION("cube")
    {
        dhat = sqrt(2.0);
        mesh_name = "cube.obj";
    }
    SECTION("two cubes far")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-far.obj";
        all_vertices_on_surface = false;
    }
    SECTION("two cubes close")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-close.obj";
        all_vertices_on_surface = false;
    }
    // WARNING: The bunny takes too long in debug.
    // SECTION("bunny")
    // {
    //     dhat = 1e-2;
    //     mesh_name = "bunny.obj";
    // }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;

    Collisions collisions;
    collisions.set_use_convergent_formulation(use_convergent_formulation);
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(vertices, edges, faces);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(vertices, edges, faces);
        vertices = mesh.vertices(vertices);
    }
    collisions.build(mesh, vertices, dhat, /*dmin=*/0, method);
    CAPTURE(dhat, method, all_vertices_on_surface);
    CHECK(collisions.size() > 0);

    BarrierPotential barrier_potential(dhat);

    // -------------------------------------------------------------------------
    // Gradient
    // -------------------------------------------------------------------------

    const Eigen::VectorXd grad_b =
        barrier_potential.gradient(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return barrier_potential(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(fd::flatten(vertices), f, fgrad_b);
    }

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    CHECK(fd::compare_gradient(grad_b, fgrad_b));

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    Eigen::MatrixXd hess_b =
        barrier_potential.hessian(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return barrier_potential.gradient(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_jacobian(fd::flatten(vertices), f, fhess_b);
    }

    REQUIRE(hess_b.squaredNorm() > 1e-3);
    CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
}

TEST_CASE(
    "Barrier potential convergent formulation",
    "[potential][barrier_potential][convergent]")
{
    const bool use_convergent_formulation = GENERATE(false, true);
    const double dhat = 1e-3;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    SECTION("2D Edge-Vertex")
    {
        //        .
        // .-------.-------.
        vertices.resize(4, 2);
        vertices.row(0) << 0, 1e-4;
        vertices.row(1) << -1, 0;
        vertices.row(2) << 1e-4, 0;
        vertices.row(3) << 1, 0;

        edges.resize(2, 2);
        edges.row(0) << 1, 2;
        edges.row(1) << 2, 3;

        CHECK(
            point_point_distance(vertices.row(0), vertices.row(2))
            < dhat * dhat);
    }
    SECTION("3D Face-Vertex")
    {
        vertices.resize(5, 3);
        vertices.row(0) << 0, 1e-4, 0;
        vertices.row(1) << -1, 0, 0;
        vertices.row(2) << 1e-4, 0, -1;
        vertices.row(3) << 1e-4, 0, 1;
        vertices.row(4) << 1, 0, 0;

        faces.resize(2, 3);
        faces.row(0) << 1, 2, 3;
        faces.row(1) << 2, 3, 4;

        igl::edges(faces, edges);

        CHECK(
            point_edge_distance(
                vertices.row(0), vertices.row(2), vertices.row(3))
            < dhat * dhat);
    }
    SECTION("3D Edge-Edge")
    {
        vertices.resize(5, 3);
        //
        vertices.row(0) << 0, 1e-4, -1;
        vertices.row(1) << 0, 1e-4, 1;
        //
        vertices.row(2) << -1e-4, 0, 0;
        vertices.row(3) << -1, 0, 0;
        vertices.row(4) << 1, 0, 0;

        edges.resize(3, 2);
        edges.row(0) << 0, 1;
        edges.row(1) << 3, 2;
        edges.row(2) << 2, 4;

        CHECK(
            point_edge_distance(
                vertices.row(2), vertices.row(0), vertices.row(1))
            < dhat * dhat);
    }
    SECTION("3D Edge-Edge Parallel")
    {
        vertices.resize(5, 3);
        //
        vertices.row(0) << -0.5, 1e-5, -1e-3;
        vertices.row(1) << 0.5, 1e-5, 1e-3;
        //
        vertices.row(2) << -1, -1e-5, 0;
        vertices.row(3) << 0, -1e-5, 0;
        vertices.row(4) << 1, -1e-5, 0;

        edges.resize(3, 2);
        edges.row(0) << 0, 1;
        edges.row(1) << 2, 3;
        edges.row(2) << 3, 4;

        CHECK(
            point_edge_distance(
                vertices.row(3), vertices.row(0), vertices.row(1))
            < dhat * dhat);
    }

    const CollisionMesh mesh(vertices, edges, faces);

    Collisions collisions;
    collisions.set_use_convergent_formulation(use_convergent_formulation);

    collisions.build(mesh, vertices, dhat);
    CHECK(collisions.size() > 0);

    BarrierPotential barrier_potential(dhat);

    const Eigen::VectorXd grad_b =
        barrier_potential.gradient(collisions, mesh, vertices);

    // const Eigen::MatrixXd force = -fd::unflatten(grad_b, vertices.cols());
    // std::cout << "force:\n" << force << std::endl;

    // if (use_convergent_formulation) {
    //     constexpr double eps = std::numeric_limits<double>::epsilon();
    //     CHECK(grad_b(0) == Catch::Approx(0).margin(eps));
    //     CHECK(grad_b(2 * vertices.cols()) == Catch::Approx(0).margin(eps));
    // } else {
    //     CHECK(grad_b(0) != 0);
    //     CHECK(grad_b(2 * vertices.cols()) != 0);
    // }

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_V = fd::unflatten(x, mesh.dim());

        Collisions fd_collisions;
        fd_collisions.set_use_convergent_formulation(
            use_convergent_formulation);

        fd_collisions.build(mesh, fd_V, dhat);

        return barrier_potential(collisions, mesh, fd_V);
    };
    Eigen::VectorXd fgrad_b;
    fd::finite_gradient(fd::flatten(vertices), f, fgrad_b);

    CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

TEST_CASE(
    "Barrier potential shape derivative",
    "[potential][barrier_potential][shape_derivative]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    tests::load_mesh("cube.obj", vertices, edges, faces);

    const bool use_convergent_formulation = GENERATE(false);
    const double dhat = 1e-1;

    // Stack cube on top of itself
    edges.conservativeResize(edges.rows() * 2, edges.cols());
    edges.bottomRows(edges.rows() / 2) =
        edges.topRows(edges.rows() / 2).array() + vertices.rows();

    faces.conservativeResize(faces.rows() * 2, faces.cols());
    faces.bottomRows(faces.rows() / 2) =
        faces.topRows(faces.rows() / 2).array() + vertices.rows();

    vertices.conservativeResize(vertices.rows() * 2, vertices.cols());
    vertices.bottomRows(vertices.rows() / 2) =
        vertices.topRows(vertices.rows() / 2);
    vertices.bottomRows(vertices.rows() / 2).col(1).array() += 1 + 0.1 * dhat;

    // Rest positions
    Eigen::MatrixXd rest_positions = vertices;
    rest_positions.bottomRows(vertices.rows() / 2).col(1).array() += 1.0;

    // Displacements
    const Eigen::MatrixXd displacements = vertices - rest_positions;

    const int ndof = vertices.size();

    // ------------------------------------------------------------------------

    CollisionMesh mesh(rest_positions, edges, faces);
    mesh.init_area_jacobians();

    Candidates candidates;
    candidates.build(mesh, vertices, dhat);

    Collisions collisions;
    collisions.set_use_convergent_formulation(use_convergent_formulation);
    collisions.set_are_shape_derivatives_enabled(true);
    collisions.build(candidates, mesh, vertices, dhat);
    REQUIRE(collisions.ee_collisions.size() > 0);

    BarrierPotential barrier_potential(dhat);

    for (int i = 0; i < collisions.size(); i++) {
        if (use_convergent_formulation)
            break;

        std::vector<Eigen::Triplet<double>> triplets;
        barrier_potential.shape_derivative(
            collisions[i], collisions[i].vertex_ids(edges, faces),
            collisions[i].dof(rest_positions, edges, faces),
            collisions[i].dof(vertices, edges, faces), triplets);
        Eigen::SparseMatrix<double> JF_wrt_X_sparse(ndof, ndof);
        JF_wrt_X_sparse.setFromTriplets(triplets.begin(), triplets.end());
        const Eigen::MatrixXd JF_wrt_X = Eigen::MatrixXd(JF_wrt_X_sparse);

        auto F_X = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            // TODO: Recompute weight based on x
            assert(use_convergent_formulation == false);
            // Recompute eps_x based on x
            double prev_eps_x = -1;
            if (collisions.is_edge_edge(i)) {
                EdgeEdgeCollision& c =
                    dynamic_cast<EdgeEdgeCollision&>(collisions[i]);
                prev_eps_x = c.eps_x;
                c.eps_x = edge_edge_mollifier_threshold(
                    x.segment<3>(3 * edges(c.edge0_id, 0)),
                    x.segment<3>(3 * edges(c.edge0_id, 1)),
                    x.segment<3>(3 * edges(c.edge1_id, 0)),
                    x.segment<3>(3 * edges(c.edge1_id, 1)));
            }

            const Eigen::MatrixXd positions =
                fd::unflatten(x, rest_positions.cols()) + displacements;
            const VectorMax12d dof = collisions[i].dof(positions, edges, faces);

            Eigen::VectorXd grad = Eigen::VectorXd::Zero(ndof);
            local_gradient_to_global_gradient(
                barrier_potential.gradient(collisions[i], dof),
                collisions[i].vertex_ids(edges, faces), vertices.cols(), grad);

            // Restore eps_x
            if (collisions.is_edge_edge(i)) {
                REQUIRE(prev_eps_x >= 0);
                dynamic_cast<EdgeEdgeCollision&>(collisions[i]).eps_x =
                    prev_eps_x;
            }

            return grad;
        };

        Eigen::MatrixXd fd_JF_wrt_X;
        fd::finite_jacobian(fd::flatten(rest_positions), F_X, fd_JF_wrt_X);

        CHECKED_ELSE(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X))
        {
            tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
        }
    }

    // ------------------------------------------------------------------------

    const Eigen::MatrixXd JF_wrt_X =
        barrier_potential.shape_derivative(collisions, mesh, vertices);

    Eigen::MatrixXd sum = Eigen::MatrixXd::Zero(ndof, ndof);
    for (int i = 0; i < collisions.size(); i++) {
        std::vector<Eigen::Triplet<double>> triplets;
        barrier_potential.shape_derivative(
            collisions[i], collisions[i].vertex_ids(edges, faces),
            collisions[i].dof(rest_positions, edges, faces),
            collisions[i].dof(vertices, edges, faces), triplets);
        Eigen::SparseMatrix<double> JF_wrt_X_sparse(ndof, ndof);
        JF_wrt_X_sparse.setFromTriplets(triplets.begin(), triplets.end());
        sum += Eigen::MatrixXd(JF_wrt_X_sparse);
    }
    CHECK(fd::compare_jacobian(JF_wrt_X, sum));

    auto F_X = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_X = fd::unflatten(x, rest_positions.cols());
        const Eigen::MatrixXd fd_V = fd_X + displacements;

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());

        // WARNING: This breaks the tests because EE distances are C0 when edges
        // are parallel
        // Collisions fd_collisions;
        // fd_collisions.set_use_convergent_formulation(
        //     collisions.use_convergent_formulation());
        // fd_collisions.build(fd_mesh, fd_V, dhat);

        return barrier_potential.gradient(collisions, fd_mesh, fd_V);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(rest_positions), F_X, fd_JF_wrt_X);

    CHECKED_ELSE(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X))
    {
        tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }
}

TEST_CASE(
    "Barrier potential shape derivative (sim data)",
    "[potential][barrier_potential][shape_derivative]")
{
    nlohmann::json data;
    {
        std::ifstream input(tests::DATA_DIR / "shape_derivative_data.json");
        REQUIRE(input.good());

        data = nlohmann::json::parse(input, nullptr, false);
        REQUIRE(!data.is_discarded());
    }

    // Parameters
    double dhat = data["dhat"];

    // Mesh
    const Eigen::MatrixXi edges = data["boundary_edges"];
    Eigen::MatrixXd X = data["boundary_nodes_pos"];
    Eigen::MatrixXd vertices = data["displaced"];

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(
        X, edges, /*faces=*/Eigen::MatrixXi());
    mesh.init_area_jacobians();
    REQUIRE(mesh.are_area_jacobians_initialized());

    X = mesh.vertices(X);
    vertices = mesh.vertices(vertices);
    const Eigen::MatrixXd U = vertices - X;

    Collisions collisions;
    const bool use_convergent_formulation = GENERATE(true, false);
    collisions.set_use_convergent_formulation(use_convergent_formulation);
    collisions.set_are_shape_derivatives_enabled(true);
    collisions.build(mesh, vertices, dhat);

    BarrierPotential barrier_potential(dhat);

    const Eigen::MatrixXd JF_wrt_X =
        barrier_potential.shape_derivative(collisions, mesh, vertices);

    auto F_X = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_X = fd::unflatten(x, X.cols());
        const Eigen::MatrixXd fd_V = fd_X + U;

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());

        Collisions fd_collision_set;
        fd_collision_set.set_use_convergent_formulation(
            collisions.use_convergent_formulation());
        fd_collision_set.build(fd_mesh, fd_V, dhat);

        return barrier_potential.gradient(fd_collision_set, fd_mesh, fd_V);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);

    CHECKED_ELSE(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X))
    {
        tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }
}

// -- Benchmarking ------------------------------------------------------------

TEST_CASE(
    "Benchmark barrier potential", "[!benchmark][potential][barrier_potential]")
{
    const bool use_convergent_formulation = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name;
    SECTION("cube")
    {
        dhat = sqrt(2.0);
        mesh_name = "cube.obj";
    }
    SECTION("bunny")
    {
        dhat = 1e-2;
        mesh_name = "bunny.obj";
    }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    REQUIRE(success);

    const CollisionMesh mesh(vertices, edges, faces);

    Collisions collisions;
    collisions.set_use_convergent_formulation(use_convergent_formulation);
    collisions.build(mesh, vertices, dhat);
    CAPTURE(mesh_name, dhat);
    CHECK(collisions.size() > 0);

    BarrierPotential barrier_potential(dhat);

    BENCHMARK("Compute barrier potential")
    {
        return barrier_potential(collisions, mesh, vertices);
    };
    BENCHMARK("Compute barrier potential gradient")
    {
        return barrier_potential.gradient(collisions, mesh, vertices);
    };
    BENCHMARK("Compute barrier potential hessian")
    {
        return barrier_potential.hessian(collisions, mesh, vertices);
    };
    BENCHMARK("Compute barrier potential hessian with PSD projection")
    {
        return barrier_potential.hessian(collisions, mesh, vertices, true);
    };
    BENCHMARK("Compute compute_minimum_distance")
    {
        return collisions.compute_minimum_distance(mesh, vertices);
    };
}

TEST_CASE(
    "Benchmark barrier potential shape derivative",
    "[!benchmark][potential][barrier_potential][shape_derivative]")
{
    nlohmann::json data;
    {
        std::ifstream input(tests::DATA_DIR / "shape_derivative_data.json");
        REQUIRE(input.good());

        data = nlohmann::json::parse(input, nullptr, false);
        REQUIRE(!data.is_discarded());
    }

    // Parameters
    double dhat = data["dhat"];

    // Mesh
    const Eigen::MatrixXi edges = data["boundary_edges"];
    Eigen::MatrixXd X = data["boundary_nodes_pos"];
    Eigen::MatrixXd vertices = data["displaced"];

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(
        X, edges, /*faces=*/Eigen::MatrixXi());
    mesh.init_area_jacobians();
    REQUIRE(mesh.are_area_jacobians_initialized());

    X = mesh.vertices(X);
    vertices = mesh.vertices(vertices);
    const Eigen::MatrixXd U = vertices - X;

    Collisions collisions;
    const bool use_convergent_formulation = GENERATE(true, false);
    collisions.set_use_convergent_formulation(use_convergent_formulation);
    collisions.set_are_shape_derivatives_enabled(true);
    collisions.build(mesh, vertices, dhat);

    BarrierPotential barrier_potential(dhat);

    Eigen::SparseMatrix<double> JF_wrt_X;

    BENCHMARK("Shape Derivative")
    {
        JF_wrt_X =
            barrier_potential.shape_derivative(collisions, mesh, vertices);
    };
}
