#include <catch2/catch_all.hpp>

#include "test_utils.hpp"

#include <ipc/ipc.hpp>
#include <ipc/config.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>

using namespace ipc;

TEST_CASE("Dummy test for IPC compilation", "[!benchmark][ipc]")
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

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = load_mesh(mesh_name, V, E, F);
    REQUIRE(success);

    const CollisionMesh mesh(V, E, F);

    CollisionConstraints collision_constraints;
    collision_constraints.set_use_convergent_formulation(
        use_convergent_formulation);
    collision_constraints.build(mesh, V, dhat);
    CAPTURE(mesh_name, dhat);
    CHECK(collision_constraints.size() > 0);

    BENCHMARK("Compute barrier potential")
    {
        return collision_constraints.compute_potential(mesh, V, dhat);
    };
    BENCHMARK("Compute barrier potential gradient")
    {
        return collision_constraints.compute_potential_gradient(mesh, V, dhat);
    };
    BENCHMARK("Compute barrier potential hessian")
    {
        return collision_constraints.compute_potential_hessian(mesh, V, dhat);
    };
    BENCHMARK("Compute barrier potential hessian with PSD projection")
    {
        return collision_constraints.compute_potential_hessian(
            mesh, V, dhat, true);
    };
    BENCHMARK("Compute compute_minimum_distance")
    {
        return collision_constraints.compute_minimum_distance(mesh, V);
    };
}

TEST_CASE("Test IPC full gradient", "[ipc][gradient]")
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
    // SECTION("bunny")
    // {
    //     dhat = 1e-2;
    //     mesh_name = "bunny.obj";
    // }

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = load_mesh(mesh_name, V, E, F);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;

    CollisionConstraints collision_constraints;
    collision_constraints.set_use_convergent_formulation(
        use_convergent_formulation);
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(V, E, F);
        collision_constraints.build(mesh, V, dhat, /*dmin=*/0, method);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(V, E, F);
        V = mesh.vertices(V);
        collision_constraints.build(mesh, V, dhat, /*dmin=*/0, method);
    }
    CAPTURE(dhat, method, all_vertices_on_surface);
    CHECK(collision_constraints.size() > 0);

    const Eigen::VectorXd grad_b =
        collision_constraints.compute_potential_gradient(mesh, V, dhat);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return collision_constraints.compute_potential(
            mesh, fd::unflatten(x, V.cols()), dhat);
    };
    Eigen::VectorXd fgrad_b;
    fd::finite_gradient(fd::flatten(V), f, fgrad_b);

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

TEST_CASE("Test IPC full hessian", "[ipc][hessian]")
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

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = load_mesh(mesh_name, V, E, F);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;

    CollisionConstraints collision_constraints;
    collision_constraints.set_use_convergent_formulation(
        use_convergent_formulation);
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(V, E, F);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(V, E, F);
        V = mesh.vertices(V);
    }
    collision_constraints.build(mesh, V, dhat, /*dmin=*/0, method);
    CAPTURE(dhat, method, all_vertices_on_surface);
    REQUIRE(collision_constraints.size() > 0);

    Eigen::MatrixXd hess_b =
        collision_constraints.compute_potential_hessian(mesh, V, dhat);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return collision_constraints.compute_potential_gradient(
            mesh, fd::unflatten(x, V.cols()), dhat);
    };
    Eigen::MatrixXd fhess_b;
    fd::finite_jacobian(fd::flatten(V), f, fhess_b);

    REQUIRE(hess_b.squaredNorm() > 1e-3);
    CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
}

TEST_CASE("Test IPC shape derivative", "[ipc][shape_opt]")
{
    nlohmann::json data;
    {
        std::ifstream input(TEST_DATA_DIR + "shape_derivative_data.json");
        REQUIRE(input.good());

        data = nlohmann::json::parse(input, nullptr, false);
        REQUIRE(!data.is_discarded());
    }

    // Parameters
    double dhat = data["dhat"];

    // Mesh
    Eigen::MatrixXd X, V;
    from_json(data["boundary_nodes_pos"], X);
    from_json(data["displaced"], V);

    Eigen::MatrixXi E;
    from_json(data["boundary_edges"], E);

    CollisionMesh mesh =
        CollisionMesh::build_from_full_mesh(X, E, /*faces=*/Eigen::MatrixXi());
    mesh.init_area_jacobians();
    REQUIRE(mesh.are_area_jacobians_initialized());

    X = mesh.vertices(X);
    V = mesh.vertices(V);
    const Eigen::MatrixXd U = V - X;

    CollisionConstraints collision_constraints;
    const bool use_convergent_formulation = GENERATE(true, false);
    collision_constraints.set_use_convergent_formulation(
        use_convergent_formulation);
    collision_constraints.set_are_shape_derivatives_enabled(true);
    collision_constraints.build(mesh, V, dhat);

    const Eigen::MatrixXd JF_wrt_X =
        collision_constraints.compute_shape_derivative(mesh, V, dhat);

    auto F_X = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_X = fd::unflatten(x, X.cols());
        const Eigen::MatrixXd fd_V = fd_X + U;

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());

        CollisionConstraints fd_constraint_set;
        fd_constraint_set.set_use_convergent_formulation(
            collision_constraints.use_convergent_formulation());
        fd_constraint_set.build(fd_mesh, fd_V, dhat);

        return fd_constraint_set.compute_potential_gradient(
            fd_mesh, fd_V, dhat);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);
    CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
    if (!fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X)) {
        print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }
}

TEST_CASE("Benchmark IPC shape derivative", "[ipc][shape_opt][!benchmark]")
{
    nlohmann::json data;
    {
        std::ifstream input(TEST_DATA_DIR + "shape_derivative_data.json");
        REQUIRE(input.good());

        data = nlohmann::json::parse(input, nullptr, false);
        REQUIRE(!data.is_discarded());
    }

    // Parameters
    double dhat = data["dhat"];

    // Mesh
    Eigen::MatrixXd X, V;
    from_json(data["boundary_nodes_pos"], X);
    from_json(data["displaced"], V);

    Eigen::MatrixXi E;
    from_json(data["boundary_edges"], E);

    CollisionMesh mesh =
        CollisionMesh::build_from_full_mesh(X, E, /*faces=*/Eigen::MatrixXi());
    mesh.init_area_jacobians();
    REQUIRE(mesh.are_area_jacobians_initialized());

    X = mesh.vertices(X);
    V = mesh.vertices(V);
    const Eigen::MatrixXd U = V - X;

    CollisionConstraints collision_constraints;
    const bool use_convergent_formulation = GENERATE(true, false);
    collision_constraints.set_use_convergent_formulation(
        use_convergent_formulation);
    collision_constraints.set_are_shape_derivatives_enabled(true);
    collision_constraints.build(mesh, V, dhat);

    Eigen::SparseMatrix<double> JF_wrt_X;

    BENCHMARK("Shape Derivative")
    {
        JF_wrt_X =
            collision_constraints.compute_shape_derivative(mesh, V, dhat);
    };
}

TEST_CASE("Test convergent formulation", "[ipc][convergent]")
{
    const bool use_convergent_formulation = GENERATE(false, true);
    const double dhat = 1e-3;

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    SECTION("2D Edge-Vertex")
    {
        //        .
        // .-------.-------.
        V.resize(4, 2);
        V.row(0) << 0, 1e-4;
        V.row(1) << -1, 0;
        V.row(2) << 1e-4, 0;
        V.row(3) << 1, 0;

        E.resize(2, 2);
        E.row(0) << 1, 2;
        E.row(1) << 2, 3;

        CHECK(point_point_distance(V.row(0), V.row(2)) < dhat * dhat);
    }
    SECTION("3D Face-Vertex")
    {
        V.resize(5, 3);
        V.row(0) << 0, 1e-4, 0;
        V.row(1) << -1, 0, 0;
        V.row(2) << 1e-4, 0, -1;
        V.row(3) << 1e-4, 0, 1;
        V.row(4) << 1, 0, 0;

        F.resize(2, 3);
        F.row(0) << 1, 2, 3;
        F.row(1) << 2, 3, 4;

        igl::edges(F, E);

        CHECK(point_edge_distance(V.row(0), V.row(2), V.row(3)) < dhat * dhat);
    }
    // SECTION("3D Edge-Edge")
    // {
    //     V.resize(5, 3);
    //     //
    //     V.row(0) << 0, 1e-4, -1;
    //     V.row(1) << 0, 1e-4, 0.9;
    //     //
    //     V.row(2) << 1e-4, 0, 0;
    //     V.row(3) << -0.33, 0, 0;
    //     V.row(4) << 0.5, 0, 0;

    //     E.resize(3, 2);
    //     E.row(0) << 0, 1;
    //     E.row(1) << 3, 2;
    //     E.row(2) << 2, 4;

    //     CHECK(point_edge_distance(V.row(2), V.row(0), V.row(1)) < dhat *
    //     dhat);
    // }
    // SECTION("3D Edge-Edge 2")
    // {
    //     V.resize(5, 3);
    //     //
    //     V.row(0) << 0, 1e-4, -1e-4;
    //     V.row(1) << 0, 1e-4, -1;
    //     //
    //     V.row(2) << 1e-4, 0, 0;
    //     V.row(3) << -0.33, 0, 0;
    //     V.row(4) << 0.5, 0, 0;

    //     E.resize(3, 2);
    //     E.row(0) << 0, 1;
    //     E.row(1) << 3, 2;
    //     E.row(2) << 2, 4;

    //     CHECK(point_edge_distance(V.row(2), V.row(0), V.row(1)) < dhat *
    //     dhat);
    // }

    const CollisionMesh mesh(V, E, F);

    CollisionConstraints collision_constraints;
    collision_constraints.set_use_convergent_formulation(
        use_convergent_formulation);

    collision_constraints.build(mesh, V, dhat);
    CHECK(collision_constraints.size() > 0);

    const Eigen::VectorXd grad_b =
        collision_constraints.compute_potential_gradient(mesh, V, dhat);

    const Eigen::MatrixXd force = -fd::unflatten(grad_b, V.cols());
    std::cout << "force:\n" << force << std::endl;

    if (use_convergent_formulation) {
        constexpr double eps = std::numeric_limits<double>::epsilon();
        CHECK(grad_b(0) == Catch::Approx(0).margin(eps));
        CHECK(grad_b(2 * V.cols()) == Catch::Approx(0).margin(eps));
        // CHECK(grad_b(3 * V.cols()) == 0);
    } else {
        CHECK(grad_b(0) != 0);
        CHECK(grad_b(2 * V.cols()) != 0);
        // CHECK(grad_b(3 * V.cols()) != 0);
    }

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return collision_constraints.compute_potential(
            mesh, fd::unflatten(x, V.cols()), dhat);
    };
    Eigen::VectorXd fgrad_b;
    fd::finite_gradient(fd::flatten(V), f, fgrad_b);

    CHECK(fd::compare_gradient(grad_b, fgrad_b));
}