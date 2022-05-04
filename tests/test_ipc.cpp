#include <catch2/catch.hpp>

#include <iostream>

#include <finitediff.hpp>

#include <ipc/ipc.hpp>

#include "test_utils.hpp"

using namespace ipc;

TEST_CASE("Dummy test for IPC compilation", "[!benchmark][ipc]")
{
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

    CollisionMesh mesh(V, E, F);

    Constraints constraint_set;
    construct_constraint_set(mesh, V, dhat, constraint_set);
    CAPTURE(mesh_name, dhat);
    CHECK(constraint_set.size() > 0);

    BENCHMARK("Compute barrier potential")
    {
        double b =
            ipc::compute_barrier_potential(mesh, V, constraint_set, dhat);
    };
    BENCHMARK("Compute barrier potential gradient")
    {
        Eigen::VectorXd grad_b = ipc::compute_barrier_potential_gradient(
            mesh, V, constraint_set, dhat);
    };
    BENCHMARK("Compute barrier potential hessian")
    {
        Eigen::MatrixXd hess_b = ipc::compute_barrier_potential_hessian(
            mesh, V, constraint_set, dhat);
    };
    BENCHMARK("Compute compute_minimum_distance")
    {
        double min_dist =
            ipc::compute_minimum_distance(mesh, V, constraint_set);
    };
}

TEST_CASE("Test IPC full gradient", "[ipc][gradient]")
{
    double dhat = -1;
    std::string mesh_name = "";
    bool all_vertices_on_surface = true;

    BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

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

    Constraints constraint_set;
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(V, E, F);
        ipc::construct_constraint_set(
            mesh, V, dhat, constraint_set, /*dmin=*/0, method);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(V, E, F);
        V = mesh.vertices(V);
        ipc::construct_constraint_set(
            mesh, V, dhat, constraint_set, /*dmin=*/0, method);
    }
    CAPTURE(dhat, method, all_vertices_on_surface);
    CHECK(constraint_set.size() > 0);

    Eigen::VectorXd grad_b =
        ipc::compute_barrier_potential_gradient(mesh, V, constraint_set, dhat);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return ipc::compute_barrier_potential(
            mesh, fd::unflatten(x, V.cols()), constraint_set, dhat);
    };
    Eigen::VectorXd fgrad_b;
    fd::finite_gradient(fd::flatten(V), f, fgrad_b);

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

TEST_CASE("Test IPC full hessian", "[ipc][hessian]")
{
    double dhat = -1;
    std::string mesh_name = "";
    bool all_vertices_on_surface = true;

    BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

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

    Constraints constraint_set;
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(V, E, F);
        ipc::construct_constraint_set(
            mesh, V, dhat, constraint_set, /*dmin=*/0, method);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(V, E, F);
        V = mesh.vertices(V);
        ipc::construct_constraint_set(
            mesh, V, dhat, constraint_set, /*dmin=*/0, method);
    }
    CAPTURE(dhat, method, all_vertices_on_surface);
    CHECK(constraint_set.size() > 0);

    Eigen::MatrixXd hess_b = ipc::compute_barrier_potential_hessian(
        mesh, V, constraint_set, dhat, /*project_to_psd=*/false);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return ipc::compute_barrier_potential_gradient(
            mesh, fd::unflatten(x, V.cols()), constraint_set, dhat);
    };
    Eigen::MatrixXd fhess_b;
    fd::finite_jacobian(fd::flatten(V), f, fhess_b);

    REQUIRE(hess_b.squaredNorm() > 1e-3);
    CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
}
