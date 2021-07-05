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

    Constraints constraint_set;
    construct_constraint_set(/*V_rest=*/V, V, E, F, dhat, constraint_set);
    CAPTURE(mesh_name, dhat);
    CHECK(constraint_set.num_constraints() > 0);

    BENCHMARK("Compute barrier potential")
    {
        double b =
            ipc::compute_barrier_potential(V, E, F, constraint_set, dhat);
    };
    BENCHMARK("Compute barrier potential gradient")
    {
        Eigen::VectorXd grad_b = ipc::compute_barrier_potential_gradient(
            V, E, F, constraint_set, dhat);
    };
    BENCHMARK("Compute barrier potential hessian")
    {
        Eigen::MatrixXd hess_b = ipc::compute_barrier_potential_hessian(
            V, E, F, constraint_set, dhat);
    };
    BENCHMARK("Compute compute_minimum_distance")
    {
        double min_dist =
            ipc::compute_minimum_distance(V, E, F, constraint_set);
    };
}

TEST_CASE("Test IPC full gradient", "[ipc][gradient]")
{
    double dhat = -1;
    std::string mesh_name;

    SECTION("cube")
    {
        dhat = sqrt(2.0);
        mesh_name = "cube.obj";
    }
    SECTION("two cubes far")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-far.obj";
    }
    SECTION("two cubes close")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-close.obj";
    }
    // SECTION("bunny")
    // {
    //     dhat = 1e-2;
    //     mesh_name = "bunny.obj";
    // }

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = load_mesh(mesh_name, V, E, F);
    REQUIRE(success);

    Constraints constraint_set;
    ipc::construct_constraint_set(/*V_rest=*/V, V, E, F, dhat, constraint_set);
    CAPTURE(mesh_name, dhat);
    CHECK(constraint_set.num_constraints() > 0);

    Eigen::VectorXd grad_b =
        ipc::compute_barrier_potential_gradient(V, E, F, constraint_set, dhat);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return ipc::compute_barrier_potential(
            unflatten(x, V.cols()), E, F, constraint_set, dhat);
    };
    Eigen::VectorXd fgrad_b;
    fd::finite_gradient(flatten(V), f, fgrad_b);

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

TEST_CASE("Test IPC full hessian", "[ipc][hessian]")
{
    double dhat = -1;
    std::string mesh_name = "blah.obj";
    bool ignore_codimensional_vertices = false;

    // SECTION("cube")
    // {
    //     dhat = sqrt(2.0);
    //     mesh_name = "cube.obj";
    // }
    SECTION("two cubes far")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-far.obj";
        ignore_codimensional_vertices = true;
    }
    SECTION("two cubes close")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-close.obj";
        ignore_codimensional_vertices = true;
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
    REQUIRE(success);

    Constraints constraint_set;
    ipc::construct_constraint_set(
        /*V_rest=*/V, V, E, F, dhat, constraint_set,
        ignore_codimensional_vertices);
    CAPTURE(mesh_name, dhat);
    CHECK(constraint_set.num_constraints() > 0);

    Eigen::MatrixXd hess_b = ipc::compute_barrier_potential_hessian(
        V, E, F, constraint_set, dhat, /*project_to_psd=*/false);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return ipc::compute_barrier_potential_gradient(
            unflatten(x, V.cols()), E, F, constraint_set, dhat);
    };
    Eigen::MatrixXd fhess_b;
    fd::finite_jacobian(flatten(V), f, fhess_b);

    REQUIRE(hess_b.squaredNorm() > 1e-3);
    CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
}
