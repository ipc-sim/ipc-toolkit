#include <catch2/catch.hpp>

#include <vector>

#include <finitediff.hpp>
#include <igl/edges.h>

#include <ipc/ipc.hpp>
#include <ipc/friction/friction.hpp>
#include <ipc/utils/logger.hpp>

#include "friction_data_generator.hpp"
#include "../test_utils.hpp"

using namespace ipc;

void check_friction_force_jacobian(
    CollisionMesh mesh,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const Constraints& constraints,
    const double mu,
    const double epsv_times_h,
    const double dhat,
    const double barrier_stiffness,
    const bool recompute_constraints)
{
    REQUIRE(constraints.compute_shape_derivatives);

    const Eigen::MatrixXd& X = mesh.vertices_at_rest();
    double distance_t0 = compute_minimum_distance(mesh, X + Ut, constraints);
    double distance_t1 = compute_minimum_distance(mesh, X + U, constraints);
    // CHECK((distance_t0 < dhat || distance_t1 < dhat));
    if (distance_t0 == 0 || distance_t1 == 0) {
        return;
    }

    CAPTURE(
        mu, epsv_times_h, dhat, barrier_stiffness,
        constraints.vv_constraints.size(), constraints.ev_constraints.size(),
        constraints.ee_constraints.size(), constraints.fv_constraints.size());

    FrictionConstraints friction_constraints;
    construct_friction_constraint_set(
        mesh, X + Ut, constraints, dhat, barrier_stiffness, mu,
        friction_constraints);
    CHECK(friction_constraints.size());

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JPA_wrt_X(mesh.num_vertices(), mesh.ndof());
    for (int i = 0; i < mesh.num_vertices(); i++) {
        JPA_wrt_X.row(i) = Eigen::VectorXd(mesh.vertex_area_gradient(i));
    }
    auto PA_X = [&](const Eigen::VectorXd& x) {
        CollisionMesh fd_mesh(
            fd::unflatten(x, X.cols()), mesh.edges(), mesh.faces());
        return fd_mesh.vertex_areas();
    };
    Eigen::MatrixXd fd_JPA_wrt_X;
    fd::finite_jacobian(fd::flatten(X), PA_X, fd_JPA_wrt_X);
    CHECK(fd::compare_jacobian(JPA_wrt_X, fd_JPA_wrt_X));
    if (!fd::compare_jacobian(JPA_wrt_X, fd_JPA_wrt_X)) {
        print_compare_nonzero(JPA_wrt_X, fd_JPA_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JEA_wrt_X(mesh.num_edges(), mesh.ndof());
    for (int i = 0; i < mesh.num_edges(); i++) {
        JEA_wrt_X.row(i) = Eigen::VectorXd(mesh.edge_area_gradient(i));
    }
    auto EA_X = [&](const Eigen::VectorXd& x) {
        CollisionMesh fd_mesh(
            fd::unflatten(x, X.cols()), mesh.edges(), mesh.faces());
        return fd_mesh.edge_areas();
    };
    Eigen::MatrixXd fd_JEA_wrt_X;
    fd::finite_jacobian(fd::flatten(X), EA_X, fd_JEA_wrt_X);
    CHECK(fd::compare_jacobian(JEA_wrt_X, fd_JEA_wrt_X));
    if (!fd::compare_jacobian(JEA_wrt_X, fd_JEA_wrt_X)) {
        print_compare_nonzero(JEA_wrt_X, fd_JEA_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_X = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraints, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::X);

    auto F_X = [&](const Eigen::VectorXd& x) {
        Eigen::MatrixXd fd_X = fd::unflatten(x, X.cols());

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());

        FrictionConstraints fd_friction_constraints;
        if (recompute_constraints) {
            Constraints fd_constraints;
            fd_constraints.use_convergent_formulation =
                constraints.use_convergent_formulation;
            fd_constraints.compute_shape_derivatives = true;
            fd_constraints.build(fd_mesh, fd_X + Ut, dhat);

            construct_friction_constraint_set(
                fd_mesh, fd_X + Ut, fd_constraints, dhat, barrier_stiffness, mu,
                fd_friction_constraints);
        } else {
            fd_friction_constraints = friction_constraints;
        }

        return compute_friction_force(
            fd_mesh, fd_X, Ut, U, fd_friction_constraints, dhat,
            barrier_stiffness, epsv_times_h);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);
    CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
    if (!fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X)) {
        print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_Ut = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraints, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::Ut);

    auto F_Ut = [&](const Eigen::VectorXd& ut) {
        Eigen::MatrixXd fd_Ut = fd::unflatten(ut, Ut.cols());

        FrictionConstraints fd_friction_constraints;
        if (recompute_constraints) {
            Constraints fd_constraints;
            fd_constraints.use_convergent_formulation =
                constraints.use_convergent_formulation;
            fd_constraints.compute_shape_derivatives = true;
            fd_constraints.build(mesh, X + fd_Ut, dhat);

            construct_friction_constraint_set(
                mesh, X + fd_Ut, fd_constraints, dhat, barrier_stiffness, mu,
                fd_friction_constraints);
        } else {
            fd_friction_constraints = friction_constraints;
        }

        return compute_friction_force(
            mesh, X, fd_Ut, U, friction_constraints, dhat, barrier_stiffness,
            epsv_times_h);
    };
    Eigen::MatrixXd fd_JF_wrt_Ut;
    fd::finite_jacobian(fd::flatten(Ut), F_Ut, fd_JF_wrt_Ut);
    CHECK(fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut));
    if (!fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut)) {
        print_compare_nonzero(JF_wrt_Ut, fd_JF_wrt_Ut);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_U = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraints, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::U);

    auto F_U = [&](const Eigen::VectorXd& u) {
        return compute_friction_force(
            mesh, X, Ut, fd::unflatten(u, U.cols()), friction_constraints, dhat,
            barrier_stiffness, epsv_times_h);
    };
    Eigen::MatrixXd fd_JF_wrt_U;
    fd::finite_jacobian(fd::flatten(U), F_U, fd_JF_wrt_U);
    CHECK(fd::compare_jacobian(JF_wrt_U, fd_JF_wrt_U));
    if (!fd::compare_jacobian(JF_wrt_U, fd_JF_wrt_U)) {
        print_compare_nonzero(JF_wrt_U, fd_JF_wrt_U);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd hess_D = compute_friction_potential_hessian(
        mesh, X + Ut, X + U, friction_constraints, epsv_times_h, false);

    auto grad = [&](const Eigen::VectorXd& u) {
        return compute_friction_potential_gradient(
            mesh, X + Ut, X + fd::unflatten(u, U.cols()), friction_constraints,
            epsv_times_h);
    };
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(fd::flatten(U), grad, fd_hessian);
    CHECK(fd::compare_jacobian(hess_D, fd_hessian));
    if (!fd::compare_jacobian(hess_D, fd_hessian)) {
        print_compare_nonzero(hess_D, fd_hessian);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd force = compute_friction_force(
        mesh, X, Ut, U, friction_constraints, dhat, barrier_stiffness,
        epsv_times_h);
    Eigen::VectorXd grad_D = compute_friction_potential_gradient(
        mesh, X + Ut, X + U, friction_constraints, epsv_times_h);
    CHECK(fd::compare_gradient(-force, grad_D));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd jac_force = compute_friction_force_jacobian(
        mesh, X, Ut, U, friction_constraints, dhat, barrier_stiffness,
        epsv_times_h, FrictionConstraint::DiffWRT::U);
    CHECK(fd::compare_jacobian(-jac_force, hess_D));
}

TEST_CASE("Test friction force jacobian", "[friction][force-jacobian][thisone]")
{
    const int x_case = GENERATE(0, 1);
    FrictionData data = friction_data_generator();
    const auto& [V0, V1, E, F, constraints, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;
    REQUIRE(constraints.compute_shape_derivatives);

    Eigen::MatrixXd X, Ut, U;
    switch (x_case) {
    case 0:
        X = V0;
        break;
    case 1:
    default:
        X = V0 - (V1 - V0);
        break;
    }
    Ut = V0 - X;
    U = V1 - X;

    CollisionMesh mesh(X, E, F);

    check_friction_force_jacobian(
        mesh, Ut, U, constraints, mu, epsv_times_h, dhat, barrier_stiffness,
        false);
}

TEST_CASE(
    "Test friction force jacobian on real data",
    "[friction][force-jacobian][real-data]")
{
    bool use_convergent_formulation = GENERATE(true, false);

    std::string scene;
    bool is_2D = true;
    double mu, dhat, kappa, epsv_dt;
    SECTION("point-plane")
    {
        scene = "point-plane";
        mu = 0.5;
        dhat = 0.1;
        kappa = 141;
        epsv_dt = 5e-6;
    }
    SECTION("square-circle")
    {
        scene = "square-circle";
        mu = 0.5;
        dhat = 1e-3;
        // kappa = 67873353;
        kappa = 67873353 / 10;
        epsv_dt = 1e-4;
    }
    SECTION("square-circle-dense")
    {
        scene = "square-circle-dense";
        mu = 0.5;
        dhat = 1e-2;
        kappa = GENERATE(8.6e9, 1e6);
        epsv_dt = 1.5e-5;
    }
    // SECTION("square-incline")
    // {
    //     scene = "square-incline";
    //     mu = 0.5;
    //     dhat = 0.1;
    //     kappa = 141;
    //     epsv_dt = 5e-6;
    // }

    CAPTURE(scene, mu, dhat, kappa, epsv_dt, use_convergent_formulation);

    Eigen::MatrixXd X, Ut, U;
    Eigen::MatrixXi E, F;
    {
        X = loadMarketXd(fmt::format(
            "{}friction-force-jacobian/{}/X.mtx", TEST_DATA_DIR, scene));
        Ut = loadMarketXd(fmt::format(
            "{}friction-force-jacobian/{}/Ut.mtx", TEST_DATA_DIR, scene));
        Ut = fd::unflatten(Ut, X.cols());
        U = loadMarketXd(fmt::format(
            "{}friction-force-jacobian/{}/U.mtx", TEST_DATA_DIR, scene));
        U = fd::unflatten(U, X.cols());
        if (is_2D) {
            E = loadMarketXi(fmt::format(
                "{}friction-force-jacobian/{}/F.mtx", TEST_DATA_DIR, scene));
        } else {
            F = loadMarketXi(fmt::format(
                "{}friction-force-jacobian/{}/F.mtx", TEST_DATA_DIR, scene));
            igl::edges(F, E);
        }
    }

    std::vector<bool> is_on_surface =
        CollisionMesh::construct_is_on_surface(X.rows(), E);
    CollisionMesh mesh(is_on_surface, X, E, F);

    X = mesh.vertices(X);
    if (Ut.rows() != X.rows()) {
        Ut = mesh.vertices(Ut);
    }
    if (U.rows() != X.rows()) {
        U = mesh.vertices(U);
    }

    Constraints constraints;
    constraints.use_convergent_formulation = use_convergent_formulation;
    constraints.compute_shape_derivatives = true;
    constraints.build(mesh, X + Ut, dhat);

    REQUIRE(constraints.compute_shape_derivatives);

    CHECK(compute_minimum_distance(mesh, X + Ut, constraints) != 0);
    CHECK(compute_minimum_distance(mesh, X + U, constraints) != 0);

    check_friction_force_jacobian(
        mesh, Ut, U, constraints, mu, epsv_dt, dhat, kappa, true);
}
