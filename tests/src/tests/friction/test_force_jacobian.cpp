#include <tests/config.hpp>
#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/ipc.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>


#include <finitediff.hpp>
#include <igl/edges.h>

using namespace ipc;

void check_friction_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const NormalCollisions& collisions,
    const double mu,
    const double epsv_times_h,
    const double dhat,
    const double barrier_stiffness,
    const bool recompute_collisions)
{
    REQUIRE(collisions.enable_shape_derivatives());

    const Eigen::MatrixXd& X = mesh.rest_positions();
    double distance_t0 = collisions.compute_minimum_distance(mesh, X + Ut);
    double distance_t1 = collisions.compute_minimum_distance(mesh, X + U);
    // CHECK((distance_t0 < dhat || distance_t1 < dhat));
    if (distance_t0 == 0 || distance_t1 == 0) {
        return;
    }

    const Eigen::MatrixXd velocities = U - Ut;

    CAPTURE(
        mu, epsv_times_h, dhat, barrier_stiffness,
        collisions.vv_collisions.size(), collisions.ev_collisions.size(),
        collisions.ee_collisions.size(), collisions.fv_collisions.size());

    TangentialCollisions tangential_collisions;
    tangential_collisions.build(
        mesh, X + Ut, collisions, BarrierPotential(dhat), barrier_stiffness,
        mu);
    CHECK(tangential_collisions.size());

    const FrictionPotential D(epsv_times_h);

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

    CHECKED_ELSE(fd::compare_jacobian(JPA_wrt_X, fd_JPA_wrt_X))
    {
        tests::print_compare_nonzero(JPA_wrt_X, fd_JPA_wrt_X);
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

    CHECKED_ELSE(fd::compare_jacobian(JEA_wrt_X, fd_JEA_wrt_X))
    {
        tests::print_compare_nonzero(JEA_wrt_X, fd_JEA_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_X = D.force_jacobian(
        tangential_collisions, mesh, X, Ut, velocities, BarrierPotential(dhat),
        barrier_stiffness, FrictionPotential::DiffWRT::REST_POSITIONS);

    auto F_X = [&](const Eigen::VectorXd& x) {
        Eigen::MatrixXd fd_X = fd::unflatten(x, X.cols());

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());
        fd_mesh.init_area_jacobians();

        TangentialCollisions fd_friction_collisions;
        if (recompute_collisions) {
            NormalCollisions fd_collisions;
            fd_collisions.set_use_area_weighting(
                collisions.use_area_weighting());
            fd_collisions.set_use_improved_max_approximator(
                collisions.use_improved_max_approximator());
            fd_collisions.set_enable_shape_derivatives(true);
            fd_collisions.build(fd_mesh, fd_X + Ut, dhat);

            fd_friction_collisions.build(
                fd_mesh, fd_X + Ut, fd_collisions, BarrierPotential(dhat),
                barrier_stiffness, mu);
        } else {
            fd_friction_collisions = tangential_collisions;
        }

        return D.force(
            fd_friction_collisions, fd_mesh, fd_X, Ut, velocities,
            BarrierPotential(dhat), barrier_stiffness);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(fd::flatten(X), F_X, fd_JF_wrt_X);

    CHECKED_ELSE(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X))
    {
        tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_Ut = D.force_jacobian(
        tangential_collisions, mesh, X, Ut, velocities, BarrierPotential(dhat),
        barrier_stiffness, FrictionPotential::DiffWRT::LAGGED_DISPLACEMENTS);

    auto F_Ut = [&](const Eigen::VectorXd& ut) {
        Eigen::MatrixXd fd_Ut = fd::unflatten(ut, Ut.cols());

        TangentialCollisions fd_friction_collisions;
        if (recompute_collisions) {
            NormalCollisions fd_collisions;
            fd_collisions.set_use_area_weighting(
                collisions.use_area_weighting());
            fd_collisions.set_use_improved_max_approximator(
                collisions.use_improved_max_approximator());
            fd_collisions.set_enable_shape_derivatives(true);
            fd_collisions.build(mesh, X + fd_Ut, dhat);

            fd_friction_collisions.build(
                mesh, X + fd_Ut, fd_collisions, BarrierPotential(dhat),
                barrier_stiffness, mu);
        } else {
            fd_friction_collisions = tangential_collisions;
        }

        return D.force(
            tangential_collisions, mesh, X, fd_Ut, velocities,
            BarrierPotential(dhat), barrier_stiffness);
    };
    Eigen::MatrixXd fd_JF_wrt_Ut;
    fd::finite_jacobian(fd::flatten(Ut), F_Ut, fd_JF_wrt_Ut);

    CHECKED_ELSE(fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut))
    {
        tests::print_compare_nonzero(JF_wrt_Ut, fd_JF_wrt_Ut);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_V = D.force_jacobian(
        tangential_collisions, mesh, X, Ut, velocities, BarrierPotential(dhat),
        barrier_stiffness, FrictionPotential::DiffWRT::VELOCITIES);

    auto F_V = [&](const Eigen::VectorXd& v) {
        return D.force(
            tangential_collisions, mesh, X, Ut,
            fd::unflatten(v, velocities.cols()), BarrierPotential(dhat),
            barrier_stiffness);
    };
    Eigen::MatrixXd fd_JF_wrt_V;
    fd::finite_jacobian(fd::flatten(velocities), F_V, fd_JF_wrt_V);

    CHECKED_ELSE(fd::compare_jacobian(JF_wrt_V, fd_JF_wrt_V))
    {
        tests::print_compare_nonzero(JF_wrt_V, fd_JF_wrt_V);
    }

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::MatrixXd hess_D =
        D.hessian(tangential_collisions, mesh, velocities);

    auto grad = [&](const Eigen::VectorXd& v) {
        return D.gradient(
            tangential_collisions, mesh, fd::unflatten(v, velocities.cols()));
    };
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(fd::flatten(velocities), grad, fd_hessian);

    CHECKED_ELSE(fd::compare_jacobian(hess_D, fd_hessian))
    {
        tests::print_compare_nonzero(hess_D, fd_hessian);
    }

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::VectorXd force = D.force(
        tangential_collisions, mesh, X, Ut, velocities, BarrierPotential(dhat),
        barrier_stiffness);
    const Eigen::VectorXd grad_D =
        D.gradient(tangential_collisions, mesh, velocities);
    CHECK(fd::compare_gradient(-force, grad_D));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd jac_force = D.force_jacobian(
        tangential_collisions, mesh, X, Ut, velocities, BarrierPotential(dhat),
        barrier_stiffness, FrictionPotential::DiffWRT::VELOCITIES);
    CHECK(fd::compare_jacobian(-jac_force, hess_D));
}

TEST_CASE("Friction force jacobian", "[friction][force-jacobian]")
{
    const int x_case = GENERATE(0, 1);
    FrictionData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;
    REQUIRE(collisions.enable_shape_derivatives());

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
    mesh.init_area_jacobians();

    check_friction_force_jacobian(
        mesh, Ut, U, collisions, mu, epsv_times_h, dhat, barrier_stiffness,
        false);
}

TEST_CASE(
    "Friction force jacobian on real data",
    "[friction][force-jacobian][real-data]")
{
    bool use_area_weighting = GENERATE(true, false);
    bool use_improved_max_approximator = GENERATE(true, false);

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
        kappa = 67873353 / 10 * dhat;
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

    CAPTURE(
        scene, mu, dhat, kappa, epsv_dt, use_area_weighting,
        use_improved_max_approximator);

    Eigen::MatrixXd X, Ut, U;
    Eigen::MatrixXi E, F;
    {
        const auto dir = tests::DATA_DIR / "friction-force-jacobian" / scene;
        X = tests::loadMarketXd((dir / "X.mtx").string());
        Ut = tests::loadMarketXd((dir / "Ut.mtx").string());
        Ut = fd::unflatten(Ut, X.cols());
        U = tests::loadMarketXd((dir / "U.mtx").string());
        U = fd::unflatten(U, X.cols());
        if (is_2D) {
            E = tests::loadMarketXi((dir / "F.mtx").string());
        } else {
            F = tests::loadMarketXi((dir / "F.mtx").string());
            igl::edges(F, E);
        }
    }

    std::vector<bool> is_on_surface =
        CollisionMesh::construct_is_on_surface(X.rows(), E);
    CollisionMesh mesh(is_on_surface, X, E, F);
    mesh.init_area_jacobians();

    X = mesh.vertices(X);
    if (Ut.rows() != X.rows()) {
        Ut = mesh.vertices(Ut);
    }
    if (U.rows() != X.rows()) {
        U = mesh.vertices(U);
    }

    NormalCollisions collisions;
    collisions.set_use_area_weighting(use_area_weighting);
    collisions.set_use_improved_max_approximator(use_improved_max_approximator);
    collisions.set_enable_shape_derivatives(true);
    collisions.build(mesh, X + Ut, dhat);

    REQUIRE(collisions.enable_shape_derivatives());

    CHECK(collisions.compute_minimum_distance(mesh, X + Ut) != 0);
    CHECK(collisions.compute_minimum_distance(mesh, X + U) != 0);

    check_friction_force_jacobian(
        mesh, Ut, U, collisions, mu, epsv_dt, dhat, kappa, true);
}
