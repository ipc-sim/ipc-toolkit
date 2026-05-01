#include <tests/config.hpp>
#include <tests/dof_layout.hpp>
#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/ipc.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/smooth_contact/smooth_contact_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>

using namespace ipc;

/// How tangential collisions get scalar μ before update_lagged_* (isotropic
/// uses build defaults; matchstick sets ellipse axes then lagged refresh).
enum class FrictionJacobianMuSetup { Isotropic, MatchstickAxes };

void check_friction_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const NormalCollisions& collisions,
    const double mu_s,
    const double mu_k,
    const double epsv_times_h,
    const double dhat,
    const double barrier_stiffness,
    const bool recompute_collisions,
    const FrictionJacobianMuSetup mu_setup = FrictionJacobianMuSetup::Isotropic)
{
    REQUIRE(collisions.enable_shape_derivatives());

    const BarrierPotential B(dhat, barrier_stiffness);
    const FrictionPotential D(epsv_times_h);

    const Eigen::MatrixXd& X = mesh.rest_positions();

    // Ensure Ut and U match the mesh size (map as displacements, not positions)
    Eigen::MatrixXd Ut_mesh =
        Ut.rows() == mesh.num_vertices() ? Ut : mesh.map_displacements(Ut);
    Eigen::MatrixXd U_mesh =
        U.rows() == mesh.num_vertices() ? U : mesh.map_displacements(U);

    double distance_t0 = collisions.compute_minimum_distance(mesh, X + Ut_mesh);
    double distance_t1 = collisions.compute_minimum_distance(mesh, X + U_mesh);
    // CHECK((distance_t0 < dhat || distance_t1 < dhat));
    if (distance_t0 == 0 || distance_t1 == 0) {
        return;
    }

    // V = (X + U) - (X + Ut) = U - Ut
    const Eigen::MatrixXd velocities = U_mesh - Ut_mesh;

    CAPTURE(
        mu_s, mu_k, epsv_times_h, dhat, barrier_stiffness,
        collisions.vv_collisions.size(), collisions.ev_collisions.size(),
        collisions.ee_collisions.size(), collisions.fv_collisions.size());

    TangentialCollisions tangential_collisions;
    tangential_collisions.build(mesh, X + Ut_mesh, collisions, B, mu_s, mu_k);
    CHECK(!tangential_collisions.empty());

    if (mu_setup == FrictionJacobianMuSetup::MatchstickAxes) {
        for (size_t i = 0; i < tangential_collisions.size(); ++i) {
            tangential_collisions[i].mu_s_aniso = Eigen::Vector2d(0.75, 0.35);
            tangential_collisions[i].mu_k_aniso = Eigen::Vector2d(0.55, 0.28);
        }
    }
    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, X, Ut_mesh, velocities);

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JPA_wrt_X(mesh.num_vertices(), mesh.ndof());
    for (int i = 0; i < mesh.num_vertices(); i++) {
        JPA_wrt_X.row(i) = Eigen::VectorXd(mesh.vertex_area_gradient(i));
    }
    auto PA_X = [&](const Eigen::VectorXd& x) {
        CollisionMesh fd_mesh(
            tests::unflatten(x, X.cols()), mesh.edges(), mesh.faces());
        return fd_mesh.vertex_areas();
    };
    Eigen::MatrixXd fd_JPA_wrt_X;
    fd::finite_jacobian(tests::flatten(X), PA_X, fd_JPA_wrt_X);

    CHECK(fd::compare_jacobian(JPA_wrt_X, fd_JPA_wrt_X));
    if (!fd::compare_jacobian(JPA_wrt_X, fd_JPA_wrt_X)) {
        tests::print_compare_nonzero(JPA_wrt_X, fd_JPA_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JEA_wrt_X(mesh.num_edges(), mesh.ndof());
    for (int i = 0; i < mesh.num_edges(); i++) {
        JEA_wrt_X.row(i) = Eigen::VectorXd(mesh.edge_area_gradient(i));
    }
    auto EA_X = [&](const Eigen::VectorXd& x) {
        CollisionMesh fd_mesh(
            tests::unflatten(x, X.cols()), mesh.edges(), mesh.faces());
        return fd_mesh.edge_areas();
    };
    Eigen::MatrixXd fd_JEA_wrt_X;
    fd::finite_jacobian(tests::flatten(X), EA_X, fd_JEA_wrt_X);

    CHECK(fd::compare_jacobian(JEA_wrt_X, fd_JEA_wrt_X));
    if (!fd::compare_jacobian(JEA_wrt_X, fd_JEA_wrt_X)) {
        tests::print_compare_nonzero(JEA_wrt_X, fd_JEA_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_X = D.force_jacobian(
        tangential_collisions, mesh, X, Ut_mesh, velocities, B,
        FrictionPotential::DiffWRT::REST_POSITIONS);

    auto F_X = [&](const Eigen::VectorXd& x) {
        Eigen::MatrixXd fd_X = tests::unflatten(x, X.cols());

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());
        fd_mesh.init_area_jacobians();

        // Ensure Ut_mesh and velocities match fd_mesh size
        // Since fd_X is created from X (which is filtered), fd_mesh should have
        // the same number of vertices as the original mesh, so Ut_mesh and
        // velocities should already match. But check to be safe.
        assert(fd_mesh.num_vertices() == mesh.num_vertices());
        assert(Ut_mesh.rows() == fd_mesh.num_vertices());
        assert(velocities.rows() == fd_mesh.num_vertices());

        TangentialCollisions fd_friction_collisions;
        if (recompute_collisions) {
            NormalCollisions fd_collisions;
            fd_collisions.set_use_area_weighting(
                collisions.use_area_weighting());
            fd_collisions.set_collision_set_type(
                collisions.collision_set_type());
            fd_collisions.set_enable_shape_derivatives(true);
            fd_collisions.build(fd_mesh, fd_X + Ut_mesh, dhat);

            fd_friction_collisions.build(
                fd_mesh, fd_X + Ut_mesh, fd_collisions, B, mu_s, mu_k);
            fd_friction_collisions
                .update_lagged_anisotropic_friction_coefficients(
                    fd_mesh, fd_X, Ut_mesh, velocities);
        } else {
            fd_friction_collisions = tangential_collisions;
        }

        return D.force(
            fd_friction_collisions, fd_mesh, fd_X, Ut_mesh, velocities, B);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(tests::flatten(X), F_X, fd_JF_wrt_X);

    CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
    if (!fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X)) {
        tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    }

    ///////////////////////////////////////////////////////////////////////////
    Eigen::MatrixXd JF_wrt_Ut = D.force_jacobian(
        tangential_collisions, mesh, X, Ut_mesh, velocities, B,
        FrictionPotential::DiffWRT::LAGGED_DISPLACEMENTS);

    auto F_Ut = [&](const Eigen::VectorXd& ut) {
        Eigen::MatrixXd fd_Ut = tests::unflatten(ut, Ut_mesh.cols());

        TangentialCollisions fd_friction_collisions;
        if (recompute_collisions) {
            NormalCollisions fd_collisions;
            fd_collisions.set_use_area_weighting(
                collisions.use_area_weighting());
            fd_collisions.set_collision_set_type(
                collisions.collision_set_type());
            fd_collisions.set_enable_shape_derivatives(true);
            fd_collisions.build(mesh, X + fd_Ut, dhat);

            fd_friction_collisions.build(
                mesh, X + fd_Ut, fd_collisions, B, mu_s, mu_k);
            fd_friction_collisions
                .update_lagged_anisotropic_friction_coefficients(
                    mesh, X, fd_Ut, velocities);
        } else {
            fd_friction_collisions = tangential_collisions;
        }

        return D.force(fd_friction_collisions, mesh, X, fd_Ut, velocities, B);
    };
    Eigen::MatrixXd fd_JF_wrt_Ut;
    fd::finite_jacobian(tests::flatten(Ut_mesh), F_Ut, fd_JF_wrt_Ut);

    CHECK(fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut));
    if (!fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut)) {
        tests::print_compare_nonzero(JF_wrt_Ut, fd_JF_wrt_Ut);
    }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_V = D.force_jacobian(
        tangential_collisions, mesh, X, Ut_mesh, velocities, B,
        FrictionPotential::DiffWRT::VELOCITIES);

    auto F_V = [&](const Eigen::VectorXd& v) {
        return D.force(
            tangential_collisions, mesh, X, Ut_mesh,
            tests::unflatten(v, velocities.cols()), B);
    };
    Eigen::MatrixXd fd_JF_wrt_V;
    fd::finite_jacobian(tests::flatten(velocities), F_V, fd_JF_wrt_V);

    CHECK(fd::compare_jacobian(JF_wrt_V, fd_JF_wrt_V));
    if (!fd::compare_jacobian(JF_wrt_V, fd_JF_wrt_V)) {
        tests::print_compare_nonzero(JF_wrt_V, fd_JF_wrt_V);
    }

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::MatrixXd hess_D =
        D.hessian(tangential_collisions, mesh, velocities);

    auto grad = [&](const Eigen::VectorXd& v) {
        return D.gradient(
            tangential_collisions, mesh,
            tests::unflatten(v, velocities.cols()));
    };
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(tests::flatten(velocities), grad, fd_hessian);

    CHECK(fd::compare_jacobian(hess_D, fd_hessian));
    if (!fd::compare_jacobian(hess_D, fd_hessian)) {
        tests::print_compare_nonzero(hess_D, fd_hessian);
    }

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::VectorXd force =
        D.force(tangential_collisions, mesh, X, Ut_mesh, velocities, B);
    const Eigen::VectorXd grad_D =
        D.gradient(tangential_collisions, mesh, velocities);
    CHECK(fd::compare_gradient(-force, grad_D));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd jac_force = D.force_jacobian(
        tangential_collisions, mesh, X, Ut_mesh, velocities, B,
        FrictionPotential::DiffWRT::VELOCITIES);
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
        mesh, Ut, U, collisions, 0.5 * mu, mu, epsv_times_h, dhat,
        barrier_stiffness, false);
}

TEST_CASE(
    "Friction force jacobian lagged matchstick mu",
    "[friction][force-jacobian][anisotropic]")
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
        mesh, Ut, U, collisions, 0.5 * mu, mu, epsv_times_h, dhat,
        barrier_stiffness, false, FrictionJacobianMuSetup::MatchstickAxes);
}

TEST_CASE(
    "Friction force jacobian on real data",
    "[friction][force-jacobian][real-data]")
{
    bool use_area_weighting = GENERATE(true, false);
    const NormalCollisions::CollisionSetType collision_set_type = GENERATE(
        NormalCollisions::CollisionSetType::IPC,
        NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX,
        NormalCollisions::CollisionSetType::OGC);

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
        collision_set_type);

    Eigen::MatrixXd X, Ut, U;
    Eigen::MatrixXi E, F;
    {
        const auto dir = tests::DATA_DIR / "friction-force-jacobian" / scene;
        X = tests::loadMarketXd((dir / "X.mtx").string());
        Ut = tests::loadMarketXd((dir / "Ut.mtx").string());
        Ut = fd::unflatten(Ut, X.cols()); // data files are always RowMajor
        U = tests::loadMarketXd((dir / "U.mtx").string());
        U = fd::unflatten(U, X.cols()); // data files are always RowMajor
        if (is_2D) {
            E = tests::loadMarketXi((dir / "F.mtx").string());
        } else {
            F = tests::loadMarketXi((dir / "F.mtx").string());
            igl::edges(F, E);
        }
    }

    std::vector<bool> is_on_surface =
        CollisionMesh::construct_is_on_surface(X.rows(), E);
    CollisionMesh mesh(
        is_on_surface, std::vector<bool>(X.rows(), false), X, E, F);
    mesh.init_area_jacobians();

    X = mesh.vertices(X);
    if (Ut.rows() != X.rows()) {
        Ut = mesh.map_displacements(Ut);
    }
    if (U.rows() != X.rows()) {
        U = mesh.map_displacements(U);
    }

    NormalCollisions collisions;
    collisions.set_use_area_weighting(use_area_weighting);
    collisions.set_collision_set_type(collision_set_type);
    collisions.set_enable_shape_derivatives(true);
    collisions.build(mesh, X + Ut, dhat);

    REQUIRE(collisions.enable_shape_derivatives());

    CHECK(collisions.compute_minimum_distance(mesh, X + Ut) != 0);
    CHECK(collisions.compute_minimum_distance(mesh, X + U) != 0);

    check_friction_force_jacobian(
        mesh, Ut, U, collisions, 0.5 * mu, mu, epsv_dt, dhat, kappa, true);
}

void check_smooth_friction_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const SmoothCollisions& collisions,
    const double mu,
    const double epsv_times_h,
    const SmoothContactParameters& params,
    const double barrier_stiffness,
    const bool recompute_collisions)
{
    const int dim = mesh.dim();
    const double dhat = params.dhat;
    const Eigen::MatrixXd& X = mesh.rest_positions();

    // Ensure Ut and U match the mesh size
    Eigen::MatrixXd Ut_mesh =
        Ut.rows() == mesh.num_vertices() ? Ut : mesh.map_displacements(Ut);
    Eigen::MatrixXd U_mesh =
        U.rows() == mesh.num_vertices() ? U : mesh.map_displacements(U);

    double distance_t0 = collisions.compute_minimum_distance(mesh, X + Ut_mesh);
    double distance_t1 = collisions.compute_minimum_distance(mesh, X + U_mesh);
    // CHECK((distance_t0 < dhat || distance_t1 < dhat));
    if (distance_t0 == 0 || distance_t1 == 0) {
        return;
    }

    Eigen::MatrixXd velocities = U_mesh - Ut_mesh;

    CAPTURE(mu, epsv_times_h, dhat, barrier_stiffness, collisions.size());

    TangentialCollisions friction_collisions;
    friction_collisions.build(
        mesh, X + Ut_mesh, collisions, params, barrier_stiffness,
        Eigen::VectorXd::Ones(mesh.num_vertices()) * mu,
        Eigen::VectorXd::Ones(mesh.num_vertices()) * mu);
    CHECK(!friction_collisions.empty());
    friction_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, X, Ut_mesh, velocities);

    const FrictionPotential D(epsv_times_h);

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::VectorXd force = D.smooth_contact_force(
        friction_collisions, mesh, X, Ut_mesh, velocities);
    const Eigen::VectorXd grad_D =
        D.gradient(friction_collisions, mesh, velocities);
    CHECK((force + grad_D).norm() <= 1e-8 * force.norm());

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::MatrixXd hess_D =
        D.hessian(friction_collisions, mesh, velocities);

    auto grad = [&](const Eigen::VectorXd& v) {
        return D.gradient(
            friction_collisions, mesh, tests::unflatten(v, velocities.cols()));
    };
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        tests::flatten(velocities), grad, fd_hessian, fd::AccuracyOrder::FOURTH,
        1e-6 * dhat);
    // CHECK(fd::compare_jacobian(hess_D, fd_hessian));
    // if (!fd::compare_jacobian(hess_D, fd_hessian)) {
    //     tests::print_compare_nonzero(hess_D, fd_hessian);
    // }
    CHECK(
        (hess_D.norm() == 0
         || (hess_D - fd_hessian).norm() <= 1e-7 * hess_D.norm()));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd jac_force = D.smooth_contact_force_jacobian(
        friction_collisions, mesh, X, Ut_mesh, velocities, params,
        FrictionPotential::DiffWRT::VELOCITIES);
    CHECK((hess_D + jac_force).norm() <= 1e-7 * hess_D.norm());

    ///////////////////////////////////////////////////////////////////////////

    auto create_smooth_collision = [&](const CollisionMesh& fd_mesh,
                                       const Eigen::MatrixXd&
                                           fd_lagged_positions) {
        SmoothCollisions fd_collisions;
        assert(friction_collisions.size() == 1);

        auto cc = friction_collisions[0].smooth_collision;
        std::shared_ptr<SmoothCollision> fd_cc;
        if (dim == 3) {
            if (cc->type() == CollisionType::EDGE_EDGE) {
                fd_cc = std::make_shared<SmoothCollisionTemplate<Edge3, Edge3>>(
                    (*cc)[0], (*cc)[1],
                    PrimitiveDistType<Edge3, Edge3>::type::AUTO, fd_mesh,
                    params, dhat, fd_lagged_positions);
            } else if (cc->type() == CollisionType::EDGE_VERTEX) {
                fd_cc =
                    std::make_shared<SmoothCollisionTemplate<Edge3, Point3>>(
                        (*cc)[0], (*cc)[1],
                        PrimitiveDistType<Edge3, Point3>::type::AUTO, fd_mesh,
                        params, dhat, fd_lagged_positions);
            } else if (cc->type() == CollisionType::VERTEX_VERTEX) {
                fd_cc =
                    std::make_shared<SmoothCollisionTemplate<Point3, Point3>>(
                        (*cc)[0], (*cc)[1],
                        PrimitiveDistType<Point3, Point3>::type::AUTO, fd_mesh,
                        params, dhat, fd_lagged_positions);
            } else if (cc->type() == CollisionType::FACE_VERTEX) {
                fd_cc = std::make_shared<SmoothCollisionTemplate<Face, Point3>>(
                    (*cc)[0], (*cc)[1],
                    PrimitiveDistType<Face, Point3>::type::AUTO, fd_mesh,
                    params, dhat, fd_lagged_positions);
            }

            fd_collisions.collisions.push_back(fd_cc);
        } else {
            if (cc->type() == CollisionType::EDGE_VERTEX) {
                fd_cc =
                    std::make_shared<SmoothCollisionTemplate<Edge2, Point2>>(
                        (*cc)[0], (*cc)[1],
                        PrimitiveDistType<Edge2, Point2>::type::AUTO, fd_mesh,
                        params, dhat, fd_lagged_positions);
            } else if (cc->type() == CollisionType::VERTEX_VERTEX) {
                fd_cc =
                    std::make_shared<SmoothCollisionTemplate<Point2, Point2>>(
                        (*cc)[0], (*cc)[1],
                        PrimitiveDistType<Point2, Point2>::type::AUTO, fd_mesh,
                        params, dhat, fd_lagged_positions);
            }

            fd_collisions.collisions.push_back(fd_cc);
        }

        return fd_collisions;
    };

    ///////////////////////////////////////////////////////////////////////////

    // test contact force norm derivative
    //{
    //    const Eigen::MatrixXd lagged_positions = X + Ut;
    //    Eigen::VectorXd normal_force_jacobian =
    //    Eigen::VectorXd::Zero(X.size());
    //    {
    //        auto cc = create_smooth_collision(mesh, lagged_positions);
    //        SmoothContactPotential<SmoothCollisions> potential(params);
    //        Eigen::VectorXd g = potential.gradient(cc, mesh,
    //        lagged_positions); Eigen::SparseMatrix<double> h =
    //        potential.hessian(cc, mesh, lagged_positions);
    //        normal_force_jacobian = (h * g) / g.norm();
    //    }

    //    // finite difference
    //    auto F_X = [&](const Eigen::VectorXd& x) {
    //        Eigen::MatrixXd fd_X = fd::unflatten(x, dim);
    //        const Eigen::MatrixXd fd_lagged_positions = fd_X + Ut;

    //        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());
    //        auto fd_cc = create_smooth_collision(fd_mesh,
    //        fd_lagged_positions);

    //        SmoothContactPotential<SmoothCollisions> potential(params);
    //        return potential.gradient(fd_cc, fd_mesh,
    //        fd_lagged_positions).norm();
    //    };

    //    Eigen::VectorXd fd_normal_force_jacobian;
    //    fd::finite_gradient(fd::flatten(X), F_X, fd_normal_force_jacobian);
    //    CHECK((normal_force_jacobian - fd_normal_force_jacobian).norm() <=
    //    1e-7 * std::max(normal_force_jacobian.norm(), 1e-8));
    //}

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_X = D.smooth_contact_force_jacobian(
        friction_collisions, mesh, X, Ut_mesh, velocities, params,
        FrictionPotential::DiffWRT::REST_POSITIONS);

    auto F_X = [&](const Eigen::VectorXd& x) {
        Eigen::MatrixXd fd_X = tests::unflatten(x, X.cols());
        Eigen::MatrixXd fd_lagged_positions = fd_X + Ut_mesh;

        CollisionMesh fd_mesh(fd_X, mesh.edges(), mesh.faces());

        // Ensure Ut_mesh and velocities match fd_mesh size
        // Since fd_X is created from X (which is filtered), fd_mesh should have
        // the same number of vertices as the original mesh
        assert(fd_mesh.num_vertices() == mesh.num_vertices());
        assert(Ut_mesh.rows() == fd_mesh.num_vertices());
        assert(velocities.rows() == fd_mesh.num_vertices());

        auto fd_collisions =
            create_smooth_collision(fd_mesh, fd_lagged_positions);

        TangentialCollisions fd_friction_collisions;
        fd_friction_collisions.build(
            fd_mesh, fd_lagged_positions, fd_collisions, params,
            barrier_stiffness,
            Eigen::VectorXd::Ones(fd_mesh.num_vertices()) * mu,
            Eigen::VectorXd::Ones(fd_mesh.num_vertices()) * mu);
        fd_friction_collisions.update_lagged_anisotropic_friction_coefficients(
            fd_mesh, fd_X, Ut_mesh, velocities);

        return D.smooth_contact_force(
            fd_friction_collisions, fd_mesh, fd_X, Ut_mesh, velocities);
    };
    Eigen::MatrixXd fd_JF_wrt_X;
    fd::finite_jacobian(
        tests::flatten(X), F_X, fd_JF_wrt_X, fd::AccuracyOrder::FOURTH,
        1e-6 * dhat);
    // CHECK(fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X));
    // if (!fd::compare_jacobian(JF_wrt_X, fd_JF_wrt_X)) {
    //     tests::print_compare_nonzero(JF_wrt_X, fd_JF_wrt_X);
    // }
    CHECK(
        (JF_wrt_X - fd_JF_wrt_X).norm()
        <= 1e-7 * std::max(JF_wrt_X.norm(), 1e-8));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_Ut = D.smooth_contact_force_jacobian(
        friction_collisions, mesh, X, Ut_mesh, velocities, params,
        FrictionPotential::DiffWRT::LAGGED_DISPLACEMENTS);

    auto F_Ut = [&](const Eigen::VectorXd& ut) {
        Eigen::MatrixXd fd_Ut = tests::unflatten(ut, Ut_mesh.cols());
        Eigen::MatrixXd fd_lagged_positions = X + fd_Ut;

        auto fd_collisions = create_smooth_collision(mesh, fd_lagged_positions);

        TangentialCollisions fd_friction_collisions;
        fd_friction_collisions.build(
            mesh, fd_lagged_positions, fd_collisions, params, barrier_stiffness,
            Eigen::VectorXd::Ones(mesh.num_vertices()) * mu,
            Eigen::VectorXd::Ones(mesh.num_vertices()) * mu);
        fd_friction_collisions.update_lagged_anisotropic_friction_coefficients(
            mesh, X, fd_Ut, velocities);

        return D.smooth_contact_force(
            fd_friction_collisions, mesh, X, fd_Ut, velocities);
    };
    Eigen::MatrixXd fd_JF_wrt_Ut;
    fd::finite_jacobian(
        tests::flatten(Ut_mesh), F_Ut, fd_JF_wrt_Ut, fd::AccuracyOrder::FOURTH,
        1e-6 * dhat);
    // CHECK(fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut));
    // if (!fd::compare_jacobian(JF_wrt_Ut, fd_JF_wrt_Ut)) {
    //     tests::print_compare_nonzero(JF_wrt_Ut, fd_JF_wrt_Ut);
    // }
    CHECK(
        (JF_wrt_Ut - fd_JF_wrt_Ut).norm()
        <= 1e-7 * std::max(JF_wrt_Ut.norm(), 1e-8));

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd JF_wrt_V = D.smooth_contact_force_jacobian(
        friction_collisions, mesh, X, Ut_mesh, velocities, params,
        FrictionPotential::DiffWRT::VELOCITIES);

    auto F_V = [&](const Eigen::VectorXd& v) {
        return D.smooth_contact_force(
            friction_collisions, mesh, X, Ut_mesh,
            tests::unflatten(v, velocities.cols()));
    };
    Eigen::MatrixXd fd_JF_wrt_V;
    fd::finite_jacobian(
        tests::flatten(velocities), F_V, fd_JF_wrt_V, fd::AccuracyOrder::FOURTH,
        1e-6 * dhat);
    CHECK(
        (JF_wrt_V.norm() == 0
         || (fd_JF_wrt_V - JF_wrt_V).norm() <= 1e-7 * JF_wrt_V.norm()));

    if (mesh.dim() == 3) {
        for (size_t i = 0; i < friction_collisions.size(); ++i) {
            friction_collisions[i].mu_s_aniso = Eigen::Vector2d(0.8, 0.45);
            friction_collisions[i].mu_k_aniso = Eigen::Vector2d(0.65, 0.32);
        }
        friction_collisions.update_lagged_anisotropic_friction_coefficients(
            mesh, X, Ut_mesh, velocities);

        Eigen::MatrixXd JF_wrt_V_aniso = D.smooth_contact_force_jacobian(
            friction_collisions, mesh, X, Ut_mesh, velocities, params,
            FrictionPotential::DiffWRT::VELOCITIES, 0.0, false);

        auto F_V_aniso = [&](const Eigen::VectorXd& v) {
            return D.smooth_contact_force(
                friction_collisions, mesh, X, Ut_mesh,
                fd::unflatten(v, velocities.cols()), 0.0, false);
        };
        Eigen::MatrixXd fd_JF_wrt_V_aniso;
        fd::finite_jacobian(
            fd::flatten(velocities), F_V_aniso, fd_JF_wrt_V_aniso,
            fd::AccuracyOrder::FOURTH, 1e-6 * dhat);
        const double fd_err_V_aniso =
            (fd_JF_wrt_V_aniso - JF_wrt_V_aniso).norm();
        const bool aniso_jacobian_fd_ok =
            (JF_wrt_V_aniso.norm() == 0.0
             || fd_err_V_aniso <= 1e-7 * JF_wrt_V_aniso.norm());
        CHECK(aniso_jacobian_fd_ok);
    }
}

TEST_CASE(
    "Smooth friction force jacobian 2D", "[friction-smooth][force-jacobian]")
{
    SmoothFrictionData data = smooth_friction_data_generator_2d();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, params, barrier_stiffness] =
        data;

    Eigen::MatrixXd X, Ut, U;
    X = V0;
    Ut = V0 - X;
    U = V1 - X;

    CollisionMesh mesh(X, E, F);

    check_smooth_friction_force_jacobian(
        mesh, Ut, U, collisions, mu, epsv_times_h, params, barrier_stiffness,
        false);
}

TEST_CASE(
    "Smooth friction force no_mu and no_contact_force_multiplier",
    "[friction-smooth][force][no-mu]")
{
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(NDEBUG)
    SKIP(
        "'Smooth friction force no_mu and no_contact_force_multiplier' test is "
        "skipped in debug mode");
#endif

    SmoothFrictionData data = smooth_friction_data_generator_3d();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, params, barrier_stiffness] =
        data;

    Eigen::MatrixXd X = V0;
    Eigen::MatrixXd Ut = V0 - X;
    Eigen::MatrixXd U = V1 - X;
    CollisionMesh mesh(X, E, F);

    TangentialCollisions friction_collisions;
    friction_collisions.build(
        mesh, X + Ut, collisions, params, barrier_stiffness,
        Eigen::VectorXd::Ones(mesh.num_vertices()) * mu,
        Eigen::VectorXd::Ones(mesh.num_vertices()) * mu);

    if (friction_collisions.empty()) {
        return;
    }

    Eigen::MatrixXd velocities = U - Ut;
    friction_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, X, Ut, velocities);

    const FrictionPotential D(epsv_times_h);

    // Batch smooth_contact_force with no_mu=false then no_mu=true
    const Eigen::VectorXd force_default = D.smooth_contact_force(
        friction_collisions, mesh, X, Ut, velocities, 0.0, false);
    const Eigen::VectorXd force_no_mu = D.smooth_contact_force(
        friction_collisions, mesh, X, Ut, velocities, 0.0, true);

    CHECK(force_default.array().isFinite().all());
    CHECK(force_no_mu.array().isFinite().all());
    // With no_mu=true, mu is effectively 1; with no_mu=false, mu is applied
    // So magnitudes can differ (e.g. no_mu force larger when mu < 1)
    if (force_default.norm() > 1e-12) {
        CHECK(force_no_mu.norm() > 1e-12);
    }

    // Single-collision: no_contact_force_multiplier=true uses 1.0 instead of N
    const auto& collision = friction_collisions[0];
    const auto rest = collision.dof(X, E, F);
    const auto lagged = collision.dof(Ut, E, F);
    const auto vel = collision.dof(velocities, E, F);

    const Eigen::VectorXd local_force_N =
        D.smooth_contact_force(collision, rest, lagged, vel, false, false);
    const Eigen::VectorXd local_force_no_N =
        D.smooth_contact_force(collision, rest, lagged, vel, false, true);

    CHECK(local_force_N.array().isFinite().all());
    CHECK(local_force_no_N.array().isFinite().all());
    const double N = collision.normal_force_magnitude;
    if (N > 1e-10 && local_force_N.norm() > 1e-12) {
        // F_no_N = F_N / N (formula uses 1.0 instead of N)
        CHECK(
            (local_force_no_N - local_force_N / N).norm()
            <= 1e-8 * local_force_N.norm());
    }

    // Cover batch smooth_contact_force_jacobian with no_mu=true
    const Eigen::SparseMatrix<double> jac_no_mu =
        D.smooth_contact_force_jacobian(
            friction_collisions, mesh, X, Ut, velocities, params,
            FrictionPotential::DiffWRT::VELOCITIES, 0.0, true);
    CHECK(jac_no_mu.size() > 0);
}

TEST_CASE(
    "Smooth friction force jacobian 3D", "[friction-smooth][force-jacobian]")
{
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(NDEBUG)
    SKIP("'Smooth friction force jacobian 3D' test is skipped in debug mode");
#endif

    SmoothFrictionData data = smooth_friction_data_generator_3d();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, params, barrier_stiffness] =
        data;

    Eigen::MatrixXd X, Ut, U;
    X = V0;
    Ut = V0 - X;
    U = V1 - X;

    CollisionMesh mesh(X, E, F);

    check_smooth_friction_force_jacobian(
        mesh, Ut, U, collisions, mu, epsv_times_h, params, barrier_stiffness,
        false);
}