#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/smooth_contact/smooth_contact_potential.hpp>
#include <ipc/distance/line_line.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <igl/readCSV.h>
#include <ipc/ipc.hpp>

using namespace ipc;

#if defined(NDEBUG) && !defined(WIN32)
std::string tagsopt = "[smooth_potential]";
#else
std::string tagsopt = "[.][smooth_potential]";
#endif

TEST_CASE(
    "Number of contact pairs",
    tagsopt)
{
    const BroadPhaseMethod method{0};

    double dhat = -1;
    std::string mesh_name = "";
    bool all_vertices_on_surface = true;

    SECTION("two cubes close")
    {
        dhat = 5e-3;
        mesh_name = (tests::DATA_DIR / "step_1000_surf_contact.obj").string();
        all_vertices_on_surface = true;
    }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(vertices, edges, faces);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(vertices, edges, faces);
        vertices = mesh.vertices(vertices);
    }

    {
        Collisions collisions;

        collisions.build(mesh, vertices, dhat, /*dmin=*/0, method);
        CHECK(collisions.size() > 0);
        std::cout << "IPC number of pairs " << collisions.size() << "\n";
    }

    {
        Collisions collisions;
        collisions.set_use_convergent_formulation(true);

        collisions.build(mesh, vertices, dhat, /*dmin=*/0, method);
        CHECK(collisions.size() > 0);
        std::cout << "CIPC number of pairs " << collisions.size() << "\n";
    }

    {
        SmoothCollisions<3> collisions;

        ParameterType param(dhat, 0.8, 0, 1, 0, 2);
        collisions.build(mesh, vertices, param, false, method);
        CHECK(collisions.size() > 0);
        std::cout << "OIPC number of pairs (only tangent) " << collisions.size() << "\n";
    }

    {
        SmoothCollisions<3> collisions;

        ParameterType param(dhat, 1, 0, 0, 0.1, 2);
        collisions.build(mesh, vertices, param, false, method);
        CHECK(collisions.size() > 0);
        std::cout << "OIPC number of pairs (only normal) " << collisions.size() << "\n";
    }

    {
        SmoothCollisions<3> collisions;

        ParameterType param(dhat, 0.8, 0, 0, 0.1, 2);
        collisions.build(mesh, vertices, param, false, method);
        CHECK(collisions.size() > 0);
        std::cout << "OIPC number of pairs (both) " << collisions.size() << "\n";
    }

    {
        SmoothCollisions<3> collisions;

        ParameterType param(dhat, 0.5, 0, 0, 0.1, 2);
        collisions.build(mesh, vertices, param, false, method);
        CHECK(collisions.size() > 0);
        std::cout << "OIPC number of pairs (both) " << collisions.size() << "\n";
    }

    {
        SmoothCollisions<3> collisions;

        ParameterType param(dhat, 0.1, 0, 0, 0.1, 2);
        collisions.build(mesh, vertices, param, false, method);
        CHECK(collisions.size() > 0);
        std::cout << "OIPC number of pairs (both) " << collisions.size() << "\n";
    }
}

TEST_CASE(
    "Smooth barrier potential full gradient and hessian 3D",
    tagsopt)
{
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID;
    const bool adaptive_dhat = GENERATE(true, false);
    double dhat = -1;
    std::string mesh_name = "";
    bool all_vertices_on_surface = true;
    // energy is zero
    // SECTION("cube")
    // {
    //     dhat = sqrt(2.0);
    //     mesh_name = "cube.obj";
    // }
    SECTION("two cubes far")
    {
        dhat = 1;
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

    double min_dist_ratio = 1.5;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    vertices += Eigen::MatrixXd::Random(vertices.rows(), vertices.cols()) * 1e-3;
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;

    SmoothCollisions<3> collisions;
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(vertices, edges, faces);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(vertices, edges, faces);
        vertices = mesh.vertices(vertices);
    }

    ParameterType param(dhat, 0.85, 0.5, 0.95, 0.6, 2);
    param.set_adaptive_dhat_ratio(min_dist_ratio);
    collisions.compute_adaptive_dhat(mesh, vertices, param, method);
    collisions.build(mesh, vertices, param, adaptive_dhat, method);
    CAPTURE(dhat, method, adaptive_dhat, all_vertices_on_surface);
    CHECK(collisions.size() > 0);
    CHECK(!has_intersections(mesh, vertices));

    SmoothContactPotential<SmoothCollisions<3>> potential(param);
    std::cout << "energy: " << potential(collisions, mesh, vertices) << "\n";

    // -------------------------------------------------------------------------
    // Gradient
    // -------------------------------------------------------------------------

    const Eigen::VectorXd grad_b =
        potential.gradient(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(fd::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND, 1e-8);
    }

    // REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error " << (grad_b - fgrad_b).norm() / grad_b.norm() << ", norms " << grad_b.norm() << " " << fgrad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() / grad_b.norm() < 1e-5);

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    Eigen::MatrixXd hess_b =
        potential.hessian(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_jacobian(fd::flatten(vertices), f, fhess_b, fd::AccuracyOrder::SECOND, 1e-8);
    }

    REQUIRE(hess_b.squaredNorm() > 1e-8);
    std::cout << "hess relative error " << (hess_b - fhess_b).norm() / hess_b.norm() << ", norms " << hess_b.norm() << " " << fhess_b.norm() << "\n";
    CHECK((hess_b - fhess_b).norm() / hess_b.norm() < 1e-5);
}

TEST_CASE(
    "Smooth barrier potential real sim 2D C^2",
    "[smooth_potential]")
{
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID;
    const bool adaptive_dhat = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name = "";
    SECTION("debug1")
    {
        mesh_name = (tests::DATA_DIR / "nonlinear_solve_iter020.obj").string();
        dhat = 3e-2;
    }

    double min_dist_ratio = 1.5;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = igl::readCSV(mesh_name + "-v.csv", vertices);
    success = success && igl::readCSV(mesh_name + "-e.csv", edges);
    CAPTURE(mesh_name);
    REQUIRE(success);

    // std::cout << "\n" <<  vertices << "\n" << edges << "\n";

    CollisionMesh mesh;
    ParameterType param(dhat, 0.9, -0.05, 0.95, 0.05, 1);
    param.set_adaptive_dhat_ratio(min_dist_ratio);
    SmoothCollisions<2> collisions;
    mesh = CollisionMesh(vertices, edges, faces);
    collisions.compute_adaptive_dhat(mesh, vertices, param, method);
    collisions.build(mesh, vertices, param, adaptive_dhat, method);
    CAPTURE(dhat, method, adaptive_dhat);
    CHECK(collisions.size() > 0);
    std::cout << "smooth collision candidate size " << collisions.size() << "\n";

    CHECK(!has_intersections(mesh, vertices));

    SmoothContactPotential<SmoothCollisions<2>> potential(param);
    std::cout << "energy: " << potential(collisions, mesh, vertices) << "\n";

    // -------------------------------------------------------------------------
    // Gradient
    // -------------------------------------------------------------------------

    const Eigen::VectorXd grad_b =
        potential.gradient(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(fd::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND, 1e-8);
    }

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error " << (grad_b - fgrad_b).norm() / grad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() < 1e-7 * grad_b.norm());
    // CHECK(fd::compare_gradient(grad_b, fgrad_b));

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    Eigen::MatrixXd hess_b =
        potential.hessian(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_jacobian(fd::flatten(vertices), f, fhess_b, fd::AccuracyOrder::SECOND, 1e-8);
    }

    REQUIRE(hess_b.squaredNorm() > 1e-3);
    std::cout << "hess relative error " << (hess_b - fhess_b).norm() / hess_b.norm() << "\n";
    CHECK((hess_b - fhess_b).norm() < 1e-7 * hess_b.norm());
    // CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
}

TEST_CASE(
    "Smooth barrier potential real sim 2D C^1",
    "[smooth_potential]")
{
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID;
    const bool adaptive_dhat = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name = "";
    SECTION("debug2")
    {
        mesh_name = (tests::DATA_DIR / "simple_2d.obj").string();
        dhat = 0.1;
    }

    double min_dist_ratio = 1.5;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = igl::readCSV(mesh_name + "-v.csv", vertices);
    success = success && igl::readCSV(mesh_name + "-e.csv", edges);
    CAPTURE(mesh_name);
    REQUIRE(success);

    // std::cout << "\n" <<  vertices << "\n" << edges << "\n";

    CollisionMesh mesh;
    ParameterType param(dhat, 0.9, -0.05, 0.95, 0.05, 1);
    param.set_adaptive_dhat_ratio(min_dist_ratio);
    SmoothCollisions<2> collisions;
    mesh = CollisionMesh(vertices, edges, faces);
    collisions.compute_adaptive_dhat(mesh, vertices, param, method);
    collisions.build(mesh, vertices, param, adaptive_dhat, method);
    CAPTURE(dhat, method, adaptive_dhat);
    CHECK(collisions.size() > 0);
    std::cout << "smooth collision candidate size " << collisions.size() << "\n";

    CHECK(!has_intersections(mesh, vertices));

    SmoothContactPotential<SmoothCollisions<2>> potential(param);
    std::cout << "energy: " << potential(collisions, mesh, vertices) << "\n";

    // -------------------------------------------------------------------------
    // Gradient
    // -------------------------------------------------------------------------

    const Eigen::VectorXd grad_b =
        potential.gradient(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(fd::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND, 1e-8);
    }

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error " << (grad_b - fgrad_b).norm() / grad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() < 1e-7 * grad_b.norm());
    // CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

TEST_CASE(
    "Benchmark on OIPC",
    tagsopt)
{
    const BroadPhaseMethod method{0};

    double dhat = -1;
    std::string mesh_name = "";
    bool all_vertices_on_surface = true;

    SECTION("mat-twist")
    {
        dhat = 1e-3;
        mesh_name = (tests::DATA_DIR / "step_1000_surf_contact.obj").string();
        all_vertices_on_surface = true;
    }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(vertices, edges, faces);
    } else {
        mesh = CollisionMesh::build_from_full_mesh(vertices, edges, faces);
        vertices = mesh.vertices(vertices);
    }

    {
        igl::Timer timer;
        timer.start();
        Collisions collisions;

        collisions.build(mesh, vertices, dhat, /*dmin=*/0, method);
        CHECK(collisions.size() > 0);
        // std::cout << "IPC number of pairs " << collisions.size() << "\n";

        timer.stop();
        std::cout << "IPC build time " << timer.getElapsedTime() << " s\n";
        timer.start();

        BarrierPotential barrier_potential(dhat);
        const Eigen::VectorXd grad_b =
            barrier_potential.gradient(collisions, mesh, vertices);

        timer.stop();
        std::cout << "IPC grad time " << timer.getElapsedTime() << " s\n";
        timer.start();

        barrier_potential.hessian(collisions, mesh, vertices);

        timer.stop();
        std::cout << "IPC hess time " << timer.getElapsedTime() << " s\n";
        timer.start();
    }

    {
        igl::Timer timer;
        timer.start();
        SmoothCollisions<3> collisions;

        ParameterType param(dhat, 0.8, 0, 0, 0.1, 2);
        collisions.build(mesh, vertices, param, false, method);
        CHECK(collisions.size() > 0);
        // std::cout << "OIPC number of pairs (both) " << collisions.size() << "\n";

        timer.stop();
        std::cout << "OIPC build time " << timer.getElapsedTime() << " s\n";
        timer.start();

        SmoothContactPotential<SmoothCollisions<3>> potential(param);
        const Eigen::VectorXd grad_b =
            potential.gradient(collisions, mesh, vertices);

        timer.stop();
        std::cout << "OIPC grad time " << timer.getElapsedTime() << " s\n";

        potential.hessian(collisions, mesh, vertices);

        timer.stop();
        std::cout << "OIPC hess time " << timer.getElapsedTime() << " s\n";
        timer.start();
    }
}

TEST_CASE(
    "Benchmark autogen code", "[!benchmark]")
{
    ipc::Vector3d ea0, ea1, eb0, eb1;
    ea0 << -0.9, 0, 0;
    ea1 <<  1.05, 0, 0;
    eb0 << 0, -1.1, 1.02;
    eb1 << 0,  0.99, 1.01;

    BENCHMARK("autogen")
    {
        ipc::line_line_distance_hessian(ea0, ea1, eb0, eb1);
    };

    BENCHMARK("manual")
    {
        ipc::line_line_closest_point_direction_hessian(ea0, ea1, eb0, eb1);
    };
}
