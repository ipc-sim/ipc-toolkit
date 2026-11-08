#include <tests/config.hpp>
#include <tests/dof_layout.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/ipc.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/gcp/gcp_potential.hpp>
#include <ipc/distance/line_line.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <igl/readCSV.h>

using namespace ipc;

TEST_CASE("Smooth barrier potential codim", "[gcp_potential]")
{
    const auto method = make_default_broad_phase();
    double dhat = 2;
    std::string mesh_name;

    Eigen::MatrixXd vertices(4, 2);
    Eigen::MatrixXi edges(2, 2), faces;

    vertices << -1, 0, 0, 0, 1, 0, 1.5, 0.2;
    edges << 0, 1, 1, 2;

    CollisionMesh mesh;

    GCPCollisions collisions;
    mesh = CollisionMesh(
        std::vector<bool>(vertices.rows(), true),
        std::vector<bool>(vertices.rows(), false), vertices, edges, faces);

    GCPParameters params(dhat, 0.85, 0.5, 0.95, 0.6, 2);
    collisions.build(mesh, vertices, params, false, method.get());
    CAPTURE(dhat, method);
    CHECK(!collisions.empty());
    CHECK(!has_intersections(mesh, vertices));

    GCPPotential potential(params);
    std::cout << "energy: " << potential(collisions, mesh, vertices) << "\n";

    // -------------------------------------------------------------------------
    // Minimum distance
    // -------------------------------------------------------------------------

    CHECK(
        collisions.compute_minimum_distance(mesh, vertices)
        <= collisions.compute_active_minimum_distance(mesh, vertices)
            * (1. + 1e-15));

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
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(
            tests::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    // REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error "
              << (grad_b - fgrad_b).norm() / grad_b.norm() << ", norms "
              << grad_b.norm() << " " << fgrad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() / grad_b.norm() < 1e-5);

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    Eigen::MatrixXd hess_b = potential.hessian(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_jacobian(
            tests::flatten(vertices), f, fhess_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    REQUIRE(hess_b.squaredNorm() > 1e-8);
    std::cout << "hess relative error "
              << (hess_b - fhess_b).norm() / hess_b.norm() << ", norms "
              << hess_b.norm() << " " << fhess_b.norm() << "\n";
    CHECK((hess_b - fhess_b).norm() / hess_b.norm() < 1e-5);
}

TEST_CASE(
    "GCP rest-shape measure keeps quadrature weights constant",
    "[gcp_potential][rest_shape_measure]")
{
    const auto method = make_default_broad_phase();
    const bool orientable = GENERATE(true, false);

    double dhat = -1;
    double perturbation = 0;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool all_vertices_on_surface = true;

    SECTION("2D edges and a codimensional point")
    {
        dhat = 2;
        perturbation = 5e-2;
        vertices.resize(4, 2);
        vertices << -1, 0, 0, 0, 1, 0, 1.5, 0.2;
        edges.resize(2, 2);
        edges << 0, 1, 1, 2;
    }
    SECTION("3D two cubes")
    {
        dhat = 1e-1;
        perturbation = 1e-3;
        all_vertices_on_surface = false;
        REQUIRE(
            tests::load_mesh("two-cubes-close.ply", vertices, edges, faces));
    }

    const std::vector<bool> is_on_surface = all_vertices_on_surface
        ? std::vector<bool>(vertices.rows(), true)
        : CollisionMesh::construct_is_on_surface(vertices.rows(), edges);
    const CollisionMesh mesh(
        is_on_surface, std::vector<bool>(vertices.rows(), orientable), vertices,
        edges, faces);

    // The mesh was constructed from the rest shape, so we deform it a little
    // to make the rest and deformed measures differ.
    const Eigen::MatrixXd rest = mesh.rest_positions();
    const Eigen::MatrixXd deformed =
        rest + perturbation * Eigen::MatrixXd::Random(rest.rows(), rest.cols());

    auto make_params = [&](const bool use_rest_shape_measure) {
        GCPParameters params(dhat, 0.85, 0.5, 0.95, 0.6, 2);
        params.use_rest_shape_measure = use_rest_shape_measure;
        return params;
    };
    auto energy = [&](const bool use_rest_shape_measure,
                      Eigen::ConstRef<Eigen::MatrixXd> V) {
        const GCPParameters params = make_params(use_rest_shape_measure);
        GCPCollisions collisions;
        collisions.build(mesh, V, params, false, method.get());
        REQUIRE(!collisions.empty());
        return GCPPotential(params)(collisions, mesh, V);
    };

    // The rest-shape measure only changes the quadrature weights, so the two
    // formulations have to agree when the mesh is in its rest configuration.
    const double rest_energy = energy(true, rest);
    CHECK(rest_energy > 0);
    CHECK(rest_energy == Catch::Approx(energy(false, rest)).epsilon(1e-10));

    // Once the mesh deforms, the deformed measure changes but the rest-shape
    // measure does not, so the energies must differ.
    CHECK(energy(true, deformed) != Catch::Approx(energy(false, deformed)));

    // The rest-shape energy must be consistently differentiated: since the
    // measure is constant, the analytic derivatives only come from the
    // smooth contact terms, and have to match finite differences of the energy.
    const GCPParameters params = make_params(true);
    GCPCollisions collisions;
    collisions.build(mesh, deformed, params, false, method.get());
    REQUIRE(!collisions.empty());
    const GCPPotential potential(params);

    const Eigen::VectorXd grad = potential.gradient(collisions, mesh, deformed);
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        tests::flatten(deformed),
        [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, tests::unflatten(x, deformed.cols()));
        },
        fgrad, fd::AccuracyOrder::SECOND, 1e-8);
    REQUIRE(grad.norm() > 0);
    CHECK((grad - fgrad).norm() / grad.norm() < 1e-5);

    const Eigen::MatrixXd hess = potential.hessian(collisions, mesh, deformed);
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        tests::flatten(deformed),
        [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, tests::unflatten(x, deformed.cols()));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-8);
    REQUIRE(hess.norm() > 0);
    CHECK((hess - fhess).norm() / hess.norm() < 1e-5);
}

#if defined(NDEBUG) && !defined(WIN32)
static std::string tagsopt = "[gcp_potential]";
#else
static std::string tagsopt = "[.][gcp_potential]";
#endif

TEST_CASE("Smooth barrier potential full gradient and hessian 3D", tagsopt)
{
    const auto method = make_default_broad_phase();
    const bool adaptive_dhat = GENERATE(true, false);
    const bool orientable = GENERATE(true, false);
    double dhat = -1;
    std::string mesh_name;
    bool all_vertices_on_surface = true;

    SECTION("two cubes far")
    {
        dhat = 1;
        mesh_name = "two-cubes-far.ply";
        all_vertices_on_surface = false;
    }
    SECTION("two cubes close")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-close.ply";
        all_vertices_on_surface = false;
    }

    double min_dist_ratio = 1.5;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    vertices +=
        Eigen::MatrixXd::Random(vertices.rows(), vertices.cols()) * 1e-3;
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh;

    GCPCollisions collisions;
    if (all_vertices_on_surface) {
        mesh = CollisionMesh(
            std::vector<bool>(vertices.rows(), true),
            std::vector<bool>(vertices.rows(), orientable), vertices, edges,
            faces);
    } else {
        mesh = CollisionMesh(
            ipc::CollisionMesh::construct_is_on_surface(vertices.rows(), edges),
            std::vector<bool>(vertices.rows(), orientable), vertices, edges,
            faces);

        vertices = mesh.vertices(vertices);
    }

    GCPParameters params(dhat, 0.85, 0.5, 0.95, 0.6, 2);
    params.set_adaptive_dhat_ratio(min_dist_ratio);
    collisions.compute_adaptive_dhat(mesh, vertices, params, method.get());
    collisions.build(mesh, vertices, params, adaptive_dhat, method.get());
    CAPTURE(dhat, method, adaptive_dhat, all_vertices_on_surface);
    CHECK(!collisions.empty());
    CHECK(!has_intersections(mesh, vertices));

    GCPPotential potential(params);
    std::cout << "energy: " << potential(collisions, mesh, vertices) << "\n";

    // -------------------------------------------------------------------------
    // Minimum distance
    // -------------------------------------------------------------------------

    CHECK(
        collisions.compute_minimum_distance(mesh, vertices)
        <= collisions.compute_active_minimum_distance(mesh, vertices)
            * (1. + 1e-15));

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
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(
            tests::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    // REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error "
              << (grad_b - fgrad_b).norm() / grad_b.norm() << ", norms "
              << grad_b.norm() << " " << fgrad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() / grad_b.norm() < 1e-5);

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    Eigen::MatrixXd hess_b = potential.hessian(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_jacobian(
            tests::flatten(vertices), f, fhess_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    REQUIRE(hess_b.squaredNorm() > 1e-8);
    std::cout << "hess relative error "
              << (hess_b - fhess_b).norm() / hess_b.norm() << ", norms "
              << hess_b.norm() << " " << fhess_b.norm() << "\n";
    CHECK((hess_b - fhess_b).norm() / hess_b.norm() < 1e-5);
}

TEST_CASE("Smooth barrier potential real sim 2D C^2", "[gcp_potential]")
{
    const auto method = make_default_broad_phase();
    const bool adaptive_dhat = GENERATE(true, false);
    const bool orientable = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name;
    SECTION("debug1")
    {
        mesh_name =
            (tests::DATA_DIR / "gcp" / "nonlinear_solve_iter020.obj").string();
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
    GCPParameters params(dhat, 0.9, -0.05, 0.95, 0.05, 1);
    params.set_adaptive_dhat_ratio(min_dist_ratio);
    GCPCollisions collisions;
    mesh = CollisionMesh(
        std::vector<bool>(vertices.rows(), true),
        std::vector<bool>(vertices.rows(), orientable), vertices, edges, faces);
    collisions.compute_adaptive_dhat(mesh, vertices, params, method.get());
    collisions.build(mesh, vertices, params, adaptive_dhat, method.get());
    CAPTURE(dhat, method, adaptive_dhat);
    CHECK(!collisions.empty());
    std::cout << "smooth collision candidate size " << collisions.size()
              << "\n";

    CHECK(!has_intersections(mesh, vertices));

    GCPPotential potential(params);
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
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(
            tests::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error "
              << (grad_b - fgrad_b).norm() / grad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() < 1e-7 * grad_b.norm());
    // CHECK(fd::compare_gradient(grad_b, fgrad_b));

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    Eigen::MatrixXd hess_b = potential.hessian(collisions, mesh, vertices);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_b;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_jacobian(
            tests::flatten(vertices), f, fhess_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    REQUIRE(hess_b.squaredNorm() > 1e-3);
    std::cout << "hess relative error "
              << (hess_b - fhess_b).norm() / hess_b.norm() << "\n";
    CHECK((hess_b - fhess_b).norm() < 1e-7 * hess_b.norm());
    // CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
}

TEST_CASE("Smooth barrier potential real sim 2D C^1", "[gcp_potential]")
{
    const auto method = make_default_broad_phase();
    const bool adaptive_dhat = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name;
    SECTION("debug2")
    {
        mesh_name = (tests::DATA_DIR / "gcp" / "simple_2d.obj").string();
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
    GCPParameters params(dhat, 0.9, -0.05, 0.95, 0.05, 1);
    params.set_adaptive_dhat_ratio(min_dist_ratio);
    GCPCollisions collisions;
    mesh = CollisionMesh(vertices, edges, faces);
    collisions.compute_adaptive_dhat(mesh, vertices, params, method.get());
    collisions.build(mesh, vertices, params, adaptive_dhat, method.get());
    CAPTURE(dhat, method, adaptive_dhat);
    CHECK(!collisions.empty());
    std::cout << "smooth collision candidate size " << collisions.size()
              << "\n";
    std::cout << collisions.to_string(mesh, vertices, params) << "\n";

    CHECK(!has_intersections(mesh, vertices));

    GCPPotential potential(params);
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
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        };
        fd::finite_gradient(
            tests::flatten(vertices), f, fgrad_b, fd::AccuracyOrder::SECOND,
            1e-8);
    }

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    std::cout << "grad relative error "
              << (grad_b - fgrad_b).norm() / grad_b.norm() << "\n";
    CHECK((grad_b - fgrad_b).norm() < 1e-7 * grad_b.norm());
    // CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

TEST_CASE("Benchmark autogen code", "[!benchmark]")
{
    Eigen::Vector3d ea0, ea1, eb0, eb1;
    ea0 << -0.9, 0, 0;
    ea1 << 1.05, 0, 0;
    eb0 << 0, -1.1, 1.02;
    eb1 << 0, 0.99, 1.01;

    BENCHMARK("autogen")
    {
        ipc::line_line_distance_hessian(ea0, ea1, eb0, eb1);
    };

    BENCHMARK("manual")
    {
        ipc::line_line_closest_point_direction_hessian(ea0, ea1, eb0, eb1);
    };
}
