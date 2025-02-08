#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>

#include <igl/edges.h>

using namespace ipc;

TEST_CASE("Initial barrier stiffness", "[stiffness][adaptive_stiffness]")
{
    const double bbox_diagonal = 1.0;
    const ClampedLogBarrier barrier;
    const double dhat = 1e-3;
    const double average_mass = 1.0;
    const Eigen::VectorXd grad_energy = Eigen::VectorXd::Constant(1, 100);
    const Eigen::VectorXd grad_barrier = Eigen::VectorXd::Constant(1, -100);

    double max_barrier_stiffness;
    const double kappa = initial_barrier_stiffness(
        bbox_diagonal, barrier, dhat, average_mass, grad_energy, grad_barrier,
        max_barrier_stiffness);

    const double expected_kappa = 1e11 * average_mass
        / (4 * 1e-16 * barrier_second_derivative(1e-16, dhat * dhat));

    CHECK(kappa == Catch::Approx(expected_kappa));
    CHECK(max_barrier_stiffness == Catch::Approx(100 * expected_kappa));
}

TEST_CASE("Update barrier stiffness", "[stiffness][adaptive_stiffness]")
{
    double prev_min_distance = 1e-10 * 1e-10;
    double min_distance = 1e-11 * 1e-11;
    double max_barrier_stiffness = 10.0;
    double barrier_stiffness = 1.0;
    double bbox_diagonal = 1.0;

    double kappa = update_barrier_stiffness(
        prev_min_distance, min_distance, max_barrier_stiffness,
        barrier_stiffness, bbox_diagonal);

    REQUIRE(kappa == Catch::Approx(2.0));
}

TEST_CASE("Semi-implicit stiffness", "[stiffness]")
{
    Eigen::MatrixXd vertices(4, 3);
    Eigen::MatrixXi edges, faces(1, 3);

    const double d = GENERATE(range(0.1, 1.0, 0.1));
    vertices << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.333, 0.333, d;
    faces << 0, 1, 2;
    igl::edges(faces, edges);

    const CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions collisions;
    collisions.fv_collisions.emplace_back(0, 3);

    const double m = GENERATE(range(1.0, 10.0, 1.0));
    const Eigen::VectorXd vertex_masses =
        Eigen::VectorXd::Constant(vertices.rows(), m);

    const double dmin = 0;

    Eigen::SparseMatrix<double> hess(vertices.size(), vertices.size()); // zero

    const Eigen::VectorXd kappa = semi_implicit_stiffness(
        mesh, vertices, collisions, vertex_masses, hess, dmin);

    REQUIRE(kappa.size() == 1);
    REQUIRE(kappa(0) == Catch::Approx(m / (d * d)));
}