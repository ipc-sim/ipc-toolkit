#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <tests/utils.hpp>

#include <ipc/ipc.hpp>
#include <ipc/adhesion/adhesion.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/normal_adhesion_potential.hpp>
#include <ipc/potentials/tangential_adhesion_potential.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

using namespace ipc;

TEST_CASE("Normal adhesion potential", "[potential][adhesion]")
{
    const double dhat_p = 1e-3;
    const double dhat_a = 2 * dhat_p;
    const double Y = 1e3;
    const double eps_c = 0.5;

    // Setup a mesh

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    tests::load_mesh("cube.ply", vertices, edges, faces);

    const double gap =
        GENERATE_COPY(0.1, 0.5, 0.999, 1.0, 1.001, 1.5, 1.999) * dhat_p;

    const size_t n = vertices.rows();
    vertices.conservativeResize(2 * n, Eigen::NoChange);
    vertices.bottomRows(n) = vertices.topRows(n);
    vertices.bottomRows(n).col(1).array() += 1 + gap;
    bool rotated = GENERATE(false, true);
    if (rotated) {
        Eigen::Matrix3d R =
            Eigen::AngleAxis<double>(
                Eigen::AngleAxisd(igl::PI / 4, Eigen::Vector3d::UnitY()))
                .toRotationMatrix();
        vertices.bottomRows(n) =
            (vertices.bottomRows(n) * R.transpose()).eval();
    }

    edges.conservativeResize(2 * edges.rows(), Eigen::NoChange);
    edges.bottomRows(edges.rows() / 2) =
        edges.topRows(edges.rows() / 2).array() + n;

    faces.conservativeResize(2 * faces.rows(), Eigen::NoChange);
    faces.bottomRows(faces.rows() / 2) =
        faces.topRows(faces.rows() / 2).array() + n;

    CollisionMesh mesh(vertices, edges, faces);

    REQUIRE(!has_intersections(mesh, vertices));

    NormalCollisions collisions;
    const bool use_area_weighting = GENERATE(false, true);
    // TODO: Debug why the improved max approx. requires area weighting.
    const bool use_improved_max_approximator = use_area_weighting;
    collisions.set_use_area_weighting(use_area_weighting);
    collisions.set_use_improved_max_approximator(use_improved_max_approximator);
    collisions.build(mesh, vertices, dhat_a);

    REQUIRE(collisions.size() > 0);

    REQUIRE(
        collisions.compute_minimum_distance(mesh, vertices)
        == Catch::Approx(gap * gap));

    CAPTURE(
        use_area_weighting, use_improved_max_approximator, dhat_a, dhat_p, Y,
        eps_c, gap);

    // --- Check the potential gradients ---------------------------------------

    CHECK(BarrierPotential(dhat_a)(collisions, mesh, vertices) > 0);

    NormalAdhesionPotential potential(dhat_p, dhat_a, Y, eps_c);
    CHECK(potential(collisions, mesh, vertices) < 0);
    if (potential(collisions, mesh, vertices) >= 0) {
        logger().critical(collisions.to_string(mesh, vertices));
    }

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        fd::flatten(vertices),
        [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        },
        fgrad);

    Eigen::VectorXd grad = potential.gradient(collisions, mesh, vertices);
    CHECK(fd::compare_gradient(grad, fgrad));

    // -------------------------------------------------------------------------

    if (gap != dhat_p) { // The potential is only C¹ at d = d̂ₚ
        Eigen::MatrixXd hess = potential.hessian(collisions, mesh, vertices);

        Eigen::MatrixXd fhess;

        fd::finite_jacobian(
            fd::flatten(vertices),
            [&](const Eigen::VectorXd& x) {
                return potential.gradient(
                    collisions, mesh, fd::unflatten(x, vertices.cols()));
            },
            fhess);

        REQUIRE(hess.squaredNorm() > 1e-3);
        CHECK(fd::compare_hessian(hess, fhess, 1e-3));
    }

    // --- Maximum normal adhesion force magnitude -----------------------------

    CHECK(
        potential.force_magnitude(1e-3, 0, 1)
        == max_normal_adhesion_force_magnitude(
            dhat_p * dhat_p, dhat_a * dhat_a,
            Y * eps_c / (4 * (dhat_p) * (dhat_p * dhat_p - dhat_a * dhat_a))));

    CHECK(potential.force_magnitude_gradient(1e-3, VectorMax12d::Zero(12), 0, 1)
              .isZero());
}

TEST_CASE("Tangetial adhesion potential", "[potential][adhesion]")
{
    TangentialAdhesionPotential potential(1e-3);
}