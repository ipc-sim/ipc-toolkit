#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/potentials/barrier_potential.hpp>

using namespace ipc;

TEST_CASE("Barrier Potential Refactor", "[potential][barrier_potential]")
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

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    REQUIRE(success);

    const CollisionMesh mesh(vertices, edges, faces);

    CollisionConstraints contacts;
    contacts.set_use_convergent_formulation(use_convergent_formulation);
    contacts.build(mesh, vertices, dhat);
    // contacts.ee_constraints.clear(); // Remove mollified collisions

    BarrierPotential B(dhat);

    CAPTURE(mesh_name, dhat);
    CHECK(contacts.size() > 0);

    const double expected_potential =
        contacts.compute_potential(mesh, vertices, dhat);
    const double actual_potential = B(mesh, vertices, contacts);
    CHECK(actual_potential == Catch::Approx(expected_potential));

    const Eigen::VectorXd expected_gradient =
        contacts.compute_potential_gradient(mesh, vertices, dhat);
    const Eigen::VectorXd actual_gradient =
        B.gradient(mesh, vertices, contacts);
    CHECK(actual_gradient.isApprox(expected_gradient));

    Eigen::SparseMatrix<double> expected_hessian, actual_hessian;
    expected_hessian =
        contacts.compute_potential_hessian(mesh, vertices, dhat, false);
    actual_hessian = B.hessian(mesh, vertices, contacts, false);
    CHECK(actual_hessian.isApprox(expected_hessian));

    // Projected hessian
    expected_hessian =
        contacts.compute_potential_hessian(mesh, vertices, dhat, true);
    actual_hessian = B.hessian(mesh, vertices, contacts, true);
    INFO("project_hessian_to_psd = true");
    CHECK(actual_hessian.isApprox(expected_hessian));
}
