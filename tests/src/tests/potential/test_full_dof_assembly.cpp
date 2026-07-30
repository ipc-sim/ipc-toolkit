// Tests for assembling potential derivatives directly in full-mesh DOF
// (`in_full_dof=true`), which must match applying `CollisionMesh::to_full_dof`
// to the collision-DOF result.

#include "assembly_scene.hpp"

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/potentials/barrier_potential.hpp>

using namespace ipc;

TEST_CASE(
    "Full-DOF assembly matches to_full_dof", "[potential][assembly][full_dof]")
{
    // Scenes are padded with interior vertices, so the selection matrix is a
    // genuine subset selection (full_ndof == 2 * ndof), not a permutation.
#ifdef NDEBUG
    const size_t scene_index = GENERATE(size_t(0), size_t(1), size_t(3));
#else
    const size_t scene_index = 0; // two-cubes; larger scenes are slow in debug
#endif
    const auto& spec = ipc::tests::assembly_scene_specs().at(scene_index);

    const std::optional<ipc::tests::AssemblyScene> maybe_scene =
        ipc::tests::build_assembly_scene(spec);
    if (!maybe_scene.has_value()) {
        SKIP(fmt::format("Scene '{}' is unavailable.", spec.mesh_name));
    }
    const ipc::tests::AssemblyScene& scene = maybe_scene.value();
    CAPTURE(scene.label());

    const CollisionMesh& mesh = scene.mesh();
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXd& X = scene.vertices();

    REQUIRE(mesh.is_selection_dof_map());
    REQUIRE(mesh.full_ndof() > mesh.ndof());

    SECTION("gradient")
    {
        const Eigen::VectorXd grad_folded =
            potential.gradient(collisions, mesh, X, /*in_full_dof=*/true);
        const Eigen::VectorXd grad_mapped =
            mesh.to_full_dof(potential.gradient(collisions, mesh, X));

        REQUIRE(grad_folded.size() == mesh.full_ndof());
        // tbb::combinable partitions work nondeterministically, so the two
        // calls may sum in different orders; compare with a tight tolerance.
        const double scale = std::max(1.0, grad_mapped.norm());
        CHECK((grad_folded - grad_mapped).norm() <= 1e-13 * scale);
    }

    SECTION("hessian")
    {
        const PSDProjectionMethod psd =
            GENERATE(PSDProjectionMethod::NONE, PSDProjectionMethod::CLAMP);
        CAPTURE(psd);

        const Eigen::SparseMatrix<double> hess_folded =
            potential.hessian(collisions, mesh, X, psd, /*in_full_dof=*/true);
        const Eigen::SparseMatrix<double> hess_mapped =
            mesh.to_full_dof(potential.hessian(collisions, mesh, X, psd));

        REQUIRE(hess_folded.rows() == mesh.full_ndof());
        REQUIRE(hess_folded.cols() == mesh.full_ndof());

        const double scale = std::max(1.0, hess_mapped.norm());
        CHECK((hess_folded - hess_mapped).norm() <= 1e-13 * scale);
    }
}

TEST_CASE(
    "Full-DOF assembly with a non-selection displacement map",
    "[potential][assembly][full_dof]")
{
    // With a user-provided displacement map, folding by index remap is not
    // possible; in_full_dof must fall back to applying to_full_dof and still
    // produce the mapped result.
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("two-cubes-close.ply", vertices, edges, faces));

    // A diagonal (non-identity) displacement map: full vars are half-scale.
    Eigen::SparseMatrix<double> displacement_map(
        vertices.rows(), vertices.rows());
    displacement_map.setIdentity();
    displacement_map *= 0.5;

    const CollisionMesh mesh(vertices, edges, faces, displacement_map);
    REQUIRE(!mesh.is_selection_dof_map());

    const double dhat = 1e-1;
    NormalCollisions collisions;
    collisions.build(mesh, vertices, dhat);
    REQUIRE(!collisions.empty());

    const BarrierPotential potential(dhat, /*stiffness=*/1.0);

    // Both calls run the same collision-DOF assembly (whose tbb reduction is
    // order-nondeterministic) followed by to_full_dof, so compare with a
    // tight tolerance rather than exactly.
    const Eigen::VectorXd grad_full =
        potential.gradient(collisions, mesh, vertices, /*in_full_dof=*/true);
    const Eigen::VectorXd grad_mapped =
        mesh.to_full_dof(potential.gradient(collisions, mesh, vertices));
    CHECK(
        (grad_full - grad_mapped).norm()
        <= 1e-13 * std::max(1.0, grad_mapped.norm()));

    const Eigen::SparseMatrix<double> hess_full = potential.hessian(
        collisions, mesh, vertices, PSDProjectionMethod::NONE,
        /*in_full_dof=*/true);
    const Eigen::SparseMatrix<double> hess_mapped =
        mesh.to_full_dof(potential.hessian(collisions, mesh, vertices));
    CHECK(
        (hess_full - hess_mapped).norm()
        <= 1e-13 * std::max(1.0, hess_mapped.norm()));
}

TEST_CASE(
    "Full-DOF assembly with no collisions", "[potential][assembly][full_dof]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("two-cubes-far.ply", vertices, edges, faces));

    const CollisionMesh mesh =
        CollisionMesh::build_from_full_mesh(vertices, edges, faces);
    vertices = mesh.vertices(vertices);

    const NormalCollisions collisions; // empty
    const BarrierPotential potential(1e-3, /*stiffness=*/1.0);

    const Eigen::VectorXd grad =
        potential.gradient(collisions, mesh, vertices, /*in_full_dof=*/true);
    CHECK(grad.size() == mesh.full_ndof());
    CHECK(grad.norm() == 0.0);

    const Eigen::SparseMatrix<double> hess = potential.hessian(
        collisions, mesh, vertices, PSDProjectionMethod::NONE,
        /*in_full_dof=*/true);
    CHECK(hess.rows() == mesh.full_ndof());
    CHECK(hess.cols() == mesh.full_ndof());
    CHECK(hess.nonZeros() == 0);
}
