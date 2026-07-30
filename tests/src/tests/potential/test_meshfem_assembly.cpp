// Tests for the MeshFEMSparse-backed Hessian assembler: it must produce the
// same global matrix as the default triplet assembler.
//
// Only compiled to something when IPC_TOOLKIT_WITH_MESHFEM_SPARSE is enabled.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE

#include "assembly_scene.hpp"

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/utils/hessian_assembler.hpp>
#include <ipc/utils/meshfem_hessian_assembler.hpp>

#include <limits>

using namespace ipc;

TEST_CASE(
    "MeshFEM assembly matches triplet assembly",
    "[potential][assembly][meshfem]")
{
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

    const PSDProjectionMethod psd =
        GENERATE(PSDProjectionMethod::NONE, PSDProjectionMethod::CLAMP);
    const bool in_full_dof = GENERATE(false, true);
    CAPTURE(psd, in_full_dof);

    TripletHessianAssembler triplet_assembler;
    potential.assemble_hessian(
        collisions, mesh, X, triplet_assembler, psd, in_full_dof);
    const Eigen::SparseMatrix<double> expected = triplet_assembler.get_matrix();

    MeshFEMHessianAssembler meshfem_assembler;
    potential.assemble_hessian(
        collisions, mesh, X, meshfem_assembler, psd, in_full_dof);
    const Eigen::SparseMatrix<double> actual = meshfem_assembler.get_matrix();

    REQUIRE(actual.rows() == expected.rows());
    REQUIRE(actual.cols() == expected.cols());

    // Same additions in a different order; compare with a tight tolerance.
    const double scale = std::max(1.0, expected.norm());
    CHECK((actual - expected).norm() <= 1e-13 * scale);
}

TEST_CASE("MeshFEM assembly pattern reuse", "[potential][assembly][meshfem]")
{
    // A persistent assembler must produce correct results across repeated
    // assemblies, reusing the cached sparsity pattern when the contact set
    // allows it.
#ifndef NDEBUG
    SKIP("Building the bunny's collision sets is too slow in debug mode.");
#endif
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("bunny.ply", vertices, edges, faces));

    const CollisionMesh mesh =
        CollisionMesh::build_from_full_mesh(vertices, edges, faces);
    vertices = mesh.vertices(vertices);

    // Two nested contact sets: the small-dhat set is a subset of the large-
    // dhat set (fewer vertex pairs within the activation distance).
    const double dhat_large = 1e-2, dhat_small = 5e-3;
    NormalCollisions collisions_large, collisions_small;
    collisions_large.build(mesh, vertices, dhat_large);
    collisions_small.build(mesh, vertices, dhat_small);
    REQUIRE(!collisions_small.empty());
    REQUIRE(collisions_small.size() < collisions_large.size());

    const BarrierPotential potential(dhat_large, /*stiffness=*/1.0);

    const auto reference = [&](const NormalCollisions& collisions) {
        TripletHessianAssembler triplet_assembler;
        potential.assemble_hessian(
            collisions, mesh, vertices, triplet_assembler);
        return triplet_assembler.get_matrix();
    };
    const auto check_matches = [](const Eigen::SparseMatrix<double>& actual,
                                  const Eigen::SparseMatrix<double>& expected) {
        // Stale blocks assemble to explicit zeros, so compare values (the
        // difference ignores pattern mismatches), not nonZeros().
        const double scale = std::max(1.0, expected.norm());
        CHECK((actual - expected).norm() <= 1e-13 * scale);
    };

    MeshFEMHessianAssembler assembler;

    // 1. First assembly builds the pattern.
    potential.assemble_hessian(collisions_large, mesh, vertices, assembler);
    CHECK(!assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_large));

    // 2. Same collision set: the pattern must be reused.
    potential.assemble_hessian(collisions_large, mesh, vertices, assembler);
    CHECK(assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_large));

    // 3. Shrunken collision set with a permissive tolerance: reused pattern
    //    with stale (explicitly zero) blocks; values must still be correct.
    assembler.set_stale_block_tolerance(std::numeric_limits<int>::max());
    potential.assemble_hessian(collisions_small, mesh, vertices, assembler);
    CHECK(assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_small));

    // 4. Shrunken collision set with zero tolerance: rebuild.
    assembler.set_stale_block_tolerance(0);
    potential.assemble_hessian(collisions_small, mesh, vertices, assembler);
    CHECK(!assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_small));

    // 5. Grown collision set: new entries always force a rebuild, regardless
    //    of the tolerance.
    assembler.set_stale_block_tolerance(std::numeric_limits<int>::max());
    potential.assemble_hessian(collisions_large, mesh, vertices, assembler);
    CHECK(!assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_large));

    // 6. Assume-unchanged fast path: identical set, detection skipped.
    assembler.set_stale_block_tolerance(0);
    assembler.set_assume_unchanged_stencils(true);
    potential.assemble_hessian(collisions_large, mesh, vertices, assembler);
    CHECK(assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_large));

    // 7. Assume-unchanged with a differing stencil count: the assumption is
    //    disproven, so it falls back to detection (and rebuilds here, since
    //    the tolerance is zero and blocks disappeared).
    potential.assemble_hessian(collisions_small, mesh, vertices, assembler);
    CHECK(!assembler.reused_pattern());
    check_matches(assembler.get_matrix(), reference(collisions_small));
}

TEST_CASE(
    "MeshFEM assembly with no collisions", "[potential][assembly][meshfem]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("two-cubes-far.ply", vertices, edges, faces));

    const CollisionMesh mesh =
        CollisionMesh::build_from_full_mesh(vertices, edges, faces);
    vertices = mesh.vertices(vertices);

    const NormalCollisions collisions; // empty
    const BarrierPotential potential(1e-3, /*stiffness=*/1.0);

    MeshFEMHessianAssembler assembler;
    potential.assemble_hessian(collisions, mesh, vertices, assembler);
    const Eigen::SparseMatrix<double> hess = assembler.get_matrix();

    CHECK(hess.rows() == mesh.ndof());
    CHECK(hess.cols() == mesh.ndof());
    CHECK(hess.norm() == 0.0);
}

#endif // IPC_TOOLKIT_WITH_MESHFEM_SPARSE
