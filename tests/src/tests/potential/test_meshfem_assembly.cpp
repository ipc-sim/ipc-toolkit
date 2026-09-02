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
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ipc/utils/hessian_assembler.hpp>
#include <ipc/utils/meshfem_hessian_assembler.hpp>

// Deliberately included here rather than by meshfem_hessian_assembler.hpp:
// block_matrix() is declared against a forward declaration, so only callers
// that want the block matrix pay for MeshFEM's headers.
#include <MeshFEMSparse/BlockCSCHessian.hh>

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

TEST_CASE(
    "MeshFEM assembly with a degenerate stencil",
    "[potential][assembly][meshfem]")
{
    // A stencil can repeat a vertex -- e.g., a point coincident with an
    // endpoint of the edge it collides with, which the friction test fixtures
    // exercise. Several local blocks then map to one global block, which the
    // column-merge scatter cannot express, so the assembler falls back. Drive
    // the assemblers directly: the geometry that produces such a stencil is
    // degenerate, but the assembly of it must still match triplets exactly.
    const int dim = GENERATE(2, 3);
    CAPTURE(dim);

    const int ndof = 4 * dim; // four vertices
    // Edge-vertex stencil {vertex, e0, e1} where the vertex is endpoint e0.
    const std::array<index_t, 4> ids = { { 1, 1, 2, -1 } };
    const int n = 3 * dim;

    MatrixMax12d local_hess(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            local_hess(i, j) = (i + 1.0) * (j + 2.0);
        }
    }
    // Local Hessians are symmetric; the block-CSC backend stores only the
    // upper triangle, so an asymmetric input would not be a fair comparison.
    local_hess = (0.5 * (local_hess + local_hess.transpose())).eval();
    local_hess.diagonal().array() += 10.0;

    const auto stencil = [&ids](const size_t) { return ids; };

    TripletHessianAssembler triplet_assembler;
    triplet_assembler.begin(ndof, dim, 1, stencil);
    triplet_assembler.add_local_hessian(local_hess, ids);
    triplet_assembler.end();
    const Eigen::SparseMatrix<double> expected = triplet_assembler.get_matrix();

    MeshFEMHessianAssembler meshfem_assembler;
    meshfem_assembler.begin(ndof, dim, 1, stencil);
    meshfem_assembler.add_local_hessian(local_hess, ids);
    meshfem_assembler.end();
    const Eigen::SparseMatrix<double> actual = meshfem_assembler.get_matrix();

    REQUIRE(expected.norm() > 0); // guard against comparing two empty matrices
    REQUIRE(actual.rows() == expected.rows());
    REQUIRE(actual.cols() == expected.cols());
    CHECK((actual - expected).norm() <= 1e-13 * expected.norm());
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
    "MeshFEM assembly exposes the block-CSC matrix",
    "[potential][assembly][meshfem]")
{
    // A downstream user should be able to consume the native block-CSC matrix
    // (e.g., to hand it to MeshFEM's solvers) instead of paying for the Eigen
    // conversion. This also checks that the forward declaration in our header
    // is enough: MeshFEMSparse's header is included by this test, not by
    // meshfem_hessian_assembler.hpp.
    const auto& spec = ipc::tests::assembly_scene_specs().at(0); // two-cubes

    const std::optional<ipc::tests::AssemblyScene> maybe_scene =
        ipc::tests::build_assembly_scene(spec);
    if (!maybe_scene.has_value()) {
        SKIP(fmt::format("Scene '{}' is unavailable.", spec.mesh_name));
    }
    const ipc::tests::AssemblyScene& scene = maybe_scene.value();

    MeshFEMHessianAssembler assembler;
    scene.potential().assemble_hessian(
        scene.collisions(), scene.mesh(), scene.vertices(), assembler);

    const MeshFEM::BlockCSCHessianBase& H = assembler.block_matrix();
    const Eigen::SparseMatrix<double> H_eigen = assembler.get_matrix();

    REQUIRE(H.numScalarCols() == size_t(H_eigen.cols()));
    // Only the upper triangle is stored, so the block matrix holds fewer
    // scalar entries than the symmetrized Eigen matrix.
    CHECK(H.scalarNNZ() <= size_t(H_eigen.nonZeros()));

    // trace() must skip the empty block columns of a contact pattern rather
    // than reading the preceding column's storage.
    CHECK_THAT(
        H.trace(),
        Catch::Matchers::WithinRel(Eigen::MatrixXd(H_eigen).trace(), 1e-12));

    // addDiag()/setDiag() cannot work without every diagonal block present,
    // so they must reject a contact pattern instead of corrupting it.
    const Eigen::VectorXd ones =
        Eigen::VectorXd::Ones(Eigen::Index(H.numScalarCols()));
    CHECK_THROWS(const_cast<MeshFEM::BlockCSCHessianBase&>(H).addDiag(ones));

    // Exercise MeshFEM's own API on the result: y = H x must match Eigen's.
    const Eigen::VectorXd x =
        Eigen::VectorXd::Random(Eigen::Index(H.numScalarCols()));
    const Eigen::VectorXd y_expected = H_eigen * x;
    const Eigen::VectorXd y = H.apply(x);
    CHECK((y - y_expected).norm() <= 1e-12 * std::max(1.0, y_expected.norm()));
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
