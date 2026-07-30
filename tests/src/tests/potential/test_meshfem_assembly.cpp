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
