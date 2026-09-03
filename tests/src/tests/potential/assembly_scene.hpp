#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/utils/gradient_assembler.hpp>

#include <optional>
#include <string>
#include <vector>

namespace ipc::tests {

/// @brief Which strategy `assemble_gradient` uses for a given shape.
///
/// Delegates to `ipc::gradient_assembly_is_two_pass`, supplying the
/// normal-collision stencil size.
///
/// @param num_collisions Number of collisions being assembled.
/// @param out_ndof Size of the assembled gradient.
inline bool
assembly_is_two_pass(const size_t num_collisions, const size_t out_ndof)
{
    return gradient_assembly_is_two_pass(
        NormalCollision::STENCIL_SIZE * num_collisions, out_ndof);
}

/// @brief Description of a contact scene used to benchmark assembly.
struct AssemblySceneSpec {
    /// @brief Name of the mesh file relative to the test data directory.
    std::string mesh_name;
    /// @brief Barrier activation distance.
    double dhat;
    /// @brief Number of synthetic interior vertices per collision vertex.
    ///
    /// The collision meshes in the test data are surfaces, so every vertex ends
    /// up in the collision mesh and the selection matrix used by
    /// `CollisionMesh::to_full_dof` is (a permutation of) the identity. Real
    /// users embed the surface in a volumetric FE mesh whose interior nodes are
    /// absent from the collision mesh. We emulate that by padding the full mesh
    /// with vertices that no edge or face references, which
    /// `construct_is_on_surface` then excludes. A ratio of 1.0 means the full
    /// mesh has twice as many vertices as the collision mesh.
    double interior_vertex_ratio = 1.0;
    /// @brief Short label used in benchmark names.
    std::string label;
};

/// @brief A pre-built contact scene: mesh, positions, collisions, potential.
///
/// Building the collision set is deliberately *not* part of what the assembly
/// benchmarks measure, so it happens once here.
class AssemblyScene {
public:
    /// @note `CollisionMesh` and `BarrierPotential` declare destructors, which
    /// suppresses their implicit move constructors, so they are taken by
    /// const reference rather than by value.
    AssemblyScene(
        std::string label,
        const CollisionMesh& mesh,
        Eigen::MatrixXd vertices,
        NormalCollisions collisions,
        const BarrierPotential& potential);

    const std::string& label() const { return m_label; }
    const CollisionMesh& mesh() const { return m_mesh; }
    const Eigen::MatrixXd& vertices() const { return m_vertices; }
    const NormalCollisions& collisions() const { return m_collisions; }
    const BarrierPotential& potential() const { return m_potential; }

    size_t num_collisions() const { return m_collisions.size(); }
    size_t num_vertices() const { return m_mesh.num_vertices(); }
    size_t full_num_vertices() const { return m_mesh.full_num_vertices(); }
    size_t ndof() const { return m_mesh.ndof(); }
    size_t full_ndof() const { return m_mesh.full_ndof(); }

    /// @brief A human-readable one-line summary of the scene's size.
    std::string stats() const;

private:
    std::string m_label;
    CollisionMesh m_mesh;
    Eigen::MatrixXd m_vertices;
    NormalCollisions m_collisions;
    BarrierPotential m_potential;
};

/// @brief Build a scene, or return nullopt if its mesh is unavailable.
///
/// Some test meshes are private, so callers must handle absence by skipping.
std::optional<AssemblyScene>
build_assembly_scene(const AssemblySceneSpec& spec);

/// @brief The standard set of scenes used by the assembly benchmarks.
///
/// Ordered by increasing collision count so a truncated run is still useful.
const std::vector<AssemblySceneSpec>& assembly_scene_specs();

} // namespace ipc::tests
