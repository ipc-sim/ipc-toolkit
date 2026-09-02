#pragma once

#include "assembly_scene.hpp"

#include <cstddef>
#include <cstdint>

namespace ipc::tests {

/// @brief Which collision type dominates a synthetic contact scene.
enum class SyntheticContactType : std::uint8_t {
    FACE_VERTEX, ///< 3D only.
    EDGE_EDGE,   ///< 3D only.
    EDGE_VERTEX, ///< 2D only.
};

/// @brief Short label for a synthetic contact type.
inline const char* to_string(const SyntheticContactType type)
{
    switch (type) {
    case SyntheticContactType::FACE_VERTEX:
        return "fv";
    case SyntheticContactType::EDGE_EDGE:
        return "ee";
    default:
        return "ev";
    }
}

/// @brief Description of a synthetic contact scene.
///
/// Unlike the mesh-based scenes in assembly_scene.hpp, synthetic scenes
/// decouple the collision count from the DOF count so benchmarks can sweep
/// them independently. All vertices are placed inside a unit box and the
/// barrier activation distance exceeds the box diagonal, so every sampled
/// collision is active regardless of which vertices it touches.
struct SyntheticContactSpec {
    /// @brief Number of collisions to sample.
    size_t num_collisions;
    /// @brief Desired collision-mesh DOF (rounded to a multiple of dim).
    size_t target_ndof;
    /// @brief Spatial dimension (2 or 3; constrains the collision type).
    int dim = 3;
    /// @brief Which collision type to sample.
    SyntheticContactType type = SyntheticContactType::FACE_VERTEX;
    /// @brief Interior padding ratio (see AssemblySceneSpec).
    double interior_vertex_ratio = 1.0;
    /// @brief If true, construct the mesh with an explicit (identity)
    /// displacement map. The map is semantically a no-op, but its presence
    /// makes CollisionMesh::is_selection_dof_map() false, which forces the
    /// map_to_full path in Potential::gradient(in_full_dof=true).
    bool with_displacement_map = false;
    /// @brief RNG seed for vertex positions and collision sampling.
    unsigned seed = 0;
};

/// @brief Build a synthetic contact scene (see SyntheticContactSpec).
AssemblyScene build_synthetic_contact_scene(const SyntheticContactSpec& spec);

} // namespace ipc::tests
