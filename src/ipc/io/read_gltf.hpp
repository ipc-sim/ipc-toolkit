#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ipc::rigid {

// Forward declarations
class RigidBodies;
struct Pose;

/// @brief Read a rigid body scene from a glTF file and return a RigidBodies object.
///
/// We only read geometry: each mesh node becomes a body named after its node
/// (or mesh), made unique with RigidBodies::make_unique_names over the mesh
/// nodes in file order, posed by its translation and rotation, with its scale
/// baked into the rest positions. Nodes scaled to zero are skipped. Every body
/// gets a default density of 1000; set physical properties on the returned
/// bodies.
///
/// The one exception is @p is_dof_fixed. Fixing exactly two rotational DOFs
/// changes the frame a body's inertia is expressed in, which is baked into the
/// body at construction, so we take the masks here instead of through a setter.
///
/// @param filename The input glTF filename.
/// @param convert_planes If true, any mesh with exactly two coplanar faces will be converted to a plane.
/// @param is_dof_fixed Per-DOF fixed flags ([position | rotation]) keyed by body name. Bodies not listed have all DOFs free.
/// @return A pair containing the RigidBodies object and a vector of initial poses for each body.
std::pair<std::shared_ptr<RigidBodies>, std::vector<Pose>> read_gltf(
    const std::string& filename,
    const bool convert_planes = false,
    const std::unordered_map<std::string, VectorMax6b>& is_dof_fixed = {});

} // namespace ipc::rigid
