#pragma once

#include <memory>
#include <string>
#include <vector>

namespace ipc::rigid {

// Forward declarations
class RigidBodies;
struct Pose;

/// @brief Read a rigid body scene from a glTF file and return a RigidBodies object.
/// @param filename The input glTF filename.
/// @param convert_planes If true, any mesh with exactly two coplanar faces will be converted to a plane.
/// @return A pair containing the RigidBodies object and a vector of initial poses for each body.
std::pair<std::shared_ptr<RigidBodies>, std::vector<Pose>>
read_gltf(const std::string& filename, const bool convert_planes = false);

} // namespace ipc::rigid
