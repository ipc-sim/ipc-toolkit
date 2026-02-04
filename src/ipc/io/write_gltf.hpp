#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <string>
#include <vector>

namespace ipc::rigid {

/// @brief Write a sequence of rigid body poses to a glTF file.
/// @param filename The output glTF filename.
/// @param bodies The rigid bodies to write.
/// @param poses A list of poses for each timestep.
/// @param timestep The time interval between each pose in seconds.
/// @param embed_buffers Whether to embed the binary buffers in the glTF file.
/// @param write_binary Whether to write a binary .glb file (true) or a text .gltf file (false).
/// @param prettyPrint Whether to pretty-print the JSON content.
/// @return True if successful, false otherwise.
bool write_gltf(
    const std::string& filename,
    const RigidBodies& bodies,
    const std::list<std::vector<Pose>>& poses,
    double timestep,
    bool embed_buffers = true,
    bool write_binary = true,
    bool prettyPrint = true);

} // namespace ipc::rigid
