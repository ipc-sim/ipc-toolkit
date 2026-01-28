#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <string>
#include <vector>

namespace ipc::rigid {

bool write_gltf(
    const std::string& filename,
    const RigidBodies& bodies,
    const std::list<std::vector<Pose>>& poses,
    double timestep,
    bool embed_buffers = true,
    bool write_binary = true,
    bool prettyPrint = true);

} // namespace ipc::rigid
