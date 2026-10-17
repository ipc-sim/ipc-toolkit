// NOTE: This is an internal header file, not meant to be used outside of the
// IPC Toolkit library. It includes ipc/utils/unordered_map_and_set.hpp, which
// pulls in Abseil, a private dependency of the IPC Toolkit library.

#pragma once

#include <ipc/esp/esp_collisions.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace ipc {

/// @brief Collisions of one quadrature point, keyed by
/// ESPCollision::get_typed_hash().
struct ESPCollisionPairMap
    : unordered_map<std::array<index_t, 3>, std::shared_ptr<ESPCollision>> { };

struct ESPCollisions::Maps {
    /// @brief collision sets for 3D quadrature

    // vertex_collisions[vi] provides the contact set for vertex vi
    unordered_map<index_t, std::unique_ptr<ESPCollisionDict<PointType::VERTEX>>>
        vertex_collisions;
    // edge_edge_collisions[(ei, ej)] provides the contact set for the closest
    // point on ei, between edge ei and ej.
    unordered_map<
        std::pair<index_t, index_t>,
        std::unique_ptr<ESPCollisionDict<PointType::EDGE>>>
        edge_edge_collisions;
    // face_collisions[fi][qi] provides the contact set for quadrature point qi
    // of face fi
    unordered_map<
        index_t,
        std::vector<std::unique_ptr<ESPCollisionDict<PointType::FACE>>>>
        face_collisions;

    /// @brief collision sets for 2D quadrature
    // edge_collisions_2d[ei][qi] provides the contact set for Gauss-Lobatto QP
    // qi on edge ei
    unordered_map<
        index_t,
        std::vector<std::unique_ptr<ESPCollisionDict<PointType::EDGE, 2>>>>
        edge_collisions_2d;
};

} // namespace ipc
