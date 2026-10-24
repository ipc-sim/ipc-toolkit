#pragma once

#include <ipc/collision_mesh.hpp>

#include <vector>

namespace ipc {

/// @brief Compute the inclusion-exclusion weight of every edge and vertex of
/// a mesh (every face has weight one).
///
/// We choose the weights so that the weights of the elements containing any
/// point of the mesh sum to one (supplemental S2). Faces take one, so an edge
/// takes 1 - (faces on it) and a vertex takes 1 - (edges at it) + (faces at
/// it), since each face at a vertex has two of its edges there. A closed mesh
/// gets the alternating signs (3D: faces +1, edges -1, vertices +1; 2D: edges
/// +1, vertices -1). A boundary edge or vertex, and an open polyline's end,
/// gets zero (S4). An edge in no face and an isolated vertex get one.
///
/// Using the closed-mesh signs everywhere instead is wrong on any other
/// input: beyond a straight boundary edge of an open sheet the face and edge
/// terms cancel, so the potential is exactly zero where it should be b(d), and
/// a free edge or isolated vertex gets a negative weight, so its potential is
/// attractive and unbounded below.
///
/// @param mesh The collision mesh.
/// @param[out] edge_weights Weight of each edge (size mesh.num_edges()).
/// @param[out] vertex_weights Weight of each vertex (size mesh.num_vertices()).
void compute_inclusion_exclusion_weights(
    const CollisionMesh& mesh,
    std::vector<int>& edge_weights,
    std::vector<int>& vertex_weights);

} // namespace ipc
