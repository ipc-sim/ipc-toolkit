#pragma once

#include <vector>

#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/collisions/vertex_vertex.hpp>
#include <ipc/collisions/edge_vertex.hpp>
#include <ipc/collisions/edge_edge.hpp>
#include <ipc/collisions/face_vertex.hpp>
#include <ipc/collisions/plane_vertex.hpp>

namespace ipc {

struct Constraints {
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;
    std::vector<PlaneVertexConstraint> pv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    bool empty() const;

    void clear();

    CollisionConstraint& operator[](size_t idx);
    const CollisionConstraint& operator[](size_t idx) const;
};

} // namespace ipc
