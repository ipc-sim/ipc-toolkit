#include <ipc/collision_constraint.hpp>

namespace ipc {

VertexVertexConstraint::VertexVertexConstraint(
    long vertex0_index, long vertex1_index)
    : vertex0_index(vertex0_index)
    , vertex1_index(vertex1_index)
{
}

bool VertexVertexConstraint::
operator==(const VertexVertexConstraint& other) const
{
    return vertex0_index == other.vertex0_index
        && vertex1_index == other.vertex1_index;
}

/// @brief Compare EdgeVertexConstraints for sorting.
bool VertexVertexConstraint::
operator<(const VertexVertexConstraint& other) const
{
    if (vertex0_index == other.vertex0_index) {
        return vertex1_index < other.vertex1_index;
    }
    return vertex0_index < other.vertex0_index;
}

EdgeVertexConstraint::EdgeVertexConstraint(long edge_index, long vertex_index)
    : EdgeVertexCandidate(edge_index, vertex_index)
{
}

EdgeVertexConstraint::EdgeVertexConstraint(const EdgeVertexCandidate& candidate)
    : EdgeVertexCandidate(candidate)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    long edge0_index, long edge1_index, double eps_x)
    : EdgeEdgeCandidate(edge0_index, edge1_index)
    , eps_x(eps_x)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const EdgeEdgeCandidate& candidate, double eps_x)
    : EdgeEdgeCandidate(candidate)
    , eps_x(eps_x)
{
}

FaceVertexConstraint::FaceVertexConstraint(long face_index, long vertex_index)
    : FaceVertexCandidate(face_index, vertex_index)
{
}
FaceVertexConstraint::FaceVertexConstraint(const FaceVertexCandidate& candidate)
    : FaceVertexCandidate(candidate)
{
}

size_t Constraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size();
}

size_t Constraints::num_constraints() const
{
    size_t num_constraints = 0;
    for (const auto& vv_constraint : vv_constraints) {
        num_constraints += vv_constraint.multiplicity;
    }
    for (const auto& ev_constraint : ev_constraints) {
        num_constraints += ev_constraint.multiplicity;
    }
    for (const auto& ee_constraint : ee_constraints) {
        num_constraints += ee_constraint.multiplicity;
    }
    for (const auto& fv_constraint : fv_constraints) {
        num_constraints += fv_constraint.multiplicity;
    }
    return num_constraints;
}

void Constraints::clear()
{
    vv_constraints.clear();
    ev_constraints.clear();
    ee_constraints.clear();
    fv_constraints.clear();
}

} // namespace ipc
