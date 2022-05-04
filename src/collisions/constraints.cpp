#include <ipc/collisions/constraints.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

size_t Constraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size() + pv_constraints.size();
}

bool Constraints::empty() const
{
    return vv_constraints.empty() && ev_constraints.empty()
        && ee_constraints.empty() && fv_constraints.empty()
        && pv_constraints.empty();
}

void Constraints::clear()
{
    vv_constraints.clear();
    ev_constraints.clear();
    ee_constraints.clear();
    fv_constraints.clear();
    pv_constraints.clear();
}

CollisionConstraint& Constraints::operator[](size_t idx)
{
    if (idx < vv_constraints.size()) {
        return vv_constraints[idx];
    }
    idx -= vv_constraints.size();
    if (idx < ev_constraints.size()) {
        return ev_constraints[idx];
    }
    idx -= ev_constraints.size();
    if (idx < ee_constraints.size()) {
        return ee_constraints[idx];
    }
    idx -= ee_constraints.size();
    if (idx < fv_constraints.size()) {
        return fv_constraints[idx];
    }
    idx -= fv_constraints.size();
    if (idx < pv_constraints.size()) {
        return pv_constraints[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

const CollisionConstraint& Constraints::operator[](size_t idx) const
{
    if (idx < vv_constraints.size()) {
        return vv_constraints[idx];
    }
    idx -= vv_constraints.size();
    if (idx < ev_constraints.size()) {
        return ev_constraints[idx];
    }
    idx -= ev_constraints.size();
    if (idx < ee_constraints.size()) {
        return ee_constraints[idx];
    }
    idx -= ee_constraints.size();
    if (idx < fv_constraints.size()) {
        return fv_constraints[idx];
    }
    idx -= fv_constraints.size();
    if (idx < pv_constraints.size()) {
        return pv_constraints[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

} // namespace ipc
