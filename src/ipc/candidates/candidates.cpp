#include "candidates.hpp"

#include <ipc/utils/save_obj.hpp>

#include <fstream>

namespace ipc {

size_t Candidates::size() const
{
    return ev_candidates.size() + ee_candidates.size() + fv_candidates.size();
}

bool Candidates::empty() const
{
    return ev_candidates.empty() && ee_candidates.empty()
        && fv_candidates.empty();
}

void Candidates::clear()
{
    ev_candidates.clear();
    ee_candidates.clear();
    fv_candidates.clear();
}

ContinuousCollisionCandidate& Candidates::operator[](size_t idx)
{
    if (idx < ev_candidates.size()) {
        return ev_candidates[idx];
    }
    idx -= ev_candidates.size();
    if (idx < ee_candidates.size()) {
        return ee_candidates[idx];
    }
    idx -= ee_candidates.size();
    if (idx < fv_candidates.size()) {
        return fv_candidates[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

const ContinuousCollisionCandidate& Candidates::operator[](size_t idx) const
{
    if (idx < ev_candidates.size()) {
        return ev_candidates[idx];
    }
    idx -= ev_candidates.size();
    if (idx < ee_candidates.size()) {
        return ee_candidates[idx];
    }
    idx -= ee_candidates.size();
    if (idx < fv_candidates.size()) {
        return fv_candidates[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

bool Candidates::save_obj(
    const std::string& filename,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::ofstream obj(filename, std::ios::out);
    if (!obj.is_open()) {
        return false;
    }
    int v_offset = 0;
    ipc::save_obj(obj, V, E, F, ev_candidates, v_offset);
    v_offset += ev_candidates.size() * 3;
    ipc::save_obj(obj, V, E, F, ee_candidates, v_offset);
    v_offset += ee_candidates.size() * 4;
    ipc::save_obj(obj, V, F, F, fv_candidates, v_offset);
    return true;
}

} // namespace ipc
