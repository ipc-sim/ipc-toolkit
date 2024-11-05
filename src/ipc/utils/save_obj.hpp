#pragma once

#include <Eigen/Core>

#include <fstream>
#include <ostream>
#include <vector>

namespace ipc {

template <typename Candidate>
void save_obj(
    std::ostream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<Candidate>& candidates,
    const int v_offset = 0);

template <typename Candidate>
bool save_obj(
    const std::string& filename,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<Candidate>& candidates)
{
    std::ofstream obj(filename);
    if (!obj.is_open()) {
        return false;
    }
    save_obj(obj, V, E, F, candidates);
    return true;
}

} // namespace ipc
