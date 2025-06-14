#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <fstream>
#include <ostream>
#include <vector>

namespace ipc {

template <typename Candidate>
void save_obj(
    std::ostream& out,
    Eigen::ConstRef<Eigen::MatrixXd> V,
    Eigen::ConstRef<Eigen::MatrixXi> E,
    Eigen::ConstRef<Eigen::MatrixXi> F,
    const std::vector<Candidate>& candidates,
    const int v_offset = 0);

template <typename Candidate>
bool save_obj(
    const std::string& filename,
    Eigen::ConstRef<Eigen::MatrixXd> V,
    Eigen::ConstRef<Eigen::MatrixXi> E,
    Eigen::ConstRef<Eigen::MatrixXi> F,
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
