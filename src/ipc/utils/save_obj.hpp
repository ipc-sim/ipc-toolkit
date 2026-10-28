#pragma once

// Deprecated: this header moved to <ipc/io/write_candidates_obj.hpp>, and
// save_obj was renamed to write_candidates_obj. We keep these forwarding
// overloads so existing code builds (with a deprecation warning).

#include <ipc/io/write_candidates_obj.hpp>

namespace ipc {

template <typename Candidate>
[[deprecated(
    "Use ipc::write_candidates_obj from <ipc/io/write_candidates_obj.hpp> instead.")]]
void save_obj(
    std::ostream& out,
    Eigen::ConstRef<Eigen::MatrixXd> V,
    Eigen::ConstRef<Eigen::MatrixXi> E,
    Eigen::ConstRef<Eigen::MatrixXi> F,
    const std::vector<Candidate>& candidates,
    const int v_offset = 0)
{
    write_candidates_obj(out, V, E, F, candidates, v_offset);
}

template <typename Candidate>
[[deprecated(
    "Use ipc::write_candidates_obj from <ipc/io/write_candidates_obj.hpp> instead.")]]
bool save_obj(
    const std::string& filename,
    Eigen::ConstRef<Eigen::MatrixXd> V,
    Eigen::ConstRef<Eigen::MatrixXi> E,
    Eigen::ConstRef<Eigen::MatrixXi> F,
    const std::vector<Candidate>& candidates)
{
    return write_candidates_obj(filename, V, E, F, candidates);
}

} // namespace ipc
