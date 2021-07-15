#include <ipc/utils/faces_to_edges.hpp>

#include <ipc/utils/unordered_map_and_set.hpp>

namespace ipc {

Eigen::MatrixXi
faces_to_edges(const Eigen::MatrixXi& F, const Eigen::MatrixXi& E)
{
    auto max_vi = E.maxCoeff();
    auto hash_edge = [&max_vi](const Eigen::RowVector2i& e) {
        return e.minCoeff() * max_vi + e.maxCoeff();
    };
    auto eq_edges = [](const Eigen::RowVector2i& ea,
                       const Eigen::RowVector2i& eb) {
        return ea.minCoeff() == eb.minCoeff() && ea.maxCoeff() == eb.maxCoeff();
    };

    unordered_map<
        Eigen::RowVector2i, int, decltype(hash_edge), decltype(eq_edges)>
        edge_map(/*bucket_count=*/E.rows(), hash_edge, eq_edges);
    for (int ei = 0; ei < E.rows(); ei++) {
        edge_map[E.row(ei)] = ei;
    }

    Eigen::MatrixXi F2E(F.rows(), F.cols());
    for (int fi = 0; fi < F.rows(); fi++) {
        for (int fj = 0; fj < F.cols(); fj++) {
            Eigen::RowVector2i e(F(fi, fj), F(fi, (fj + 1) % F.cols()));
            auto search = edge_map.find(e);
            if (search != edge_map.end()) {
                F2E(fi, fj) = search->second;
            } else {
                throw std::runtime_error("Unable to find edge!");
            }
        }
    }

    return F2E;
}

} // namespace ipc
