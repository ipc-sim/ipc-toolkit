#include "save_obj.hpp"

#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/edge_face.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

template <>
void save_obj(
    std::ostream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi&,
    const Eigen::MatrixXi&,
    const std::vector<VertexVertexCandidate>& vv_candidates,
    const int)
{
    out << "o VV\n";
    for (const auto& vv_candidate : vv_candidates) {
        out << V.row(vv_candidate.vertex0_id).format(OBJ_VERTEX_FORMAT);
        out << V.row(vv_candidate.vertex1_id).format(OBJ_VERTEX_FORMAT);
    }
}

template <>
void save_obj(
    std::ostream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeVertexCandidate>& ev_candidates,
    const int v_offset)
{
    out << "o EV\n";
    int i = v_offset + 1;
    for (const auto& ev_candidate : ev_candidates) {
        out << V.row(E(ev_candidate.edge_id, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ev_candidate.edge_id, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(ev_candidate.vertex_id).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("l {:d} {:d}\n", i, i + 1);
        i += 3;
    }
}

template <>
void save_obj(
    std::ostream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeEdgeCandidate>& ee_candidates,
    const int v_offset)
{
    out << "o EE\n";
    int i = v_offset + 1;
    for (const auto& ee_candidate : ee_candidates) {
        out << V.row(E(ee_candidate.edge0_id, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ee_candidate.edge0_id, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ee_candidate.edge1_id, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ee_candidate.edge1_id, 1)).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("l {:d} {:d}\n", i + 0, i + 1);
        out << fmt::format("l {:d} {:d}\n", i + 2, i + 3);
        i += 4;
    }
}

template <>
void save_obj(
    std::ostream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<FaceVertexCandidate>& fv_candidates,
    const int v_offset)
{
    out << "o FV\n";
    int i = v_offset + 1;
    for (const auto& fv_candidate : fv_candidates) {
        out << V.row(F(fv_candidate.face_id, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(fv_candidate.face_id, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(fv_candidate.face_id, 2)).format(OBJ_VERTEX_FORMAT);
        out << V.row(fv_candidate.vertex_id).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("f {:d} {:d} {:d}\n", i, i + 1, i + 2);
        i += 4;
    }
}

template <>
void save_obj(
    std::ostream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeFaceCandidate>& ef_candidates,
    const int v_offset)
{
    out << "o EF\n";
    int i = v_offset + 1;
    for (const auto& ef_candidate : ef_candidates) {
        out << V.row(E(ef_candidate.edge_id, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ef_candidate.edge_id, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(ef_candidate.face_id, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(ef_candidate.face_id, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(ef_candidate.face_id, 2)).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("l {:d} {:d}\n", i, i + 1);
        out << fmt::format("f {:d} {:d} {:d}\n", i + 2, i + 3, i + 4);
        i += 5;
    }
}

} // namespace ipc