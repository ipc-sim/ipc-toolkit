#include <ipc/broad_phase/collision_candidate.hpp>

#include <algorithm> // std::min/max
#include <iostream>

#include <ipc/utils/logger.hpp>

namespace ipc {

static const Eigen::IOFormat OBJ_VERTEX_FORMAT(
    Eigen::FullPrecision, Eigen::DontAlignCols, " ", "", "v ", "\n", "", "");

VertexVertexCandidate::VertexVertexCandidate(
    long vertex0_index, long vertex1_index)
    : vertex0_index(vertex0_index)
    , vertex1_index(vertex1_index)
{
}

bool VertexVertexCandidate::operator==(const VertexVertexCandidate& other) const
{
    return vertex0_index == other.vertex0_index
        && vertex1_index == other.vertex1_index;
}

bool VertexVertexCandidate::operator<(const VertexVertexCandidate& other) const
{
    if (vertex0_index == other.vertex0_index) {
        return vertex1_index < other.vertex1_index;
    }
    return vertex0_index < other.vertex0_index;
}

///////////////////////////////////////////////////////////////////////////////

EdgeVertexCandidate::EdgeVertexCandidate(long edge_index, long vertex_index)
    : edge_index(edge_index)
    , vertex_index(vertex_index)
{
}

bool EdgeVertexCandidate::operator==(const EdgeVertexCandidate& other) const
{
    return edge_index == other.edge_index && vertex_index == other.vertex_index;
}

bool EdgeVertexCandidate::operator<(const EdgeVertexCandidate& other) const
{
    if (edge_index == other.edge_index) {
        return vertex_index < other.vertex_index;
    }
    return edge_index < other.edge_index;
}

bool EdgeVertexCandidate::ccd(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double& toi,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return point_edge_ccd(
        // Point at t=0
        V0.row(vertex_index),
        // Edge at t=0
        V0.row(E(edge_index, 0)), V0.row(E(edge_index, 1)),
        // Point at t=1
        V1.row(vertex_index),
        // Edge at t=1
        V1.row(E(edge_index, 0)), V1.row(E(edge_index, 1)), //
        toi, tmax, tolerance, max_iterations, conservative_rescaling);
}

void EdgeVertexCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::cout << V0.row(E(edge_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

///////////////////////////////////////////////////////////////////////////////

EdgeEdgeCandidate::EdgeEdgeCandidate(long edge0_index, long edge1_index)
    : edge0_index(edge0_index)
    , edge1_index(edge1_index)
{
}

bool EdgeEdgeCandidate::operator==(const EdgeEdgeCandidate& other) const
{
    // (i, j) == (i, j) || (i, j) == (j, i)
    return (this->edge0_index == other.edge0_index
            && this->edge1_index == other.edge1_index)
        || (this->edge0_index == other.edge1_index
            && this->edge1_index == other.edge0_index);
}

bool EdgeEdgeCandidate::operator<(const EdgeEdgeCandidate& other) const
{
    long this_min = std::min(this->edge0_index, this->edge1_index);
    long other_min = std::min(other.edge0_index, other.edge1_index);
    if (this_min == other_min) {
        return std::max(this->edge0_index, this->edge1_index)
            < std::max(other.edge0_index, other.edge1_index);
    }
    return this_min < other_min;
}

bool EdgeEdgeCandidate::ccd(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double& toi,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return edge_edge_ccd(
        // Edge 1 at t=0
        V0.row(E(edge0_index, 0)), V0.row(E(edge0_index, 1)),
        // Edge 2 at t=0
        V0.row(E(edge1_index, 0)), V0.row(E(edge1_index, 1)),
        // Edge 1 at t=1
        V1.row(E(edge0_index, 0)), V1.row(E(edge0_index, 1)),
        // Edge 2 at t=1
        V1.row(E(edge1_index, 0)), V1.row(E(edge1_index, 1)), //
        toi, tmax, tolerance, max_iterations, conservative_rescaling);
}

void EdgeEdgeCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::cout << V0.row(E(edge0_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge0_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge1_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(E(edge1_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge0_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge0_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge1_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(E(edge1_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

///////////////////////////////////////////////////////////////////////////////

EdgeFaceCandidate::EdgeFaceCandidate(long edge_index, long face_index)
    : edge_index(edge_index)
    , face_index(face_index)
{
}

bool EdgeFaceCandidate::operator==(const EdgeFaceCandidate& other) const
{
    return edge_index == other.edge_index && face_index == other.face_index;
}

bool EdgeFaceCandidate::operator<(const EdgeFaceCandidate& other) const
{
    if (edge_index == other.edge_index) {
        return face_index < other.face_index;
    }
    return edge_index < other.edge_index;
}

///////////////////////////////////////////////////////////////////////////////

FaceVertexCandidate::FaceVertexCandidate(long face_index, long vertex_index)
    : face_index(face_index)
    , vertex_index(vertex_index)
{
}

bool FaceVertexCandidate::operator==(const FaceVertexCandidate& other) const
{
    return face_index == other.face_index && vertex_index == other.vertex_index;
}

bool FaceVertexCandidate::operator<(const FaceVertexCandidate& other) const
{
    if (face_index == other.face_index) {
        return vertex_index < other.vertex_index;
    }
    return face_index < other.face_index;
}

bool FaceVertexCandidate::ccd(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double& toi,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return point_triangle_ccd(
        // Point at t=0
        V0.row(vertex_index),
        // Triangle at t=0
        V0.row(F(face_index, 0)), V0.row(F(face_index, 1)),
        V0.row(F(face_index, 2)),
        // Point at t=1
        V1.row(vertex_index),
        // Triangle at t=1
        V1.row(F(face_index, 0)), V1.row(F(face_index, 1)),
        V1.row(F(face_index, 2)), //
        toi, tmax, tolerance, max_iterations, conservative_rescaling);
}

void FaceVertexCandidate::print_ccd_query(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    std::cout << V0.row(F(face_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(F(face_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(F(face_index, 2)).format(OBJ_VERTEX_FORMAT);
    std::cout << V0.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(F(face_index, 0)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(F(face_index, 1)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(F(face_index, 2)).format(OBJ_VERTEX_FORMAT);
    std::cout << V1.row(vertex_index).format(OBJ_VERTEX_FORMAT);
    std::cout << std::flush;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef IPC_TOOLKIT_WITH_CUDA
Candidates::Candidates(
    const std::vector<std::pair<int, int>>& overlaps,
    const std::vector<ccdgpu::Aabb> boxes)
{
    for (int i = 0; i < overlaps.size(); i++) {
        const std::pair<int, int>& overlap = overlaps[i];
        const ccdgpu::Aabb& boxA = boxes[overlap.first];
        const ccdgpu::Aabb& boxB = boxes[overlap.second];
        if (is_face(boxA) && is_vertex(boxB)) { // FV
            fv_candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        } else if (is_face(boxB) && is_vertex(boxA)) { // VF
            fv_candidates.emplace_back(boxB.ref_id, boxA.ref_id);
        } else if (is_edge(boxA) && is_edge(boxB)) { // EE
            ee_candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        } else if (is_edge(boxA) && is_vertex(boxB)) { // EV
            ev_candidates.emplace_back(boxA.ref_id, boxB.ref_id);
        } else if (is_edge(boxB) && is_vertex(boxA)) { // VE
            ev_candidates.emplace_back(boxB.ref_id, boxA.ref_id);
        } else { // VV, FF, FE, EF
            // Invalid collision candidate
        }
    }
}
#endif

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

///////////////////////////////////////////////////////////////////////////////

void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeVertexCandidate>& ev_candidates,
    const int v_offset)
{
    out << "o EV\n";
    int i = v_offset + 1;
    for (const auto& ev_candidate : ev_candidates) {
        out << V.row(E(ev_candidate.edge_index, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ev_candidate.edge_index, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(ev_candidate.vertex_index).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("l {:d} {:d}\n", i, i + 1);
        i += 3;
    }
}

void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeEdgeCandidate>& ee_candidates,
    const int v_offset)
{
    out << "o EE\n";
    int i = v_offset + 1;
    for (const auto& ee_candidate : ee_candidates) {
        out << V.row(E(ee_candidate.edge0_index, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ee_candidate.edge0_index, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ee_candidate.edge1_index, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ee_candidate.edge1_index, 1)).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("l {:d} {:d}\n", i + 0, i + 1);
        out << fmt::format("l {:d} {:d}\n", i + 2, i + 3);
        i += 4;
    }
}

void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<FaceVertexCandidate>& fv_candidates,
    const int v_offset)
{
    out << "o FV\n";
    int i = v_offset + 1;
    for (const auto& fv_candidate : fv_candidates) {
        out << V.row(F(fv_candidate.face_index, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(fv_candidate.face_index, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(fv_candidate.face_index, 2)).format(OBJ_VERTEX_FORMAT);
        out << V.row(fv_candidate.vertex_index).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("f {:d} {:d} {:d}\n", i, i + 1, i + 2);
        i += 4;
    }
}

void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeFaceCandidate>& ef_candidates,
    const int v_offset)
{
    out << "o EF\n";
    int i = v_offset + 1;
    for (const auto& ef_candidate : ef_candidates) {
        out << V.row(E(ef_candidate.edge_index, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(E(ef_candidate.edge_index, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(ef_candidate.face_index, 0)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(ef_candidate.face_index, 1)).format(OBJ_VERTEX_FORMAT);
        out << V.row(F(ef_candidate.face_index, 2)).format(OBJ_VERTEX_FORMAT);
        out << fmt::format("l {:d} {:d}\n", i, i + 1);
        out << fmt::format("f {:d} {:d} {:d}\n", i + 2, i + 3, i + 4);
        i += 5;
    }
}

} // namespace ipc
