#include <ipc/broad_phase/collision_candidate.hpp>

#include <algorithm> // std::min/max

#include <ipc/utils/logger.hpp>

namespace ipc {

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

///////////////////////////////////////////////////////////////////////////////

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

bool Candidates::save_obj(
    const std::string& filename,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F)
{
    std::ofstream obj(filename, std::ios::out);
    if (!obj.is_open()) {
        return false;
    }
    ipc::save_obj(obj, V, E, F, ev_candidates);
    ipc::save_obj(obj, V, E, F, ee_candidates);
    ipc::save_obj(obj, V, F, F, fv_candidates);
    return true;
}

///////////////////////////////////////////////////////////////////////////////

static const Eigen::IOFormat OBJ_VERTEX_FORMAT(
    Eigen::FullPrecision, Eigen::DontAlignCols, " ", "", "v ", "\n", "", "");

void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeVertexCandidate>& ev_candidates)
{
    out << "o EV\n";
    int i = 1;
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
    const std::vector<EdgeEdgeCandidate>& ee_candidates)
{
    out << "o EE\n";
    int i = 1;
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
    const std::vector<FaceVertexCandidate>& fv_candidates)
{
    out << "o FV\n";
    int i = 1;
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
    const std::vector<EdgeFaceCandidate>& ef_candidates)
{
    out << "o EF\n";
    int i = 1;
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
