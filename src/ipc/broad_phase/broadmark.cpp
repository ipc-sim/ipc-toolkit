#include <ipc/broad_phase/broadmark.hpp>

namespace ipc {

template <class T>
void Broadmark<T>::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(edges, faces);
    num_vertices = vertices.rows();
    interface.ConstructBoxes(vertices, vertices, edges, faces, boxes);
    interface.CalcOverlaps(
        vertices, vertices, edges, faces, boxes, /*init=*/false);
}

template <class T>
void Broadmark<T>::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(edges, faces);
    num_vertices = vertices_t0.rows();
    interface.ConstructBoxes(vertices_t0, vertices_t1, edges, faces, boxes);
    interface.CalcOverlaps(
        vertices_t0, vertices_t1, edges, faces, boxes, /*init=*/false);
}

template <class T> void Broadmark<T>::clear()
{
    num_vertices = 0;
    interface.Clear();
    boxes.clear();
}

template <class T>
void Broadmark<T>::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

template <class T>
void Broadmark<T>::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    // interface.FilterOverlaps(num_vertices, edges, faces);
    // candidates = interface.candidates.ee_candidates;
}

template <class T>
void Broadmark<T>::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    // interface.FilterOverlaps(num_vertices, edges, faces);
    // candidates = interface.candidates.fv_candidates;
}

template <class T>
void Broadmark<T>::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    // interface.FilterOverlaps(num_vertices, edges, faces);
    // candidates = interface.candidates.ef_candidates;
}

template <class T>
void Broadmark<T>::detect_collision_candidates(
    int dim, Candidates& candidates) const
{
    interface.FilterOverlaps(num_vertices, edges, faces, candidates);
    // interface.m_broadPhase = candidates.size();
}

// Explicitly instantiate the template for SAP type
template class Broadmark<SAP>;
template class Broadmark<Grid_3D>;
template class Broadmark<GPU_LBVH>;

} // namespace ipc