#include "broadmark.hpp"
#ifdef IPC_TOOLKIT_WITH_BROADMARK

#include <Broadphase/Algorithms/GPU/Bullet3GPUAlgorithms.h>
#include <Broadphase/Algorithms/BF/BF.h>
#include <Broadphase/Algorithms/DBVT/DBVT.h>
#include <Broadphase/Algorithms/Grid/Grid_ND.h>
#include <Broadphase/Algorithms/Grid/Grid_ND_Parallel.h>
#include <Broadphase/Algorithms/Grid/Grid_ND_SAP.h>
#include <Broadphase/Algorithms/KD/KD.h>
#include <Broadphase/Algorithms/SAP/SAP.h>
#include <Broadphase/Algorithms/SAP/SAP_Parallel.h>
#include <Broadphase/Algorithms/SAP/SAP_SIMD_Parallel.h>
#include <Broadphase/Algorithms/Tracy/Tracy.h>
#include <Broadphase/Algorithms/Tracy/Tracy_Parallel.h>
#include <Broadphase/Algorithms/iSAP/AxisSweep.h>
#include <Broadphase/Algorithms/CGAL/CGAL.h>

namespace ipc {

void CopyMeshBroadPhase::copy_mesh(
    const Eigen::MatrixXi& p_edges, const Eigen::MatrixXi& p_faces)
{
    edges = p_edges;
    faces = p_faces;
}

bool CopyMeshBroadPhase::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const long e0i = edges(ei, 0), e1i = edges(ei, 1);

    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool CopyMeshBroadPhase::can_edges_collide(size_t eai, size_t ebi) const
{
    const long ea0i = edges(eai, 0), ea1i = edges(eai, 1);
    const long eb0i = edges(ebi, 0), eb1i = edges(ebi, 1);

    const bool share_endpoint =
        ea0i == eb0i || ea0i == eb1i || ea1i == eb0i || ea1i == eb1i;

    return !share_endpoint
        && (can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
            || can_vertices_collide(ea1i, eb0i)
            || can_vertices_collide(ea1i, eb1i));
}

bool CopyMeshBroadPhase::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const long f0i = faces(fi, 0), f1i = faces(fi, 1), f2i = faces(fi, 2);

    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool CopyMeshBroadPhase::can_edge_face_collide(size_t ei, size_t fi) const
{
    const long e0i = edges(ei, 0), e1i = edges(ei, 1);
    const long f0i = faces(fi, 0), f1i = faces(fi, 1), f2i = faces(fi, 2);

    const bool share_endpoint = e0i == f0i || e0i == f1i || e0i == f2i
        || e1i == f0i || e1i == f1i || e1i == f2i;

    return !share_endpoint
        && (can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
            || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
            || can_vertices_collide(e1i, f1i)
            || can_vertices_collide(e1i, f2i));
}

// ============================================================================

template <class T>
void Broadmark<T>::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& p_edges,
    const Eigen::MatrixXi& p_faces,
    double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(p_edges, p_faces);
    BroadPhase::build(vertices, edges, faces, inflation_radius);
    num_vertices = vertices.rows();
    std::vector<ipc::AABB> ipc_aabbs;
    ipc::combine_aabbs(edge_boxes, face_boxes, vertex_boxes, ipc_aabbs);
    ipc::to_broadmark_aabbs(ipc_aabbs, boxes);
    ipc_aabbs.clear();
    interface.CalcOverlaps(vertices, edges, faces, boxes, /*init=*/true);
}

template <class T>
void Broadmark<T>::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& p_edges,
    const Eigen::MatrixXi& p_faces,
    double inflation_radius)
{
    CopyMeshBroadPhase::copy_mesh(p_edges, p_faces);
    BroadPhase::build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
    num_vertices = vertices_t0.rows();
    std::vector<ipc::AABB> ipc_aabbs;
    ipc::combine_aabbs(edge_boxes, face_boxes, vertex_boxes, ipc_aabbs);
    ipc::to_broadmark_aabbs(ipc_aabbs, boxes);
    ipc_aabbs.clear();
    interface.CalcOverlaps(vertices_t0, edges, faces, boxes, /*init=*/true);
}

template <class T> void Broadmark<T>::clear()
{
    BroadPhase::clear();
    num_vertices = 0;
    interface.Clear();
    boxes.clear();
}

template <class T>
void Broadmark<T>::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
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
    throw std::runtime_error("Not implemented!");
}

template <class T>
void Broadmark<T>::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

template <class T>
void Broadmark<T>::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

template <class T>
void Broadmark<T>::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    throw std::runtime_error("Not implemented!");
}

template <class T>
void Broadmark<T>::detect_collision_candidates(
    int dim, Candidates& candidates) const
{
    if (dim == 2)
        throw std::runtime_error("Not implemented!");
    interface.FilterOverlaps(num_vertices, edges, faces, candidates);
    // interface.m_broadPhase = candidates.size();
}

// Explicitly instantiate the template for SAP type
template class Broadmark<SAP>;
template class Broadmark<SAP_Parallel>;
template class Broadmark<Grid_3D>;
template class Broadmark<Grid_3D_Parallel>;
template class Broadmark<DBVT_D>;
template class Broadmark<DBVT_F>;
template class Broadmark<AxisSweep>;
template class Broadmark<KD>;
template class Broadmark<Tracy>;
template class Broadmark<Tracy_Parallel>;
template class Broadmark<Grid_3D_SAP>;
template class Broadmark<CGAL_Internal>;
template class Broadmark<GPU_Grid>;
template class Broadmark<GPU_LBVH>;
template class Broadmark<GPU_SAP>;

} // namespace ipc
#endif