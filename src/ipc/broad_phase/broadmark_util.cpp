#include <ipc/broad_phase/broadmark_util.hpp>
#include "Dependencies/Bullet3/Bullet3Collision/BroadPhaseCollision/b3DynamicBvhBroadphase.h"
#include "Dependencies/Bullet2/btAxisSweep3.h"
#include "Dependencies/Bullet3/Bullet3Common/b3AlignedAllocator.h"

namespace ipc {

void InterfaceBase::Clear()
{
    m_broadPhase = 0;
    m_Objects = 0;
    m_narrowPhase = 0;
    m_truePositive = 0;
    candidates.clear();
}

template <>
void Interface<DBVT_F>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    Candidates& candidates) const
{
    auto b3BroadphasePairs = algo->m_broadphase->getOverlappingPairCache()
                                 ->getOverlappingPairArray();
    std::cout << "b3BroadphasePairs: " << b3BroadphasePairs.size() << std::endl;

    auto is_vertex = [&](int ai) { return ai < num_vertices; };
    auto is_edge = [&](int ai) {
        return (ai >= num_vertices) && (ai < (num_vertices + edges.rows()));
    };
    auto is_face = [&](int ai) { return ai >= (num_vertices + edges.rows()); };
    auto is_endpoint = [&](int vi, int fi) {
        return vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi, 2);
    };
    auto has_common_endpoint = [&](int ei, int ej) {
        return edges(ei, 0) == edges(ej, 0) || edges(ei, 0) == edges(ej, 1)
            || edges(ei, 1) == edges(ej, 0) || edges(ei, 1) == edges(ej, 1);
    };

    auto addCandidate = [&](int ai, int bi) {
        // std::cout << "ai: " << ai << "bi: " << bi << std::endl;
        auto mi = ai < bi ? ai : bi;
        ai = std::min(ai, mi);
        bi = std::max(bi, mi);
        assert(ai < bi);

        if (is_vertex(ai) && is_face(bi)
            && !is_endpoint(ai, bi - num_vertices - edges.rows())) {
            candidates.fv_candidates.emplace_back(
                bi - num_vertices - edges.rows(), ai);
        } else if (
            is_edge(ai) && is_edge(bi)
            && !has_common_endpoint(ai - num_vertices, bi - num_vertices)) {
            candidates.ee_candidates.emplace_back(
                ai - num_vertices, bi - num_vertices);
        }
    };

    for (size_t i = 0; i < b3BroadphasePairs.size(); i++) {

        // proxy increments aabb m_id+1
        // see DBVT::Initialize in DBVT.cpp
        int xi = b3BroadphasePairs[i].x - 1;
        int yi = b3BroadphasePairs[i].y - 1;
        addCandidate(xi, yi);
    }

    //  remove duplicates

    std::sort(candidates.ee_candidates.begin(), candidates.ee_candidates.end());
    std::sort(candidates.fv_candidates.begin(), candidates.fv_candidates.end());

    candidates.ee_candidates.erase(
        unique(
            candidates.ee_candidates.begin(), candidates.ee_candidates.end()),
        candidates.ee_candidates.end());
    candidates.fv_candidates.erase(
        unique(
            candidates.fv_candidates.begin(), candidates.fv_candidates.end()),
        candidates.fv_candidates.end());

    // m_broadPhase = candidates.size();
    // std::cout << "m_broadPhase: " << m_broadPhase << std::endl;
}

template <>
void Interface<DBVT_D>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    Candidates& candidates) const
{
    // auto b3BroadphasePairs = algo->m_broadphase->getOverlappingPairCache()
    //                              ->getOverlappingPairArray();
    // std::cout << "b3BroadphasePairs: " << b3BroadphasePairs.size() <<
    // std::endl;

    // auto is_vertex = [&](int ai) { return ai < num_vertices; };
    // auto is_edge = [&](int ai) {
    //     return (ai >= num_vertices) && (ai < (num_vertices + edges.rows()));
    // };
    // auto is_face = [&](int ai) { return ai >= (num_vertices + edges.rows());
    // }; auto is_endpoint = [&](int vi, int fi) {
    //     return vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi,
    //     2);
    // };
    // auto has_common_endpoint = [&](int ei, int ej) {
    //     return edges(ei, 0) == edges(ej, 0) || edges(ei, 0) == edges(ej, 1)
    //         || edges(ei, 1) == edges(ej, 0) || edges(ei, 1) == edges(ej, 1);
    // };

    // auto addCandidate = [&](int ai, int bi) {
    //     // std::cout << "ai: " << ai << "bi: " << bi << std::endl;
    //     auto mi = ai < bi ? ai : bi;
    //     ai = std::min(ai, mi);
    //     bi = std::max(bi, mi);
    //     assert(ai < bi);

    //     if (is_vertex(ai) && is_face(bi)
    //         && !is_endpoint(ai, bi - num_vertices - edges.rows())) {
    //         candidates.fv_candidates.emplace_back(
    //             bi - num_vertices - edges.rows(), ai);
    //     } else if (
    //         is_edge(ai) && is_edge(bi)
    //         && !has_common_endpoint(ai - num_vertices, bi - num_vertices)) {
    //         candidates.ee_candidates.emplace_back(
    //             ai - num_vertices, bi - num_vertices);
    //     }
    // };

    // for (size_t i = 0; i < b3BroadphasePairs.size(); i++) {
    //     // proxy increments aabb m_id+1
    //     // see DBVT::Initialize in DBVT.cpp
    //     int xi = b3BroadphasePairs[i].x - 1;
    //     int yi = b3BroadphasePairs[i].y - 1;
    //     addCandidate(xi, yi);
    // }

    // //  remove duplicates

    // std::sort(candidates.ee_candidates.begin(),
    // candidates.ee_candidates.end());
    // std::sort(candidates.fv_candidates.begin(),
    // candidates.fv_candidates.end());

    // candidates.ee_candidates.erase(
    //     unique(
    //         candidates.ee_candidates.begin(),
    //         candidates.ee_candidates.end()),
    //     candidates.ee_candidates.end());
    // candidates.fv_candidates.erase(
    //     unique(
    //         candidates.fv_candidates.begin(),
    //         candidates.fv_candidates.end()),
    //     candidates.fv_candidates.end());

    // m_broadPhase = candidates.size();
    // std::cout << "m_broadPhase: " << m_broadPhase << std::endl;
}

template <>
void Interface<AxisSweep>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    Candidates& candidates) const
{
    // auto btBroadphasePairs = algo->m_broadphase->getOverlappingPairCache()
    //                              ->getOverlappingPairArray();
    // std::cout << "btBroadphasePairs: " << btBroadphasePairs.size() <<
    // std::endl;

    // auto is_vertex = [&](int ai) { return ai < num_vertices; };
    // auto is_edge = [&](int ai) {
    //     return (ai >= num_vertices) && (ai < (num_vertices + edges.rows()));
    // };
    // auto is_face = [&](int ai) { return ai >= (num_vertices + edges.rows());
    // }; auto is_endpoint = [&](int vi, int fi) {
    //     return vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi,
    //     2);
    // };
    // auto has_common_endpoint = [&](int ei, int ej) {
    //     return edges(ei, 0) == edges(ej, 0) || edges(ei, 0) == edges(ej, 1)
    //         || edges(ei, 1) == edges(ej, 0) || edges(ei, 1) == edges(ej, 1);
    // };

    // auto addCandidate = [&](int ai, int bi) {
    //     auto mi = ai < bi ? ai : bi;
    //     ai = std::min(ai, mi);
    //     bi = std::max(bi, mi);

    //     if (is_vertex(ai) && is_face(bi)
    //         && !is_endpoint(ai, bi - num_vertices - edges.rows())) {
    //         candidates.fv_candidates.emplace_back(
    //             bi - num_vertices - edges.rows(), ai);
    //     } else if (
    //         is_edge(ai) && is_edge(bi)
    //         && !has_common_endpoint(ai - num_vertices, bi - num_vertices)) {
    //         candidates.ee_candidates.emplace_back(
    //             ai - num_vertices, bi - num_vertices);
    //     }
    // };

    // for (size_t i = 0; i < btBroadphasePairs.size(); i++) {
    //     int xi = (btBroadphasePairs[i].m_pProxy0->m_uniqueId) - 1;
    //     int yi = (btBroadphasePairs[i].m_pProxy0->m_uniqueId) - 1;
    //     addCandidate(xi, yi);
    // }

    // //  remove duplicates

    // std::sort(candidates.ee_candidates.begin(),
    // candidates.ee_candidates.end());
    // std::sort(candidates.fv_candidates.begin(),
    // candidates.fv_candidates.end());

    // candidates.ee_candidates.erase(
    //     unique(
    //         candidates.ee_candidates.begin(),
    //         candidates.ee_candidates.end()),
    //     candidates.ee_candidates.end());
    // candidates.fv_candidates.erase(
    //     unique(
    //         candidates.fv_candidates.begin(),
    //         candidates.fv_candidates.end()),
    //     candidates.fv_candidates.end());

    // m_broadPhase = candidates.size();
    // std::cout << "m_broadPhase: " << m_broadPhase << std::endl;
}

/// @brief Compute a AABB for a vertex moving through time (i.e. temporal edge).
void calculate_vertex_extents(
    const VectorMax3d& vertex_t0,
    const VectorMax3d& vertex_t1,
    VectorMax3d& lower_bound,
    VectorMax3d& upper_bound)
{
    Eigen::MatrixXd points(2, vertex_t0.size());
    points.row(0) = vertex_t0;
    points.row(1) = vertex_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_edge_extents(
    const VectorMax3d& edge_vertex0_t0,
    const VectorMax3d& edge_vertex1_t0,
    const VectorMax3d& edge_vertex0_t1,
    const VectorMax3d& edge_vertex1_t1,
    VectorMax3d& lower_bound,
    VectorMax3d& upper_bound)
{
    Eigen::MatrixXd points(4, edge_vertex0_t0.size());
    points.row(0) = edge_vertex0_t0;
    points.row(1) = edge_vertex1_t0;
    points.row(2) = edge_vertex0_t1;
    points.row(3) = edge_vertex1_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_face_extents(
    const VectorMax3d& face_vertex0_t0,
    const VectorMax3d& face_vertex1_t0,
    const VectorMax3d& face_vertex2_t0,
    const VectorMax3d& face_vertex0_t1,
    const VectorMax3d& face_vertex1_t1,
    const VectorMax3d& face_vertex2_t1,
    VectorMax3d& lower_bound,
    VectorMax3d& upper_bound)
{
    Eigen::MatrixXd points(6, face_vertex0_t0.size());
    points.row(0) = face_vertex0_t0;
    points.row(1) = face_vertex1_t0;
    points.row(2) = face_vertex2_t0;
    points.row(3) = face_vertex0_t1;
    points.row(4) = face_vertex1_t1;
    points.row(5) = face_vertex2_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

ipc::AABB addFace(
    const VectorMax3d& face_vertex0_t0,
    const VectorMax3d& face_vertex1_t0,
    const VectorMax3d& face_vertex2_t0,
    const VectorMax3d& face_vertex0_t1,
    const VectorMax3d& face_vertex1_t1,
    const VectorMax3d& face_vertex2_t1,
    const long index,
    const double inflation_radius)
// ipc::AABB& aabb)
{
    VectorMax3d lower_bound, upper_bound;
    ipc::calculate_face_extents(
        face_vertex0_t0, face_vertex1_t0, face_vertex2_t0, //
        face_vertex0_t1, face_vertex1_t1, face_vertex2_t1, //
        lower_bound, upper_bound);
    // this->addElement(
    return ipc::AABB(
        lower_bound.array() - inflation_radius,
        upper_bound.array() + inflation_radius); //,
    // index, m_faceItems); // Faces have a positive id
}

void addFaces(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& faces,
    const double inflation_radius,
    std::vector<ipc::AABB>& face_aabbs)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
    for (long i = 0; i < faces.rows(); i++) {
        ipc::AABB aabb = addFace(
            vertices_t0.row(faces(i, 0)), vertices_t0.row(faces(i, 1)),
            vertices_t0.row(faces(i, 2)), vertices_t1.row(faces(i, 0)),
            vertices_t1.row(faces(i, 1)), vertices_t1.row(faces(i, 2)), i,
            inflation_radius);
        face_aabbs.push_back(aabb);
    }
}

ipc::AABB addVertex(
    const VectorMax3d& vertex_t0,
    const VectorMax3d& vertex_t1,
    const long index,
    const double inflation_radius)
{
    VectorMax3d lower_bound, upper_bound;
    ipc::calculate_vertex_extents(
        vertex_t0, vertex_t1, lower_bound, upper_bound);
    // this->addElement(
    return ipc::AABB(
        lower_bound.array() - inflation_radius,
        upper_bound.array() + inflation_radius); //,
    // index, m_vertexItems); // Vertices have a negative id
}

void addVertices(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double inflation_radius,
    std::vector<ipc::AABB>& vertex_aabbs)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
    for (long i = 0; i < vertices_t0.rows(); i++) {
        ipc::AABB aabb = addVertex(
            vertices_t0.row(i), vertices_t1.row(i), i, inflation_radius);
        vertex_aabbs.push_back(aabb);
    }
}

ipc::AABB addEdge(
    const VectorMax3d& edge_vertex0_t0,
    const VectorMax3d& edge_vertex1_t0,
    const VectorMax3d& edge_vertex0_t1,
    const VectorMax3d& edge_vertex1_t1,
    const long index,
    const double inflation_radius)
{
    VectorMax3d lower_bound, upper_bound;
    ipc::calculate_edge_extents(
        edge_vertex0_t0, edge_vertex1_t0, edge_vertex0_t1, edge_vertex1_t1,
        lower_bound, upper_bound);
    // this->addElement(
    return ipc::AABB(
        lower_bound.array() - inflation_radius,
        upper_bound.array() + inflation_radius); //,
    // index, m_edgeItems); // Edges have a positive id
}

void addEdges(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius,
    std::vector<ipc::AABB>& edge_aabbs)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
    for (long i = 0; i < edges.rows(); i++) {
        ipc::AABB aabb = addEdge(
            vertices_t0.row(edges(i, 0)), vertices_t0.row(edges(i, 1)),
            vertices_t1.row(edges(i, 0)), vertices_t1.row(edges(i, 1)), i,
            inflation_radius);
        edge_aabbs.push_back(aabb);
    }
}

void mesh_to_aabbs(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<ipc::AABB>& edge_aabbs,
    std::vector<ipc::AABB>& face_aabbs,
    std::vector<ipc::AABB>& vertex_aabbs,
    const double inflation_radius)
{
    ipc::addEdges(V0, V1, edges, inflation_radius, edge_aabbs);
    ipc::addFaces(V0, V1, faces, inflation_radius, face_aabbs);
    ipc::addVertices(V0, V1, inflation_radius, vertex_aabbs);
}

void combine_aabbs(
    std::vector<ipc::AABB>& edge_aabbs,
    std::vector<ipc::AABB>& face_aabbs,
    std::vector<ipc::AABB>& vertex_aabbs,
    std::vector<ipc::AABB>& aabbs)
{
    aabbs.clear();
    for (auto& aabb : vertex_aabbs)
        aabbs.emplace_back(aabb);
    for (auto& aabb : edge_aabbs)
        aabbs.emplace_back(aabb);
    for (auto& aabb : face_aabbs)
        aabbs.emplace_back(aabb);
}

void to_broadmark_aabbs(
    std::vector<ipc::AABB>& ipc_aabbs,
    std::vector<broadmark::Aabb>& broadmark_aabbs)
{
    broadmark_aabbs.clear();
    for (auto& aabb : ipc_aabbs) {
        auto min = aabb.min;
        auto max = aabb.max;
        broadmark_aabbs.emplace_back(
            Vec3(min[0], min[1], min[2]), Vec3(max[0], max[1], max[2]));
    }
}

void to_aabbs(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<broadmark::Aabb>& broadmark_aabbs)
{
    std::vector<ipc::AABB> edge_aabbs;
    std::vector<ipc::AABB> face_aabbs;
    std::vector<ipc::AABB> vertex_aabbs;
    std::vector<ipc::AABB> ipc_aabbs;
    ipc::mesh_to_aabbs(
        V0, V1, edges, faces, edge_aabbs, face_aabbs, vertex_aabbs);
    ipc::combine_aabbs(edge_aabbs, face_aabbs, vertex_aabbs, ipc_aabbs);
    ipc::to_broadmark_aabbs(ipc_aabbs, broadmark_aabbs);
}

// void toSimplexAabbs(
//     const Eigen::MatrixXd& V0,
//     const Eigen::MatrixXd& V1,
//     const Eigen::MatrixXi& E,
//     const Eigen::MatrixXi& F,
//     SimplexAabbs & S
// )
// {
//     std::vector<ipc::AABB> edge_aabbs;
//     std::vector<ipc::AABB> face_aabbs;
//     std::vector<ipc::AABB> vertex_aabbs;
//     std::vector<ipc::AABB> ipc_aabbs;
//     std::vector<Aabb> broadmark_aabbs;
//     broadphase::step_to_aabbs(V0,V1,E,F,edge_aabbs,face_aabbs,vertex_aabbs);
//     broadphase::combine_aabbs(edge_aabbs, face_aabbs, vertex_aabbs,
//     ipc_aabbs); broadphase::to_broadmark_aabbs(ipc_aabbs, broadmark_aabbs);

//     std::vector<simplexType> simplices;
//     simplices.resize(0);
//     std::map<simplexType, int> simplexFirstOccurs;

//     //adding vertices
//     std::vector<simplexType> vec;
//     vec.resize(vertex_aabbs.size());
//     vec.assign(vertex_aabbs.size(), simplexType::vertex);
//     simplexFirstOccurs[simplexType::vertex] = simplices.size();
//     simplices.insert(simplices.end(), vec.begin(), vec.end());

//     // adding edges
//     vec.resize(edge_aabbs.size());
//     vec.assign(edge_aabbs.size(), simplexType::edge);
//     simplexFirstOccurs[simplexType::edge] = simplices.size();
//     simplices.insert(simplices.end(), vec.begin(), vec.end());

//     //adding faces
//     vec.resize(face_aabbs.size());
//     vec.assign(face_aabbs.size(), simplexType::face);
//     simplexFirstOccurs[simplexType::face] = simplices.size();
//     simplices.insert(simplices.end(), vec.begin(), vec.end());

//     S.aabbs = broadmark_aabbs;
//     S.simplices = simplices;
//     S.simplexFirstOccurs = simplexFirstOccurs;
// }

broadmark::Aabb buildWorldAabb(const std::vector<broadmark::Aabb>& aabbs)
{

    Vec3 worldMin = aabbs[0].m_min;
    Vec3 worldMax = aabbs[0].m_max;
    for (auto& aabb : aabbs) {
        worldMin = Vec3::Min(worldMin, aabb.m_min);
        worldMax = Vec3::Max(worldMax, aabb.m_min);
    }

    broadmark::Aabb m_worldAabb = broadmark::Aabb(worldMin, worldMax);
    return m_worldAabb;
    // get m_min and m_max
    //  Vec3 Min(const Vec3& lhs, const Vec3& rhs)
}

void growAabbs(std::vector<broadmark::Aabb>& aabbs, const Vec3& amount)
{
    for (auto& aabb : aabbs) {
        aabb.Grow(amount);
    }
}

} // namespace ipc