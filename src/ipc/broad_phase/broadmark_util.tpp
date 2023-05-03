#pragma once
#include <tbb/tbb.h>

// #include <broadphase/interface.hpp>

// bool compare(ObjectPair a, ObjectPair b);

namespace ipc {

template <class T> Interface<T>::Interface() { algo = new T(); }

template <class T>
void Interface<T>::FilterOverlaps(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
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
    // std::cout << "V: " << num_vertices
    //           << " V+E: " << num_vertices + edges.rows()
    //           << " V+E+F: " << num_vertices + edges.rows() + faces.rows()
    //           << std::endl;
    // std::cout << "m_cache: " << algo->m_cache.m_overlaps.size() << std::endl;
    // for (size_t i = 0; i < algo->m_cache.m_overlaps.size(); i++) {
    //     int ai = algo->m_cache.m_overlaps[i].m_a->m_id;
    //     int bi = algo->m_cache.m_overlaps[i].m_b->m_id;
    //     addCandidate(ai, bi);
    // }
    // std::cout << "finished adding" << std::endl;
    // // remove duplicates
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

template <class T>
void Interface<T>::ConstructBoxes(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<broadmark::Aabb>& broadmark_aabbs)
{
    ipc::to_aabbs(V0, V1, edges, faces, broadmark_aabbs);
}

template <class T>
void Interface<T>::CalcOverlaps(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<broadmark::Aabb>& broadmark_aabbs,
    bool init)
{
    this->Clear();
    std::cout << "# candidates: " << candidates.size() << std::endl;

    SceneFrame sceneframe;
    sceneframe.m_aabbs = broadmark_aabbs.data();
    broadmark::Aabb m_worldAabb = ipc::buildWorldAabb(broadmark_aabbs);
    if (algo->m_settings.m_numberOfObjects != broadmark_aabbs.size() || init) {
        InflatedSettings settings = InflatedSettings();
        settings.m_numberOfObjects = broadmark_aabbs.size();
        settings.m_worldAabb = m_worldAabb;
        settings.m_vertices = V0.rows();
        settings.m_edges = edges.rows();
        settings.m_faces = faces.rows();
        tbb::global_control control(
            tbb::global_control::max_allowed_parallelism, 64);

        settings.m_numThreads = tbb::global_control::active_value(
            tbb::global_control::max_allowed_parallelism);
        std::cout << "max threads: " << settings.m_numThreads << std::endl;

        std::cout << "Initializing.." << std::endl;
        algo->Initialize(settings, sceneframe);
    }

    algo->UpdateObjects(sceneframe);
    algo->UpdateStructures();
    algo->CleanCache();
    algo->SearchOverlaps();

    m_Objects = algo->m_settings.m_numberOfObjects;
    std::cout << "m_Objects: " << m_Objects << std::endl;
}

} // namespace ipc