#include "esp_collisions.hpp"

#include "esp_collisions_builder.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/esp/quadrature_potential.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/profile_registry.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <numeric>
#include <stdexcept> // std::out_of_range
#include <utility>

namespace ipc {
namespace {
    template <typename Candidate>
    std::vector<VertexVertexCandidate>
    element_vertex_to_vertex_vertex_candidates(
        Eigen::ConstRef<Eigen::MatrixXi> elements,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::vector<Candidate>& candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<VertexVertexCandidate> vv_candidates;
        for (const auto& [ei, vi] : candidates) {
            for (int j = 0; j < elements.cols(); j++) {
                const int vj = elements(ei, j);
                if (is_active(point_point_distance(
                        vertices.row(vi), vertices.row(vj)))) {
                    vv_candidates.emplace_back(
                        std::min(vi, vj), std::max(vi, vj));
                }
            }
        }

        // Remove duplicates
        tbb::parallel_sort(vv_candidates.begin(), vv_candidates.end());
        vv_candidates.erase(
            std::unique(vv_candidates.begin(), vv_candidates.end()),
            vv_candidates.end());

        return vv_candidates;
    }

    std::vector<VertexVertexCandidate> face_vertex_to_vertex_vertex_candidates(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::vector<FaceVertexCandidate>& fv_candidates,
        const std::function<bool(double)>& is_active)
    {
        return element_vertex_to_vertex_vertex_candidates(
            mesh.faces(), vertices, fv_candidates, is_active);
    }

    std::vector<EdgeVertexCandidate> face_vertex_to_edge_vertex_candidates(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::vector<FaceVertexCandidate>& fv_candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        for (const auto& [fi, vi] : fv_candidates) {
            for (int j = 0; j < 3; j++) {
                const int ei = mesh.faces_to_edges()(fi, j);
                const int vj = mesh.edges()(ei, 0);
                const int vk = mesh.edges()(ei, 1);
                if (is_active(point_edge_distance(
                        vertices.row(vi), vertices.row(vj),
                        vertices.row(vk)))) {
                    ev_candidates.emplace_back(ei, vi);
                }
            }
        }

        // Remove duplicates
        tbb::parallel_sort(ev_candidates.begin(), ev_candidates.end());
        ev_candidates.erase(
            std::unique(ev_candidates.begin(), ev_candidates.end()),
            ev_candidates.end());

        return ev_candidates;
    }
} // namespace

void ESPCollisions::build(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters params)
{
    assert(vertices.rows() == mesh.num_vertices());

    IPC_PROFILE_SCOPE("ho.collision_build");
    clear();

    if (mesh.dim() == 2) {
        // Ensure candidate sets are populated (ev_set/ee_set/vv_set lookups
        // require them).
        const_cast<Candidates&>(candidates).convert_candidates_to_sets();

        tbb::enumerable_thread_specific<ESPCollisionsBuilder<2>> storage {
            ESPCollisionsBuilder<2>()
        };

        // Standard mode: loop over all edges with per-QP collision dicts.
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, mesh.num_edges()),
            [&](const tbb::blocked_range<size_t>& r) {
                ESPCollisionsBuilder<2>& local_storage = storage.local();
                local_storage.build_edge_collisions(
                    mesh, vertices, candidates, params, r.begin(), r.end());
            });
        ESPCollisionsBuilder<2>::merge(storage, *this);
    } else {
        // Compute vertex mask: which vertices to process.
        std::vector<bool> vertex_mask(mesh.num_vertices(), false);

        // Standard mode: only process vertices in face-vertex candidates.
        for (const auto& candidate : candidates.fv_candidates) {
            vertex_mask[candidate.vertex_id] = true;
        }

        std::vector<index_t> vertices_to_process;
        if (params.quad_order == 0) {
            vertices_to_process.reserve(mesh.num_vertices());
            for (int i = 0; i < mesh.num_vertices(); ++i) {
                if (vertex_mask[i]) {
                    vertices_to_process.push_back(i);
                }
            }
        }

        std::vector<index_t> faces_to_process;
        if (params.quad_order > 0) {
            faces_to_process.resize(mesh.num_faces());
            std::iota(faces_to_process.begin(), faces_to_process.end(), 0);
        }

        // create builder and parallel loops
        tbb::enumerable_thread_specific<QuadratureCollisionsBuilder> storage(
            QuadratureCollisionsBuilder(mesh, candidates, params));

        if (params.quad_order == 0) {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, vertices_to_process.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    QuadratureCollisionsBuilder& local_storage =
                        storage.local();
                    local_storage.build_vertex_collisions(
                        vertices, vertices_to_process, r.begin(), r.end());
                });
        }

        if (params.quad_order > 0) {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, faces_to_process.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    QuadratureCollisionsBuilder& local_storage =
                        storage.local();
                    local_storage.build_face_collisions(
                        vertices, faces_to_process, r.begin(), r.end());
                });
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, candidates.ee_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                QuadratureCollisionsBuilder& local_storage = storage.local();
                local_storage.build_edge_edge_collisions(
                    vertices, candidates.ee_candidates, r.begin(), r.end());
            });

        QuadratureCollisionsBuilder::merge(storage, *this);
    }
    m_candidates = candidates;

    size_t n_face_dicts = 0;
    size_t vert_pairs = 0, edge_pairs = 0, face_pairs = 0;
    for (const auto& cc : vertex_collisions) {
        vert_pairs += cc.second->size();
    }
    for (const auto& cc : edge_edge_collisions) {
        edge_pairs += cc.second->size();
    }
    for (const auto& fc : face_collisions) {
        n_face_dicts += fc.second.size();
        for (const auto& dict_ptr : fc.second) {
            face_pairs += dict_ptr->size();
        }
    }
    auto& reg = ProfileRegistry::instance();
    reg.add_value(
        "ho.collision_set.vertex_dicts",
        static_cast<double>(vertex_collisions.size()));
    reg.add_value(
        "ho.collision_set.edge_dicts",
        static_cast<double>(edge_edge_collisions.size()));
    reg.add_value(
        "ho.collision_set.face_dicts", static_cast<double>(n_face_dicts));
    reg.add_value(
        "ho.collision_set.vertex_pairs", static_cast<double>(vert_pairs));
    reg.add_value(
        "ho.collision_set.edge_pairs", static_cast<double>(edge_pairs));
    reg.add_value(
        "ho.collision_set.face_pairs", static_cast<double>(face_pairs));
    reg.add_value(
        "ho.collision_set.total_pairs",
        static_cast<double>(vert_pairs + edge_pairs + face_pairs));
    reg.add_value(
        "ho.candidates.fv",
        static_cast<double>(candidates.fv_candidates.size()));
    reg.add_value(
        "ho.candidates.ee",
        static_cast<double>(candidates.ee_candidates.size()));
    reg.add_value(
        "ho.candidates.total", static_cast<double>(candidates.size()));
}

std::unique_ptr<AdaptiveSupport> ESPCollisions::compute_adaptive_dhat(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters& params)
{
    return std::make_unique<AdaptiveSupport>(mesh, vertices, params);
}

void ESPCollisions::build(
    const Candidates& _candidates,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters params,
    const AdaptiveSupport* adaptive)
{
    adaptive_dhat =
        adaptive ? std::make_unique<AdaptiveSupport>(*adaptive) : nullptr;
    this->build(_candidates, mesh, vertices, params);
}

void ESPCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters params,
    BroadPhase* broad_phase)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius =
        params.dhat / 2; // TODO use dbar for EE collisions broad phase

    {
        IPC_PROFILE_SCOPE("ho.broad_phase");
        m_candidates.build(mesh, vertices, inflation_radius, broad_phase, true);
        {
            IPC_PROFILE_SCOPE("ho.convert_sets");
            m_candidates.convert_candidates_to_sets();
        }
    }

    this->build(m_candidates, mesh, vertices, params);
}

void ESPCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters params,
    const AdaptiveSupport* adaptive,
    BroadPhase* broad_phase)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius = params.dhat / 2;

    {
        IPC_PROFILE_SCOPE("ho.broad_phase");
        m_candidates.build(mesh, vertices, inflation_radius, broad_phase, true);
        {
            IPC_PROFILE_SCOPE("ho.convert_sets");
            m_candidates.convert_candidates_to_sets();
        }
    }

    this->build(m_candidates, mesh, vertices, params, adaptive);
}

// ============================================================================
size_t ESPCollisions::size() const
{
    size_t size = 0;
    for (const auto& cc : vertex_collisions) {
        size += cc.second->size();
    }
    for (const auto& cc : edge_edge_collisions) {
        size += cc.second->size();
    }
    for (const auto& cc : face_collisions) {
        for (const auto& dict_ptr : cc.second) {
            size += dict_ptr->size();
        }
    }
    for (const auto& cc : edge_collisions_2d) {
        for (const auto& dict_ptr : cc.second) {
            size += dict_ptr->size();
        }
    }
    return size;
}
bool ESPCollisions::empty() const
{
    return vertex_collisions.empty() && edge_edge_collisions.empty()
        && face_collisions.empty() && edge_collisions_2d.empty();
}
void ESPCollisions::clear()
{
    vertex_collisions.clear();
    edge_edge_collisions.clear();
    face_collisions.clear();
    edge_collisions_2d.clear();
}

std::string ESPCollisions::to_string(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters& params) const
{
    std::stringstream ss;

    for (const auto& ccs : vertex_collisions) {
        for (int i = 0; i < (*ccs.second).size(); i++) {
            const auto& cc = (*ccs.second)[i];
            ss << "\n";
            {
                ss << fmt::format(
                    "vert [{}]: ({} {}) weight {} dist sqr {} potential {} grad {}",
                    cc.name(), cc[0], cc[1], cc.weight,
                    cc.compute_distance(vertices),
                    cc(cc.dof(vertices), params, adaptive_dhat.get()),
                    cc.gradient(cc.dof(vertices), params, adaptive_dhat.get())
                        .norm());
            }
        }
    }
    for (const auto& ccs : edge_edge_collisions) {
        for (int i = 0; i < (*ccs.second).size(); i++) {
            const auto& cc = (*ccs.second)[i];
            ss << "\n";
            {
                ss << fmt::format(
                    "edge [{}]: ({} {}) ({} {}) weight {}", cc.name(),
                    ccs.first.first, ccs.first.second, cc[0], cc[1], cc.weight);
            }
        }
    }
    for (const auto& ccs : face_collisions) {
        for (const auto& dict_ptr : ccs.second) {
            for (int i = 0; i < dict_ptr->size(); i++) {
                const auto& cc = (*dict_ptr)[i];
                ss << "\n";
                {
                    ss << fmt::format(
                        "face [{}]: ({} {}) weight {} dist sqr {} potential {} grad {}",
                        cc.name(), cc[0], cc[1], cc.weight,
                        cc.compute_distance(vertices),
                        cc(cc.dof(vertices), params, adaptive_dhat.get()),
                        cc.gradient(
                              cc.dof(vertices), params, adaptive_dhat.get())
                            .norm());
                }
            }
        }
    }

    return ss.str();
}

// NOTE: Actually distance squared
double ESPCollisions::compute_minimum_distance(
    const CollisionMesh& mesh, Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (m_candidates.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, m_candidates.size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = m_candidates[i].compute_distance(
                    m_candidates[i].dof(vertices, edges, faces));

                local_min_dist = std::min(dist, local_min_dist);
            }
        });

    return storage.combine([](double a, double b) { return std::min(a, b); });
}

std::map<size_t, size_t> ESPCollisions::edge_id_count_distribution() const
{
    unordered_map<index_t, size_t> counts;
    for (const auto& [key, _] : edge_edge_collisions) {
        counts[key.first]++;
    }

    std::map<size_t, size_t> distribution;
    for (const auto& [_, count] : counts) {
        distribution[count]++;
    }
    return distribution;
}

Eigen::VectorXd ESPCollisions::edge_collision_counts(size_t num_edges) const
{
    Eigen::VectorXd counts = Eigen::VectorXd::Zero(num_edges);
    for (const auto& [key, _] : edge_edge_collisions) {
        counts(key.first)++;
    }
    return counts;
}

} // namespace ipc
