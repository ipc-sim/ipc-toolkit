#include "esp_collisions.hpp"

#include <ipc/esp/esp_collision_maps.hpp>
#include <ipc/esp/esp_collisions_builder.hpp>
#include <ipc/esp/quadrature_potential.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <numeric>
#include <stdexcept> // std::out_of_range, std::invalid_argument
#include <utility>

namespace ipc {

ESPCollisions::ESPCollisions() : m_maps(std::make_unique<Maps>()) { }

ESPCollisions::~ESPCollisions() = default;

void ESPCollisions::build(
    const Candidates& _candidates,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters& params)
{
    assert(static_cast<size_t>(vertices.rows()) == mesh.num_vertices());

    if (mesh.dim() == 2 && params.quad_order < 1) {
        log_and_throw_error<std::invalid_argument>(
            "quad_order {} must be at least 1 in 2D.", params.quad_order);
    }

    // In 3D, ESP corrects each point's face contacts by inclusion-exclusion
    // over the nearby edges and vertices, which it finds through the
    // vertex-vertex, edge-vertex, edge-face, and face-face candidates. Without
    // them, the energy is still finite but wrong (over-counted, or even
    // attractive), so we reject candidates that were not built with all_types.
    if (mesh.dim() == 3 && !_candidates.has_all_types) {
        log_and_throw_error<std::invalid_argument>(
            "ESPCollisions::build: 3D candidates must be built with all_types = true");
    }

    clear();

    if (&_candidates != &m_candidates) {
        m_candidates = _candidates;
    }
    m_candidates.convert_candidates_to_sets();
    const Candidates& candidates = m_candidates;

    // Record the quadrature so that evaluating the set with different
    // parameters fails loudly in check_quadrature() instead of misreading the
    // per-point collision sets.
    m_dim = mesh.dim();
    m_quad_order = params.quad_order;
    m_face_quad_rule = params.quad_rule();
    m_lattice_spacing = params.lattice_spacing();
    m_max_lattice_subdivisions = params.max_lattice_subdivisions();

    if (mesh.dim() == 2) {
        tbb::enumerable_thread_specific<ESPCollisionsBuilder<2>> storage {
            ESPCollisionsBuilder<2>()
        };
        // Built once (not per thread) since it computes the mesh's
        // inclusion-exclusion weights.
        const PointPotential point_potential(mesh, candidates, params);

        // Standard mode: loop over all edges with per-QP collision dicts.
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, mesh.num_edges()),
            [&](const tbb::blocked_range<size_t>& r) {
                ESPCollisionsBuilder<2>& local_storage = storage.local();
                local_storage.build_edge_collisions(
                    point_potential, vertices, r.begin(), r.end());
            });
        ESPCollisionsBuilder<2>::merge(storage, *this);
    } else {
        // A vertex sample needs a collision dict only if something is near
        // it, so we mark every vertex of a vertex candidate.
        std::vector<bool> vertex_mask(mesh.num_vertices(), false);
        for (const auto& candidate : candidates.fv_candidates) {
            vertex_mask[candidate.vertex_id] = true;
        }
        for (const auto& candidate : candidates.ev_candidates) {
            vertex_mask[candidate.vertex_id] = true;
        }
        for (const auto& candidate : candidates.vv_candidates) {
            vertex_mask[candidate.vertex0_id] = true;
            vertex_mask[candidate.vertex1_id] = true;
        }

        // Every vertex is a sample: a face rule's corner points use the
        // shared vertex samples (ESPParameters::corner_weights()), and edges
        // in no face and isolated vertices integrate over their vertices.
        // Only the rule's other points need per-face samples.
        const bool use_face_quadrature = params.has_non_corner_points();

        std::vector<index_t> vertices_to_process;
        vertices_to_process.reserve(mesh.num_vertices());
        for (int i = 0; i < mesh.num_vertices(); ++i) {
            if (vertex_mask[i]) {
                vertices_to_process.push_back(i);
            }
        }

        std::vector<index_t> faces_to_process;
        if (use_face_quadrature) {
            faces_to_process.resize(mesh.num_faces());
            std::iota(faces_to_process.begin(), faces_to_process.end(), 0);
        }

        // create builder and parallel loops
        tbb::enumerable_thread_specific<QuadratureCollisionsBuilder> storage(
            QuadratureCollisionsBuilder(mesh, candidates, params));

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, vertices_to_process.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                QuadratureCollisionsBuilder& local_storage = storage.local();
                local_storage.build_vertex_collisions(
                    vertices, vertices_to_process, r.begin(), r.end());
            });

        if (use_face_quadrature) {
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
}

void ESPCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters& params,
    BroadPhase* broad_phase)
{
    assert(static_cast<size_t>(vertices.rows()) == mesh.num_vertices());

    double inflation_radius =
        params.dhat / 2; // TODO use dbar for EE collisions broad phase

    m_candidates.build(mesh, vertices, inflation_radius, broad_phase, true);

    this->build(m_candidates, mesh, vertices, params);
}

// ============================================================================
size_t ESPCollisions::size() const
{
    size_t size = 0;
    for (const auto& cc : m_maps->vertex_collisions) {
        size += cc.second->size();
    }
    for (const auto& cc : m_maps->edge_edge_collisions) {
        size += cc.second->size();
    }
    for (const auto& cc : m_maps->face_collisions) {
        for (const auto& dict_ptr : cc.second) {
            size += dict_ptr->size();
        }
    }
    for (const auto& cc : m_maps->edge_collisions_2d) {
        for (const auto& dict_ptr : cc.second) {
            size += dict_ptr->size();
        }
    }
    return size;
}
bool ESPCollisions::empty() const
{
    return m_maps->vertex_collisions.empty()
        && m_maps->edge_edge_collisions.empty()
        && m_maps->face_collisions.empty()
        && m_maps->edge_collisions_2d.empty();
}
void ESPCollisions::clear()
{
    m_maps->vertex_collisions.clear();
    m_maps->edge_edge_collisions.clear();
    m_maps->face_collisions.clear();
    m_maps->edge_collisions_2d.clear();
    num_quadrature_collision_pairs = 0;
    m_dim = 0;
    m_quad_order = 0;
    m_face_quad_rule.clear();
    m_lattice_spacing = 0;
    m_max_lattice_subdivisions = 0;
}

void ESPCollisions::check_quadrature(const ESPParameters& params) const
{
    if (m_dim == 2 && params.quad_order != m_quad_order) {
        log_and_throw_error<std::invalid_argument>(
            "ESPCollisions were built with quad_order {}, but the "
            "parameters use quad_order {}. Use the same ESPParameters to "
            "build the collisions and to evaluate them.",
            m_quad_order, params.quad_order);
    }

    if (m_dim == 3) {
        const FaceQuadRule& rule = params.quad_rule();
        const bool same_rule = rule.size() == m_face_quad_rule.size()
            && std::equal(rule.begin(), rule.end(), m_face_quad_rule.begin(),
                          [](const FaceQuadPoint& a, const FaceQuadPoint& b) {
                              return a.lambda == b.lambda
                                  && a.weight == b.weight;
                          });
        if (!same_rule) {
            log_and_throw_error<std::invalid_argument>(
                "ESPCollisions were built with a {}-point face quadrature "
                "rule, but the parameters use a different {}-point rule. "
                "Use the same ESPParameters to build the collisions and "
                "to evaluate them, and call set_quad_rule() before "
                "constructing ESPPotential, which copies the parameters.",
                m_face_quad_rule.size(), rule.size());
        }
        if (params.lattice_spacing() != m_lattice_spacing
            || params.max_lattice_subdivisions()
                != m_max_lattice_subdivisions) {
            log_and_throw_error<std::invalid_argument>(
                "ESPCollisions were built with lattice spacing {} (at most "
                "{} subdivisions), but the parameters use {} (at most {}). "
                "Use the same ESPParameters to build the collisions and to "
                "evaluate them.",
                m_lattice_spacing, m_max_lattice_subdivisions,
                params.lattice_spacing(), params.max_lattice_subdivisions());
        }
    }
}

std::string ESPCollisions::to_string(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPParameters& params) const
{
    std::stringstream ss;

    for (const auto& ccs : m_maps->vertex_collisions) {
        for (int i = 0; i < (*ccs.second).size(); i++) {
            const auto& cc = (*ccs.second)[i];
            ss << "\n";
            {
                ss << fmt::format(
                    "vert [{}]: ({} {}) weight {} dist sqr {} potential {} grad {}",
                    cc.name(), cc[0], cc[1], cc.weight,
                    cc.compute_distance(vertices), cc(cc.dof(vertices), params),
                    cc.gradient(cc.dof(vertices), params).norm());
            }
        }
    }
    for (const auto& ccs : m_maps->edge_edge_collisions) {
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
    // Face collisions involve the quadrature point, a virtual vertex that is
    // not a row of vertices, so only their ids and weights are printed.
    for (const auto& ccs : m_maps->face_collisions) {
        for (size_t qi = 0; qi < ccs.second.size(); qi++) {
            const auto& dict = *ccs.second[qi];
            for (int i = 0; i < dict.size(); i++) {
                const auto& cc = dict[i];
                ss << "\n";
                ss << fmt::format(
                    "face [{}]: {} point {} ({} {}) weight {}", cc.name(),
                    ccs.first, qi, cc[0], cc[1], cc.weight);
            }
        }
    }
    // Same for the 2D edge quadrature points.
    for (const auto& ccs : m_maps->edge_collisions_2d) {
        for (size_t qi = 0; qi < ccs.second.size(); qi++) {
            const auto& dict = *ccs.second[qi];
            for (int i = 0; i < dict.size(); i++) {
                const auto& cc = dict[i];
                ss << "\n";
                ss << fmt::format(
                    "edge [{}]: {} point {} ({} {}) weight {}", cc.name(),
                    ccs.first, qi, cc[0], cc[1], cc.weight);
            }
        }
    }

    return ss.str();
}

// NOTE: Actually distance squared
double ESPCollisions::compute_minimum_distance(
    const CollisionMesh& mesh, Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(static_cast<size_t>(vertices.rows()) == mesh.num_vertices());

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

} // namespace ipc
