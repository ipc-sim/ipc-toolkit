#include "collision_constraints.hpp"

#include <ipc/collisions/collision_constraints_builder.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

namespace {
    /// @brief Convert element-vertex candidates to vertex-vertex candidates
    /// @param elements Elements matrix of the mesh
    /// @param vertices Vertex positions of the mesh
    /// @param ev_candidates Element-vertex candidates
    /// @param is_active Function to determine if a candidate is active
    /// @return Vertex-vertex candidates
    template <typename Candidate>
    std::vector<VertexVertexCandidate>
    element_vertex_to_vertex_vertex_candidates(
        const Eigen::MatrixXi& elements,
        const Eigen::MatrixXd& vertices,
        const std::vector<Candidate>& candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<VertexVertexCandidate> vv_candidates;
        for (const auto& [ei, vi] : candidates) {
            for (int j = 0; j < elements.cols(); j++) {
                const int vj = elements(ei, j);
                if (is_active(point_point_distance(
                        vertices.row(vi), vertices.row(vj)))) {
                    vv_candidates.emplace_back(vi, vj);
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

    std::vector<VertexVertexCandidate> edge_vertex_to_vertex_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& ev_candidates,
        const std::function<bool(double)>& is_active)
    {
        return element_vertex_to_vertex_vertex_candidates(
            mesh.edges(), vertices, ev_candidates, is_active);
    }

    std::vector<VertexVertexCandidate> face_vertex_to_vertex_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& fv_candidates,
        const std::function<bool(double)>& is_active)
    {
        return element_vertex_to_vertex_vertex_candidates(
            mesh.faces(), vertices, fv_candidates, is_active);
    }

    std::vector<EdgeVertexCandidate> face_vertex_to_edge_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
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
                        vertices.row(vi), //
                        vertices.row(vj), vertices.row(vk)))) {
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

    std::vector<EdgeVertexCandidate> edge_edge_to_edge_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& ee_candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        for (const EdgeEdgeCandidate& ee : ee_candidates) {
            for (int i = 0; i < 2; i++) {
                const int ei = i == 0 ? ee.edge0_id : ee.edge1_id;
                const int ej = i == 0 ? ee.edge1_id : ee.edge0_id;

                const int ei0 = mesh.edges()(ei, 0);
                const int ei1 = mesh.edges()(ei, 1);

                for (int j = 0; j < 2; j++) {
                    const int vj = mesh.edges()(ej, j);
                    if (is_active(point_edge_distance(
                            vertices.row(vj), //
                            vertices.row(ei0), vertices.row(ei1)))) {
                        ev_candidates.emplace_back(ei, vj);
                    }
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

void CollisionConstraints::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const double dmin,
    const BroadPhaseMethod broad_phase_method)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius = (dhat + dmin) / 2;

    Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase_method);

    this->build(candidates, mesh, vertices, dhat, dmin);
}

void CollisionConstraints::build(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const double dmin)
{
    assert(vertices.rows() == mesh.num_vertices());

    clear();

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.
    const double offset_sqr = (dmin + dhat) * (dmin + dhat);
    auto is_active = [&](double distance_sqr) {
        return distance_sqr < offset_sqr;
    };

    tbb::enumerable_thread_specific<CollisionConstraintsBuilder> storage(
        use_convergent_formulation(), are_shape_derivatives_enabled());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.vv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_vertex_vertex_constraints(
                mesh, vertices, candidates.vv_candidates, is_active, r.begin(),
                r.end());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ev_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_edge_vertex_constraints(
                mesh, vertices, candidates.ev_candidates, is_active, r.begin(),
                r.end());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ee_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_edge_edge_constraints(
                mesh, vertices, candidates.ee_candidates, is_active, r.begin(),
                r.end());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.fv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_face_vertex_constraints(
                mesh, vertices, candidates.fv_candidates, is_active, r.begin(),
                r.end());
        });

    if (use_convergent_formulation()) {
        if (candidates.ev_candidates.size() > 0) {
            // Convert edge-vertex to vertex-vertex
            const std::vector<VertexVertexCandidate> vv_candidates =
                edge_vertex_to_vertex_vertex_candidates(
                    mesh, vertices, candidates.ev_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), vv_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_edge_vertex_negative_vertex_vertex_constraints(
                            mesh, vertices, vv_candidates, r.begin(), r.end());
                });
        }

        if (candidates.ee_candidates.size() > 0) {
            // Convert edge-edge to edge-vertex
            const auto ev_candidates = edge_edge_to_edge_vertex_candidates(
                mesh, vertices, candidates.ee_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), ev_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_edge_edge_negative_edge_vertex_constraints(
                            mesh, vertices, ev_candidates, r.begin(), r.end());
                });
        }

        if (candidates.fv_candidates.size() > 0) {
            // Convert face-vertex to edge-vertex
            const std::vector<EdgeVertexCandidate> ev_candidates =
                face_vertex_to_edge_vertex_candidates(
                    mesh, vertices, candidates.fv_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), ev_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_face_vertex_negative_edge_vertex_constraints(
                            mesh, vertices, ev_candidates, r.begin(), r.end());
                });

            // Convert face-vertex to vertex-vertex
            const std::vector<VertexVertexCandidate> vv_candidates =
                face_vertex_to_vertex_vertex_candidates(
                    mesh, vertices, candidates.fv_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), vv_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_face_vertex_positive_vertex_vertex_constraints(
                            mesh, vertices, vv_candidates, r.begin(), r.end());
                });
        }
    }

    // -------------------------------------------------------------------------

    CollisionConstraintsBuilder::merge(storage, *this);

    // logger().debug(to_string(mesh, vertices));

    for (size_t ci = 0; ci < size(); ci++) {
        CollisionConstraint& constraint = (*this)[ci];
        constraint.dmin = dmin;
    }

    if (use_convergent_formulation()) {
        // NOTE: When using the convergent formulation we want the barrier to
        // have units of Pa⋅m, so κ gets units of Pa and the barrier function
        // should have units of m. See notebooks/physical_barrier.ipynb for more
        // details.
        const double barrier_to_physical_barrier_divisor =
            dhat * std::pow(dhat + 2 * dmin, 2);

        for (size_t ci = 0; ci < size(); ci++) {
            CollisionConstraint& constraint = (*this)[ci];
            constraint.weight /= barrier_to_physical_barrier_divisor;
            if (are_shape_derivatives_enabled()) {
                constraint.weight_gradient /=
                    barrier_to_physical_barrier_divisor;
            }
        }
    }
}

void CollisionConstraints::set_use_convergent_formulation(
    const bool use_convergent_formulation)
{
    if (!empty()
        && use_convergent_formulation != m_use_convergent_formulation) {
        logger().warn(
            "Setting use_convergent_formulation after building constraints. "
            "Re-build constraints for this to have an effect.");
    }

    m_use_convergent_formulation = use_convergent_formulation;
}

void CollisionConstraints::set_are_shape_derivatives_enabled(
    const bool are_shape_derivatives_enabled)
{
    if (!empty()
        && are_shape_derivatives_enabled != m_are_shape_derivatives_enabled) {
        logger().warn(
            "Setting enable_shape_derivatives after building constraints. "
            "Re-build constraints for this to have an effect.");
    }

    m_are_shape_derivatives_enabled = are_shape_derivatives_enabled;
}

// ============================================================================

double CollisionConstraints::compute_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat) const
{
    assert(vertices.rows() == mesh.num_vertices());
    assert(dhat > 0);

    if (empty()) {
        return 0;
    }

    tbb::enumerable_thread_specific<double> storage(0);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Quadrature weight is premultiplied by compute_potential
                local_potential += (*this)[i].compute_potential(
                    vertices, mesh.edges(), mesh.faces(), dhat);
            }
        });

    return storage.combine([](double a, double b) { return a + b; });
}

Eigen::VectorXd CollisionConstraints::compute_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (empty()) {
        return Eigen::VectorXd::Zero(vertices.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    int dim = vertices.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(vertices.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_grad = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                local_gradient_to_global_gradient(
                    (*this)[i].compute_potential_gradient(
                        vertices, edges, faces, dhat),
                    (*this)[i].vertex_ids(edges, faces), dim, local_grad);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

Eigen::SparseMatrix<double> CollisionConstraints::compute_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const bool project_hessian_to_psd) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (empty()) {
        return Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = vertices.cols();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                local_hessian_to_global_triplets(
                    (*this)[i].compute_potential_hessian(
                        vertices, edges, faces, dhat, project_hessian_to_psd),
                    (*this)[i].vertex_ids(edges, faces), dim,
                    local_hess_triplets);
            }
        });

    Eigen::SparseMatrix<double> hess(vertices.size(), vertices.size());
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(
            vertices.size(), vertices.size());
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

// ============================================================================

Eigen::SparseMatrix<double> CollisionConstraints::compute_shape_derivative(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (empty()) {
        return Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
    }

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                (*this)[i].compute_shape_derivative(
                    mesh.rest_positions(), vertices, mesh.edges(), mesh.faces(),
                    dhat, local_triplets);
            }
        });

    Eigen::SparseMatrix<double> shape_derivative(
        vertices.size(), vertices.size());
    for (const auto& local_triplets : storage) {
        Eigen::SparseMatrix<double> local_shape_derivative(
            vertices.size(), vertices.size());
        local_shape_derivative.setFromTriplets(
            local_triplets.begin(), local_triplets.end());
        shape_derivative += local_shape_derivative;
    }
    return shape_derivative;
}

// ============================================================================

// NOTE: Actually distance squared
double CollisionConstraints::compute_minimum_distance(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    // Do a single block range over all constraint vectors
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist =
                    (*this)[i].compute_distance(vertices, edges, faces);

                if (dist < local_min_dist) {
                    local_min_dist = dist;
                }
            }
        });

    return storage.combine([](double a, double b) { return std::min(a, b); });
}

// ============================================================================

size_t CollisionConstraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size() + pv_constraints.size();
}

bool CollisionConstraints::empty() const
{
    return vv_constraints.empty() && ev_constraints.empty()
        && ee_constraints.empty() && fv_constraints.empty()
        && pv_constraints.empty();
}

void CollisionConstraints::clear()
{
    vv_constraints.clear();
    ev_constraints.clear();
    ee_constraints.clear();
    fv_constraints.clear();
    pv_constraints.clear();
}

CollisionConstraint& CollisionConstraints::operator[](size_t idx)
{
    if (idx < vv_constraints.size()) {
        return vv_constraints[idx];
    }
    idx -= vv_constraints.size();
    if (idx < ev_constraints.size()) {
        return ev_constraints[idx];
    }
    idx -= ev_constraints.size();
    if (idx < ee_constraints.size()) {
        return ee_constraints[idx];
    }
    idx -= ee_constraints.size();
    if (idx < fv_constraints.size()) {
        return fv_constraints[idx];
    }
    idx -= fv_constraints.size();
    if (idx < pv_constraints.size()) {
        return pv_constraints[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

const CollisionConstraint& CollisionConstraints::operator[](size_t idx) const
{
    if (idx < vv_constraints.size()) {
        return vv_constraints[idx];
    }
    idx -= vv_constraints.size();
    if (idx < ev_constraints.size()) {
        return ev_constraints[idx];
    }
    idx -= ev_constraints.size();
    if (idx < ee_constraints.size()) {
        return ee_constraints[idx];
    }
    idx -= ee_constraints.size();
    if (idx < fv_constraints.size()) {
        return fv_constraints[idx];
    }
    idx -= fv_constraints.size();
    if (idx < pv_constraints.size()) {
        return pv_constraints[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

bool CollisionConstraints::is_vertex_vertex(size_t idx) const
{
    return idx < vv_constraints.size();
}

bool CollisionConstraints::is_edge_vertex(size_t idx) const
{
    return idx >= vv_constraints.size()
        && idx < vv_constraints.size() + ev_constraints.size();
}

bool CollisionConstraints::is_edge_edge(size_t idx) const
{
    return idx >= vv_constraints.size() + ev_constraints.size()
        && idx
        < vv_constraints.size() + ev_constraints.size() + ee_constraints.size();
}

bool CollisionConstraints::is_face_vertex(size_t idx) const
{
    return idx
        >= vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        && idx < vv_constraints.size() + ev_constraints.size()
            + ee_constraints.size() + fv_constraints.size();
}

bool CollisionConstraints::is_plane_vertex(size_t idx) const
{
    return idx >= vv_constraints.size() + ev_constraints.size()
            + ee_constraints.size() + fv_constraints.size()
        && idx < vv_constraints.size() + ev_constraints.size()
            + ee_constraints.size() + fv_constraints.size()
            + pv_constraints.size();
}

std::string CollisionConstraints::to_string(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    std::stringstream ss;
    for (const auto& vv : vv_constraints) {
        ss << "\n"
           << fmt::format(
                  "vv: {} {}, w: {:g}, d: {:g}", vv.vertex0_id, vv.vertex1_id,
                  vv.weight,
                  point_point_distance(
                      vertices.row(vv.vertex0_id),
                      vertices.row(vv.vertex1_id)));
    }
    for (const auto& ev : ev_constraints) {
        ss << "\n"
           << fmt::format(
                  "ev: {}=({}, {}) {}, w: {:g}, d: {:g}", ev.edge_id,
                  mesh.edges()(ev.edge_id, 0), mesh.edges()(ev.edge_id, 1),
                  ev.vertex_id, ev.weight,
                  point_line_distance(
                      vertices.row(ev.vertex_id),
                      vertices.row(mesh.edges()(ev.edge_id, 0)),
                      vertices.row(mesh.edges()(ev.edge_id, 1))));
    }
    for (const auto& ee : ee_constraints) {
        ss << "\n"
           << fmt::format(
                  "ee: {}=({}, {}) {}=({}, {}), w: {:g}, dtype: {}, d: {:g}",
                  ee.edge0_id, mesh.edges()(ee.edge0_id, 0),
                  mesh.edges()(ee.edge0_id, 1), ee.edge1_id,
                  mesh.edges()(ee.edge1_id, 0), mesh.edges()(ee.edge1_id, 1),
                  ee.weight, int(ee.dtype),
                  edge_edge_distance(
                      vertices.row(mesh.edges()(ee.edge0_id, 0)),
                      vertices.row(mesh.edges()(ee.edge0_id, 1)),
                      vertices.row(mesh.edges()(ee.edge1_id, 0)),
                      vertices.row(mesh.edges()(ee.edge1_id, 1)), ee.dtype));
    }
    for (const auto& fv : fv_constraints) {
        ss << "\n"
           << fmt::format(
                  "fv: {}=({}, {}, {}) {}, w: {:g}, d: {:g}", fv.face_id,
                  mesh.faces()(fv.face_id, 0), mesh.faces()(fv.face_id, 1),
                  mesh.faces()(fv.face_id, 2), fv.vertex_id, fv.weight,
                  point_plane_distance(
                      vertices.row(fv.vertex_id),
                      vertices.row(mesh.faces()(fv.face_id, 0)),
                      vertices.row(mesh.faces()(fv.face_id, 1)),
                      vertices.row(mesh.faces()(fv.face_id, 2))));
    }
    return ss.str();
}

} // namespace ipc
