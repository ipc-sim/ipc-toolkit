#include "collision_constraints.hpp"

#include <ipc/collisions/collision_constraints_builder.hpp>
// #include <ipc/utils/unordered_map_and_set.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

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
        CollisionConstraintsBuilder(*this));

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

    CollisionConstraintsBuilder::merge(storage, *this);

    for (size_t ci = 0; ci < size(); ci++) {
        CollisionConstraint& constraint = (*this)[ci];
        constraint.minimum_distance = dmin;
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
    Eigen::SparseMatrix<double> shape_derivative =
        compute_potential_hessian(mesh, vertices, dhat, false);

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = vertices.cols();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_triplets = storage.local();

            // for (size_t ci = 0; ci < constraint_set.size(); ci++) {
            for (size_t ci = r.begin(); ci < r.end(); ci++) {

                const CollisionConstraint& constraint = (*this)[ci];
                const Eigen::SparseVector<double>& weight_gradient =
                    constraint.weight_gradient;
                if (weight_gradient.size() != vertices.size()) {
                    throw std::runtime_error(
                        "Shape derivative is not computed for contact constraint!");
                }

                VectorMax12d local_barrier_grad =
                    constraint.compute_potential_gradient(
                        vertices, edges, faces, dhat);
                assert(constraint.weight != 0);
                local_barrier_grad.array() /= constraint.weight;

                const std::array<long, 4> ids =
                    constraint.vertex_ids(edges, faces);
                assert(local_barrier_grad.size() % dim == 0);
                const int n_verts = local_barrier_grad.size() / dim;
                assert(ids.size() >= n_verts); // Can be extra ids

                for (int i = 0; i < n_verts; i++) {
                    for (int d = 0; d < dim; d++) {
                        using Itr = Eigen::SparseVector<double>::InnerIterator;
                        for (Itr j(weight_gradient); j; ++j) {
                            local_triplets.emplace_back(
                                ids[i] * dim + d, j.index(),
                                local_barrier_grad[dim * i + d] * j.value());
                        }
                    }
                }
            }
        });

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

} // namespace ipc
