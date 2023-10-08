#include "friction_constraints.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

void FrictionConstraints::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const CollisionConstraints& contact_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu)
{
    assert(mus.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    const auto& C_vv = contact_constraint_set.vv_constraints;
    const auto& C_ev = contact_constraint_set.ev_constraints;
    const auto& C_ee = contact_constraint_set.ee_constraints;
    const auto& C_fv = contact_constraint_set.fv_constraints;
    auto& [FC_vv, FC_ev, FC_ee, FC_fv] = *this;

    FC_vv.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        FC_vv.emplace_back(
            c_vv, vertices, edges, faces, dhat, barrier_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

        FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i));
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(
            c_ev, vertices, edges, faces, dhat, barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
        FC_ev.back().mu = blend_mu(edge_mu, mus(vi));
    }

    FC_ee.reserve(C_ee.size());
    for (const auto& c_ee : C_ee) {
        const auto& [ea0i, ea1i, eb0i, eb1i] = c_ee.vertex_ids(edges, faces);
        const Eigen::Vector3d ea0 = vertices.row(ea0i);
        const Eigen::Vector3d ea1 = vertices.row(ea1i);
        const Eigen::Vector3d eb0 = vertices.row(eb0i);
        const Eigen::Vector3d eb1 = vertices.row(eb1i);

        // Skip EE constraints that are close to parallel
        if (edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1) < c_ee.eps_x) {
            continue;
        }

        FC_ee.emplace_back(
            c_ee, vertices, edges, faces, dhat, barrier_stiffness);

        double ea_mu =
            (mus(ea1i) - mus(ea0i)) * FC_ee.back().closest_point[0] + mus(ea0i);
        double eb_mu =
            (mus(eb1i) - mus(eb0i)) * FC_ee.back().closest_point[1] + mus(eb0i);
        FC_ee.back().mu = blend_mu(ea_mu, eb_mu);
    }

    FC_fv.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
        FC_fv.emplace_back(
            c_fv, vertices, edges, faces, dhat, barrier_stiffness);
        const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_ids(edges, faces);

        double face_mu = mus(f0i)
            + FC_fv.back().closest_point[0] * (mus(f1i) - mus(f0i))
            + FC_fv.back().closest_point[1] * (mus(f2i) - mus(f0i));
        FC_fv.back().mu = blend_mu(face_mu, mus(vi));
    }
}

double FrictionConstraints::compute_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& velocity,
    const double epsv) const
{
    assert(velocity.rows() == mesh.num_vertices());
    assert(epsv > 0);

    if (empty()) {
        return 0;
    }

    tbb::enumerable_thread_specific<double> storage(0);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](tbb::blocked_range<size_t> r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Quadrature weight is premultiplied by compute_potential
                local_potential += (*this)[i].compute_potential(
                    velocity, mesh.edges(), mesh.faces(), epsv);
            }
        });

    return storage.combine([](double a, double b) { return a + b; });
}

Eigen::VectorXd FrictionConstraints::compute_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& velocity,
    const double epsv) const
{
    const int dim = velocity.cols();
    const int ndof = velocity.size();

    if (empty()) {
        return Eigen::VectorXd::Zero(ndof);
    }
    assert(epsv > 0);

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(ndof));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& global_grad = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& constraint = (*this)[i];

                const VectorMax12d local_grad =
                    constraint.compute_potential_gradient(
                        velocity, mesh.edges(), mesh.faces(), epsv);

                const std::array<long, 4> vis =
                    constraint.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(
                    local_grad, vis, dim, global_grad);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

// ============================================================================

Eigen::SparseMatrix<double> FrictionConstraints::compute_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& velocity,
    const double epsv,
    const bool project_hessian_to_psd) const
{
    const int dim = velocity.cols();
    const int ndof = velocity.size();

    if (empty()) {
        return Eigen::SparseMatrix<double>(ndof, ndof);
    }
    assert(epsv > 0);

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& constraint = (*this)[i];

                const MatrixMax12d local_hess =
                    constraint.compute_potential_hessian(
                        velocity, mesh.edges(), mesh.faces(), epsv,
                        project_hessian_to_psd);

                const std::array<long, 4> vis =
                    constraint.vertex_ids(mesh.edges(), mesh.faces());

                local_hessian_to_global_triplets(
                    local_hess, vis, dim, hess_triplets);
            }
        });

    Eigen::SparseMatrix<double> hess(ndof, ndof);
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(ndof, ndof);
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

// ============================================================================

Eigen::VectorXd FrictionConstraints::compute_force(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& lagged_displacements,
    const Eigen::MatrixXd& velocities,
    const double dhat,
    const double barrier_stiffness,
    const double epsv,
    const double dmin,
    const bool no_mu) const
{
    if (empty()) {
        return Eigen::VectorXd::Zero(velocities.size());
    }
    assert(epsv > 0);

    int dim = velocities.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(velocities.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            Eigen::VectorXd& force = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& constraint = (*this)[i];

                const VectorMax12d local_force = constraint.compute_force(
                    rest_positions, lagged_displacements, velocities,
                    mesh.edges(), mesh.faces(), dhat, barrier_stiffness, epsv,
                    dmin, no_mu);

                const std::array<long, 4> vis =
                    constraint.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(local_force, vis, dim, force);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

// ============================================================================

Eigen::SparseMatrix<double> FrictionConstraints::compute_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& lagged_displacements,
    const Eigen::MatrixXd& velocities,
    const double dhat,
    const double barrier_stiffness,
    const double epsv,
    const FrictionConstraint::DiffWRT wrt,
    const double dmin) const
{
    if (empty()) {
        return Eigen::SparseMatrix<double>(
            velocities.size(), velocities.size());
    }
    assert(epsv > 0);

    int dim = velocities.cols();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& jac_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const FrictionConstraint& constraint = (*this)[i];

                const MatrixMax12d local_force_jacobian =
                    constraint.compute_force_jacobian(
                        rest_positions, lagged_displacements, velocities, edges,
                        faces, dhat, barrier_stiffness, epsv, wrt, dmin);

                const std::array<long, 4> vis =
                    constraint.vertex_ids(mesh.edges(), mesh.faces());

                local_hessian_to_global_triplets(
                    local_force_jacobian, vis, dim, jac_triplets);
            }
        });

    Eigen::SparseMatrix<double> jacobian(velocities.size(), velocities.size());
    for (const auto& local_jac_triplets : storage) {
        Eigen::SparseMatrix<double> local_jacobian(
            velocities.size(), velocities.size());
        local_jacobian.setFromTriplets(
            local_jac_triplets.begin(), local_jac_triplets.end());
        jacobian += local_jacobian;
    }

    // if wrt == X then compute ∇ₓ w(x)
    if (wrt == FrictionConstraint::DiffWRT::REST_POSITIONS) {
        for (int i = 0; i < this->size(); i++) {
            const FrictionConstraint& constraint = (*this)[i];
            assert(constraint.weight_gradient.size() == rest_positions.size());
            if (constraint.weight_gradient.size() != rest_positions.size()) {
                throw std::runtime_error(
                    "Shape derivative is not computed for friction constraint!");
            }

            VectorMax12d local_force = constraint.compute_force(
                rest_positions, lagged_displacements, velocities, edges, faces,
                dhat, barrier_stiffness, epsv, dmin);
            assert(constraint.weight != 0);
            local_force /= constraint.weight;

            Eigen::SparseVector<double> force(rest_positions.size());
            force.reserve(local_force.size());
            local_gradient_to_global_gradient(
                local_force, constraint.vertex_ids(edges, faces), dim, force);

            jacobian += force * constraint.weight_gradient.transpose();
        }
    }

    return jacobian;
}

// ============================================================================

size_t FrictionConstraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size();
}

bool FrictionConstraints::empty() const
{
    return vv_constraints.empty() && ev_constraints.empty()
        && ee_constraints.empty() && fv_constraints.empty();
}

void FrictionConstraints::clear()
{
    vv_constraints.clear();
    ev_constraints.clear();
    ee_constraints.clear();
    fv_constraints.clear();
}

FrictionConstraint& FrictionConstraints::operator[](size_t idx)
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
    throw std::out_of_range("Friction constraint index is out of range!");
}

const FrictionConstraint& FrictionConstraints::operator[](size_t idx) const
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
    throw std::out_of_range("Friction constraint index is out of range!");
}

} // namespace ipc
