#include <ipc/friction/friction.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Sparse>

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/normal_force_magnitude.hpp>
#include <ipc/friction/relative_displacement.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/local_to_global.hpp>

namespace ipc {

void construct_friction_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    double mu,
    FrictionConstraints& friction_constraint_set)
{
    return construct_friction_constraint_set(
        mesh, V, contact_constraint_set, dhat, barrier_stiffness,
        Eigen::VectorXd::Constant(V.rows(), mu), friction_constraint_set);
}

void construct_friction_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    const Eigen::VectorXd& mus,
    FrictionConstraints& friction_constraint_set)
{
    return construct_friction_constraint_set(
        mesh, V, contact_constraint_set, dhat, barrier_stiffness, mus,
        [](double mu0, double mu1) { return (mu0 + mu1) / 2; },
        // [](double mu0, double mu1) { return mu0 * mu1; },
        // [](double mu0, double mu1) { return std::min(mu0, mu1); },
        // [](double mu0, double mu1) { return std::max(mu0, mu1); },
        friction_constraint_set);
}

void construct_friction_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu,
    FrictionConstraints& friction_constraint_set)
{
    assert(mus.size() == V.rows());

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    friction_constraint_set.clear();

    const auto& [C_vv, C_ev, C_ee, C_fv, _] = contact_constraint_set;
    auto& [FC_vv, FC_ev, FC_ee, FC_fv] = friction_constraint_set;

    FC_vv.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        FC_vv.emplace_back(c_vv, V, E, F, dhat, barrier_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_indices(E, F);

        FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i));
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(c_ev, V, E, F, dhat, barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_indices(E, F);

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
        FC_ev.back().mu = blend_mu(edge_mu, mus(vi));
    }

    FC_ee.reserve(C_ee.size());
    for (const auto& c_ee : C_ee) {
        const auto& [ea0i, ea1i, eb0i, eb1i] = c_ee.vertex_indices(E, F);
        const Eigen::Vector3d ea0 = V.row(ea0i), ea1 = V.row(ea1i),
                              eb0 = V.row(eb0i), eb1 = V.row(eb1i);

        // Skip EE constraints that are close to parallel
        if (edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1) < c_ee.eps_x) {
            continue;
        }

        FC_ee.emplace_back(c_ee, V, E, F, dhat, barrier_stiffness);

        double ea_mu =
            (mus(ea1i) - mus(ea0i)) * FC_ee.back().closest_point[0] + mus(ea0i);
        double eb_mu =
            (mus(eb1i) - mus(eb0i)) * FC_ee.back().closest_point[1] + mus(eb0i);
        FC_ee.back().mu = blend_mu(ea_mu, eb_mu);
    }

    FC_fv.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
        FC_fv.emplace_back(c_fv, V, E, F, dhat, barrier_stiffness);
        const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_indices(E, F);

        double face_mu = mus(f0i)
            + FC_fv.back().closest_point[0] * (mus(f1i) - mus(f0i))
            + FC_fv.back().closest_point[1] * (mus(f2i) - mus(f0i));
        FC_fv.back().mu = blend_mu(face_mu, mus(vi));
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// compute_friction_potential() in friction.tpp
//
///////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd compute_friction_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h)
{
    if (friction_constraint_set.empty()) {
        return Eigen::VectorXd::Zero(V0.size());
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    auto U = V1 - V0; // absolute linear dislacement of each point
    int dim = U.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(U.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), friction_constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_grad = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                local_gradient_to_global_gradient(
                    friction_constraint_set[i].compute_potential_gradient(
                        U, E, F, epsv_times_h),
                    friction_constraint_set[i].vertex_indices(E, F), dim,
                    local_grad);
            }
        });

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(U.size());
    for (const auto& local_grad : storage) {
        grad += local_grad;
    }
    return grad;
}

///////////////////////////////////////////////////////////////////////////////

Eigen::SparseMatrix<double> compute_friction_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h,
    bool project_hessian_to_psd)
{
    if (friction_constraint_set.empty()) {
        return Eigen::SparseMatrix<double>(V0.size(), V0.size());
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    auto U = V1 - V0; // absolute linear dislacement of each point
    int dim = U.cols();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), friction_constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                local_hessian_to_global_triplets(
                    friction_constraint_set[i].compute_potential_hessian(
                        U, E, F, epsv_times_h, project_hessian_to_psd),
                    friction_constraint_set[i].vertex_indices(E, F), dim,
                    local_hess_triplets);
            }
        });

    Eigen::SparseMatrix<double> hess(U.size(), U.size());
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(U.size(), U.size());
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

///////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd compute_friction_force(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const double dmin,
    const bool no_mu)
{
    if (friction_constraint_set.empty()) {
        return Eigen::VectorXd::Zero(U.size());
    }

    int dim = U.cols();
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(U.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), friction_constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_force = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                local_gradient_to_global_gradient(
                    friction_constraint_set[i].compute_force(
                        X, Ut, U, E, F, dhat, barrier_stiffness, epsv_times_h,
                        dmin, no_mu),
                    friction_constraint_set[i].vertex_indices(E, F), dim,
                    local_force);
            }
        });

    Eigen::VectorXd force = Eigen::VectorXd::Zero(U.size());
    for (const auto& local_force : storage) {
        force += local_force;
    }
    return force;
}

///////////////////////////////////////////////////////////////////////////////

Eigen::SparseMatrix<double> compute_friction_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const FrictionConstraint::DiffWRT wrt,
    const double dmin)
{
    if (friction_constraint_set.empty()) {
        return Eigen::SparseMatrix<double>(U.size(), U.size());
    }

    int dim = U.cols();
    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), friction_constraint_set.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& local_jac_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                local_hessian_to_global_triplets(
                    friction_constraint_set[i].compute_force_jacobian(
                        X, Ut, U, E, F, dhat, barrier_stiffness, epsv_times_h,
                        wrt, dmin),
                    friction_constraint_set[i].vertex_indices(E, F), dim,
                    local_jac_triplets);
            }
        });

    Eigen::SparseMatrix<double> jacobian(U.size(), U.size());
    for (const auto& local_jac_triplets : storage) {
        Eigen::SparseMatrix<double> local_jacobian(U.size(), U.size());
        local_jacobian.setFromTriplets(
            local_jac_triplets.begin(), local_jac_triplets.end());
        jacobian += local_jacobian;
    }

#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
    // if wrt == X then compute ∇ₓ w(x)
    if (wrt == FrictionConstraint::DiffWRT::X) {
        for (int i = 0; i < friction_constraint_set.size(); i++) {
            const FrictionConstraint& constraint = friction_constraint_set[i];
            assert(constraint.weight_gradient.size() == X.size());

            VectorMax12d local_force = constraint.compute_force(
                X, Ut, U, E, F, dhat, barrier_stiffness, epsv_times_h, dmin);
            assert(constraint.weight != 0);
            local_force.array() /= constraint.weight;

            Eigen::SparseVector<double> force(X.size());
            force.reserve(local_force.size());
            local_gradient_to_global_gradient(
                local_force, constraint.vertex_indices(E, F), dim, force);

            jacobian += force * constraint.weight_gradient.transpose();
        }
    }
#endif

    return jacobian;
}

} // namespace ipc
