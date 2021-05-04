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
#include <ipc/friction/relative_displacement.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/local_to_global.hpp>

namespace ipc {

double compute_normal_force_magnitude(
    double distance_squared,
    double dhat,
    double barrier_stiffness,
    double dmin = 0)
{
    double grad_b = barrier_gradient(
        distance_squared - dmin * dmin, 2 * dmin * dhat + dhat * dhat);
    grad_b *= barrier_stiffness;
    return -grad_b * 2 * sqrt(distance_squared); // / (h * h) eliminated here
}

void construct_friction_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    double mu,
    FrictionConstraints& friction_constraint_set)
{
    friction_constraint_set.vv_constraints.reserve(
        contact_constraint_set.vv_constraints.size());
    for (const auto& vv_constraint : contact_constraint_set.vv_constraints) {
        const auto& p0 = V.row(vv_constraint.vertex0_index);
        const auto& p1 = V.row(vv_constraint.vertex1_index);

        friction_constraint_set.vv_constraints.emplace_back(vv_constraint);
        // Do not initialize closest point because it is trivial
        friction_constraint_set.vv_constraints.back().tangent_basis =
            point_point_tangent_basis(p0, p1);
        friction_constraint_set.vv_constraints.back().normal_force_magnitude =
            compute_normal_force_magnitude(
                point_point_distance(p0, p1), dhat, barrier_stiffness,
                vv_constraint.minimum_distance);
        friction_constraint_set.vv_constraints.back().mu = mu;
    }

    friction_constraint_set.ev_constraints.reserve(
        contact_constraint_set.ev_constraints.size());
    for (const auto& ev_constraint : contact_constraint_set.ev_constraints) {
        const auto& p = V.row(ev_constraint.vertex_index);
        const auto& e0 = V.row(E(ev_constraint.edge_index, 0));
        const auto& e1 = V.row(E(ev_constraint.edge_index, 1));

        friction_constraint_set.ev_constraints.emplace_back(ev_constraint);

        friction_constraint_set.ev_constraints.back().closest_point.resize(1);
        friction_constraint_set.ev_constraints.back().closest_point[0] =
            point_edge_closest_point(p, e0, e1);
        friction_constraint_set.ev_constraints.back().tangent_basis =
            point_edge_tangent_basis(p, e0, e1);
        friction_constraint_set.ev_constraints.back().normal_force_magnitude =
            compute_normal_force_magnitude(
                point_edge_distance(p, e0, e1, PointEdgeDistanceType::P_E),
                dhat, barrier_stiffness, ev_constraint.minimum_distance);
        friction_constraint_set.ev_constraints.back().mu = mu;
    }

    friction_constraint_set.ee_constraints.reserve(
        contact_constraint_set.ee_constraints.size());
    for (const auto& ee_constraint : contact_constraint_set.ee_constraints) {
        const auto& ea0 = V.row(E(ee_constraint.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_constraint.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_constraint.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_constraint.edge1_index, 1));

        // Skip EE constraints that are close to parallel
        if (edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1)
            < ee_constraint.eps_x) {
            continue;
        }

        friction_constraint_set.ee_constraints.emplace_back(ee_constraint);

        friction_constraint_set.ee_constraints.back().closest_point =
            edge_edge_closest_point(ea0, ea1, eb0, eb1);
        friction_constraint_set.ee_constraints.back().tangent_basis =
            edge_edge_tangent_basis(ea0, ea1, eb0, eb1);
        friction_constraint_set.ee_constraints.back().normal_force_magnitude =
            compute_normal_force_magnitude(
                // The distance type is known because mollified PP and PE were
                // skipped above.
                edge_edge_distance(
                    ea0, ea1, eb0, eb1, EdgeEdgeDistanceType::EA_EB),
                dhat, barrier_stiffness, ee_constraint.minimum_distance);
        friction_constraint_set.ee_constraints.back().mu = mu;
    }

    friction_constraint_set.fv_constraints.reserve(
        contact_constraint_set.fv_constraints.size());
    for (const auto& fv_constraint : contact_constraint_set.fv_constraints) {
        const auto& p = V.row(fv_constraint.vertex_index);
        const auto& t0 = V.row(F(fv_constraint.face_index, 0));
        const auto& t1 = V.row(F(fv_constraint.face_index, 1));
        const auto& t2 = V.row(F(fv_constraint.face_index, 2));

        friction_constraint_set.fv_constraints.emplace_back(fv_constraint);

        friction_constraint_set.fv_constraints.back().closest_point =
            point_triangle_closest_point(p, t0, t1, t2);
        friction_constraint_set.fv_constraints.back().tangent_basis =
            point_triangle_tangent_basis(p, t0, t1, t2);
        friction_constraint_set.fv_constraints.back().normal_force_magnitude =
            compute_normal_force_magnitude(
                point_triangle_distance(
                    p, t0, t1, t2, PointTriangleDistanceType::P_T),
                dhat, barrier_stiffness, fv_constraint.minimum_distance);
        friction_constraint_set.fv_constraints.back().mu = mu;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// compute_friction_potential() in friction.tpp
//
///////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd compute_friction_potential_gradient(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h)
{
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
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h,
    bool project_hessian_to_psd)
{
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

} // namespace ipc
