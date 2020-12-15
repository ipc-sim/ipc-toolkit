#include <ipc/friction/friction.hpp>

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
#include <ipc/utils/local_hessian_to_global_triplets.hpp>

namespace ipc {

double compute_normal_force_magnitude(
    double distance_squared, double dhat_squared, double barrier_stiffness)
{
    double grad_b = barrier_gradient(distance_squared, dhat_squared);
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
    double dhat_squared = dhat * dhat;

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
                point_point_distance(p0, p1), dhat_squared, barrier_stiffness);
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
                dhat_squared, barrier_stiffness);
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
                dhat_squared, barrier_stiffness);
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
                dhat_squared, barrier_stiffness);
        friction_constraint_set.fv_constraints.back().mu = mu;
    }
}

///////////////////////////////////////////////////////////////////////////////

template <typename DerivedLocalGrad>
void add_local_grad(
    const Eigen::MatrixBase<DerivedLocalGrad>& local_grad,
    const std::vector<long>& ids,
    Eigen::VectorXd& grad)
{
    assert(ids.size() == local_grad.rows());
    int dim = local_grad.cols();
    for (int i = 0; i < local_grad.rows(); i++) {
        grad.segment(dim * ids[i], dim) += local_grad.row(i);
    }
}

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

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(U.size());

    // TODO: 2D

    for (const auto& vv_constraint : friction_constraint_set.vv_constraints) {
        add_local_grad(
            /*local_grad=*/
            compute_friction_potential_gradient(
                U.row(vv_constraint.vertex0_index),
                U.row(vv_constraint.vertex1_index), //
                vv_constraint, epsv_times_h),
            { { vv_constraint.vertex0_index, vv_constraint.vertex1_index } },
            grad);
    }

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        add_local_grad(
            /*local_grad=*/
            compute_friction_potential_gradient(
                U.row(ev_constraint.vertex_index),
                U.row(E(ev_constraint.edge_index, 0)),
                U.row(E(ev_constraint.edge_index, 1)), //
                ev_constraint, epsv_times_h),
            { { ev_constraint.vertex_index, E(ev_constraint.edge_index, 0),
                E(ev_constraint.edge_index, 1) } },
            grad);
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        add_local_grad(
            /*local_grad=*/
            compute_friction_potential_gradient(
                U.row(E(ee_constraint.edge0_index, 0)),
                U.row(E(ee_constraint.edge0_index, 1)),
                U.row(E(ee_constraint.edge1_index, 0)),
                U.row(E(ee_constraint.edge1_index, 1)), //
                ee_constraint, epsv_times_h),
            { { E(ee_constraint.edge0_index, 0),
                E(ee_constraint.edge0_index, 1),
                E(ee_constraint.edge1_index, 0),
                E(ee_constraint.edge1_index, 1) } },
            grad);
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        add_local_grad(
            /*local_grad=*/
            compute_friction_potential_gradient(
                U.row(fv_constraint.vertex_index),
                U.row(F(fv_constraint.face_index, 0)),
                U.row(F(fv_constraint.face_index, 1)),
                U.row(F(fv_constraint.face_index, 2)), //
                fv_constraint, epsv_times_h),
            { { fv_constraint.vertex_index,     //
                F(fv_constraint.face_index, 0), //
                F(fv_constraint.face_index, 1), //
                F(fv_constraint.face_index, 2) } },
            grad);
    }

    return grad;
}

Eigen::Vector2d compute_friction_potential_gradient_common(
    const FrictionConstraint& constraint,
    const Eigen::Vector3d& relative_displacement,
    double epsv_times_h)
{
    Eigen::Vector2d tangent_relative_displacement =
        constraint.tangent_basis.transpose() * relative_displacement;

    double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
        tangent_relative_displacement.squaredNorm(), epsv_times_h);

    tangent_relative_displacement *= f1_div_rel_disp_norm * constraint.mu
        * constraint.normal_force_magnitude;

    return tangent_relative_displacement;
}

///////////////////////////////////////////////////////////////////////////////

Eigen::SparseMatrix<double> compute_friction_potential_hessian(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h,
    bool project_to_psd)
{
    auto U = V1 - V0; // absolute linear dislacement of each point
    int dim = U.cols();
    int dim_sq = dim * dim;

    std::vector<Eigen::Triplet<double>> hess_triplets;
    hess_triplets.reserve(
        friction_constraint_set.ev_constraints.size() * /*3*3=*/9 * dim_sq
        + friction_constraint_set.ee_constraints.size() * /*4*4=*/16 * dim_sq
        + friction_constraint_set.fv_constraints.size() * /*4*4=*/16 * dim_sq);

    // TODO: 2D

    for (const auto& vv_constraint : friction_constraint_set.vv_constraints) {
        local_hessian_to_global_triplets(
            compute_friction_potential_hessian(
                U.row(vv_constraint.vertex0_index),
                U.row(vv_constraint.vertex1_index), //
                vv_constraint, epsv_times_h),
            { { vv_constraint.vertex0_index, vv_constraint.vertex1_index } },
            dim, hess_triplets);
    }

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        local_hessian_to_global_triplets(
            compute_friction_potential_hessian(
                U.row(ev_constraint.vertex_index),
                U.row(E(ev_constraint.edge_index, 0)),
                U.row(E(ev_constraint.edge_index, 1)), //
                ev_constraint, epsv_times_h),
            { { ev_constraint.vertex_index, //
                E(ev_constraint.edge_index, 0),
                E(ev_constraint.edge_index, 1) } },
            dim, hess_triplets);
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        local_hessian_to_global_triplets(
            compute_friction_potential_hessian(
                U.row(E(ee_constraint.edge0_index, 0)),
                U.row(E(ee_constraint.edge0_index, 1)),
                U.row(E(ee_constraint.edge1_index, 0)),
                U.row(E(ee_constraint.edge1_index, 1)), //
                ee_constraint, epsv_times_h),
            { { E(ee_constraint.edge0_index, 0),
                E(ee_constraint.edge0_index, 1),
                E(ee_constraint.edge1_index, 0),
                E(ee_constraint.edge1_index, 1) } },
            dim, hess_triplets);
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        local_hessian_to_global_triplets(
            compute_friction_potential_hessian(
                U.row(fv_constraint.vertex_index),
                U.row(F(fv_constraint.face_index, 0)),
                U.row(F(fv_constraint.face_index, 1)),
                U.row(F(fv_constraint.face_index, 2)), //
                fv_constraint, epsv_times_h),
            { { fv_constraint.vertex_index, F(fv_constraint.face_index, 0),
                F(fv_constraint.face_index, 1),
                F(fv_constraint.face_index, 2) } },
            dim, hess_triplets);
    }

    Eigen::SparseMatrix<double> hess(U.size(), U.size());
    hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
    return hess;
}

Eigen::MatrixXd compute_friction_potential_hessian_common(
    const FrictionConstraint& constraint,
    const Eigen::Vector3d& relative_displacement,
    const Eigen::MatrixXd& TT,
    const double epsv_times_h,
    bool project_to_psd,
    const int multiplicity)
{
    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    Eigen::Vector2d tangent_relative_displacement =
        constraint.tangent_basis.transpose() * relative_displacement;

    double tangent_relative_displacement_sqnorm =
        tangent_relative_displacement.squaredNorm();

    double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
        tangent_relative_displacement_sqnorm, epsv_times_h);
    double f2_term = f2_SF(tangent_relative_displacement_sqnorm, epsv_times_h);

    Eigen::MatrixXd local_hess;

    double scale =
        multiplicity * constraint.mu * constraint.normal_force_magnitude;
    if (tangent_relative_displacement_sqnorm >= epsv_times_h_squared) {
        // no SPD projection needed
        Eigen::Vector2d ubar(
            -tangent_relative_displacement[1],
            tangent_relative_displacement[0]);
        local_hess = (TT.transpose()
                      * ((scale * f1_div_rel_disp_norm
                          / tangent_relative_displacement_sqnorm)
                         * ubar))
            * (ubar.transpose() * TT);
    } else {
        double tangent_relative_displacement_norm =
            sqrt(tangent_relative_displacement_sqnorm);
        if (tangent_relative_displacement_norm == 0) {
            // no SPD projection needed
            local_hess = ((scale * f1_div_rel_disp_norm) * TT.transpose()) * TT;
        } else {
            // only need to project the inner 2x2 matrix to SPD
            Eigen::Matrix2d inner_hess =
                ((f2_term / tangent_relative_displacement_norm)
                 * tangent_relative_displacement)
                * tangent_relative_displacement.transpose();
            inner_hess.diagonal().array() += f1_div_rel_disp_norm;
            if (project_to_psd) {
                inner_hess = Eigen::project_to_psd(inner_hess);
            }
            inner_hess *= scale;

            // tensor product:
            local_hess = TT.transpose() * inner_hess * TT;
        }
    }

    return local_hess;
}

} // namespace ipc
