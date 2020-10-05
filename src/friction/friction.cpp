#include <ipc/friction/friction.hpp>

#include <Eigen/Sparse>

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
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

void compute_friction_bases(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    Constraints& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    Eigen::VectorXd& normal_force_magnitudes)
{
    double dhat_squared = dhat * dhat;

    // TODO: ignore EE constraints that are close to parallel
    friction_constraint_set = contact_constraint_set;

    closest_points.reserve(friction_constraint_set.size());
    tangent_bases.reserve(friction_constraint_set.size());
    normal_force_magnitudes.resize(friction_constraint_set.size());

    int constraint_i = 0;

    // TODO: Point-point constraints
    // for (const auto& vv_constraint : friction_constraint_set.vv_constraints)
    // {
    //     const auto& p0 = V.row(ev_constraint.vertex0_index);
    //     const auto& p1 = V.row(ev_constraint.vertex1_index);
    //
    //     Eigen::Vector1d alpha_vec;
    //     alpha_vec << -1;
    //     closestPoint.push_back(alpha_vec);
    //
    //     tangent_bases.push_back(point_point_tangent_basis(p0, p1));
    //
    //     normal_force_magnitudes[constraint_i] =
    //     compute_normal_force_magnitude(
    //         point_point_distance(p0, p1), dhat_squared, barrier_stiffness);
    //
    //     constraint_i++;
    // }

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        const auto& p = V.row(ev_constraint.vertex_index);
        const auto& e0 = V.row(E(ev_constraint.edge_index, 0));
        const auto& e1 = V.row(E(ev_constraint.edge_index, 1));

        double alpha = point_edge_closest_point(p, e0, e1);
        Eigen::Vector1d alpha_vec;
        alpha_vec << alpha;
        closest_points.push_back(alpha_vec);

        tangent_bases.push_back(point_edge_tangent_basis(p, e0, e1));

        normal_force_magnitudes[constraint_i] = compute_normal_force_magnitude(
            point_edge_distance(p, e0, e1), dhat_squared, barrier_stiffness);

        constraint_i++;
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        const auto& ea0 = V.row(E(ee_constraint.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_constraint.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_constraint.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_constraint.edge1_index, 1));

        closest_points.push_back(edge_edge_closest_point(ea0, ea1, eb0, eb1));
        tangent_bases.push_back(edge_edge_tangent_basis(ea0, ea1, eb0, eb1));
        normal_force_magnitudes[constraint_i] = compute_normal_force_magnitude(
            edge_edge_distance(ea0, ea1, eb0, eb1), dhat_squared,
            barrier_stiffness);

        constraint_i++;
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        const auto& p = V.row(fv_constraint.vertex_index);
        const auto& t0 = V.row(F(fv_constraint.face_index, 0));
        const auto& t1 = V.row(F(fv_constraint.face_index, 1));
        const auto& t2 = V.row(F(fv_constraint.face_index, 2));

        closest_points.push_back(point_triangle_closest_point(p, t0, t1, t2));
        tangent_bases.push_back(point_triangle_tangent_basis(p, t0, t1, t2));
        normal_force_magnitudes[constraint_i] = compute_normal_force_magnitude(
            point_triangle_distance(p, t0, t1, t2), dhat_squared,
            barrier_stiffness);

        constraint_i++;
    }
}

/*
double compute_friction_potential(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    const Eigen::VectorXd& normal_force_magnitudes,
    double epsv_times_h,
    double mu)
{
    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    Eigen::MatrixXd U = V1 - V0; // absolute linear dislacement of each point

    double friction_potential = 0;

    int constraint_i = 0;

    auto constraint_friction_potential = [&](const Eigen::Vector3d& rel_ui) {
        const int& ci = constraint_i;
        return normal_force_magnitudes[ci]
            * f0_SF(
                   (rel_ui.transpose() * tangent_bases[ci]).squaredNorm(),
                   epsv_times_h);
    };

    // TODO: 2D
    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        const auto& dp = U.row(ev_constraint.vertex_index);
        const auto& de0 = U.row(E(ev_constraint.edge_index, 0));
        const auto& de1 = U.row(E(ev_constraint.edge_index, 1));

        friction_potential +=
            constraint_friction_potential(point_edge_relative_displacement(
                dp, de0, de1, closest_points[constraint_i][0]));

        constraint_i++;
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        const auto& dea0 = U.row(E(ee_constraint.edge0_index, 0));
        const auto& dea1 = U.row(E(ee_constraint.edge0_index, 1));
        const auto& deb0 = U.row(E(ee_constraint.edge1_index, 0));
        const auto& deb1 = U.row(E(ee_constraint.edge1_index, 1));

        friction_potential +=
            constraint_friction_potential(edge_edge_relative_displacement(
                dea0, dea1, deb0, deb1, closest_points[constraint_i]));

        constraint_i++;
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        const auto& dp = U.row(fv_constraint.vertex_index);
        const auto& dt0 = U.row(F(fv_constraint.face_index, 0));
        const auto& dt1 = U.row(F(fv_constraint.face_index, 1));
        const auto& dt2 = U.row(F(fv_constraint.face_index, 2));

        friction_potential +=
            constraint_friction_potential(point_triangle_relative_displacement(
                dp, dt0, dt1, dt2, closest_points[constraint_i]));

        constraint_i++;
    }

    // TODO: μ per constraint
    return mu * friction_potential;
}
*/

Eigen::VectorXd compute_friction_potential_gradient(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    const Eigen::VectorXd& normal_force_magnitudes,
    double epsv_times_h,
    double mu)
{
    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    Eigen::MatrixXd U = V1 - V0; // absolute linear dislacement of each point
    int dim = U.cols();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(U.size());

    int constraint_i = 0;

    auto foo = [&](const Eigen::Vector3d& rel_ui) {
        Eigen::Vector2d tangent_relative_displacement =
            tangent_bases[constraint_i].transpose() * rel_ui;

        double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
            tangent_relative_displacement.squaredNorm(), epsv_times_h);

        tangent_relative_displacement *=
            f1_div_rel_disp_norm * mu * normal_force_magnitudes[constraint_i];

        return tangent_relative_displacement;
    };

    // TODO: 2D
    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        const auto& dp = U.row(ev_constraint.vertex_index);
        const auto& de0 = U.row(E(ev_constraint.edge_index, 0));
        const auto& de1 = U.row(E(ev_constraint.edge_index, 1));

        Eigen::Vector2d tangent_relative_displacement =
            foo(point_edge_relative_displacement(
                dp, de0, de1, closest_points[constraint_i][0]));

        Eigen::Matrix3d mesh_displacements =
            point_edge_relative_mesh_displacement(
                tangent_relative_displacement, tangent_bases[constraint_i],
                closest_points[constraint_i][0]);

        std::vector<long> ids = { { ev_constraint.vertex_index,
                                    E(ev_constraint.edge_index, 0),
                                    E(ev_constraint.edge_index, 1) } };
        for (int i = 0; i < mesh_displacements.rows(); i++) {
            grad.segment(dim * ids[i], dim) += mesh_displacements.row(i);
        }

        constraint_i++;
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        const auto& dea0 = U.row(E(ee_constraint.edge0_index, 0));
        const auto& dea1 = U.row(E(ee_constraint.edge0_index, 1));
        const auto& deb0 = U.row(E(ee_constraint.edge1_index, 0));
        const auto& deb1 = U.row(E(ee_constraint.edge1_index, 1));

        Eigen::Vector2d tangent_relative_displacement =
            foo(edge_edge_relative_displacement(
                dea0, dea1, deb0, deb1, closest_points[constraint_i]));

        Eigen::Matrix<double, 4, 3> mesh_displacements =
            edge_edge_relative_mesh_displacements(
                tangent_relative_displacement, tangent_bases[constraint_i],
                closest_points[constraint_i]);

        std::vector<long> ids = {
            { E(ee_constraint.edge0_index, 0), E(ee_constraint.edge0_index, 1),
              E(ee_constraint.edge1_index, 0), E(ee_constraint.edge1_index, 1) }
        };
        for (int i = 0; i < mesh_displacements.rows(); i++) {
            grad.segment(dim * ids[i], dim) += mesh_displacements.row(i);
        }

        constraint_i++;
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        const auto& dp = U.row(fv_constraint.vertex_index);
        const auto& dt0 = U.row(F(fv_constraint.face_index, 0));
        const auto& dt1 = U.row(F(fv_constraint.face_index, 1));
        const auto& dt2 = U.row(F(fv_constraint.face_index, 2));

        Eigen::Vector2d tangent_relative_displacement =
            foo(point_triangle_relative_displacement(
                dp, dt0, dt1, dt2, closest_points[constraint_i]));

        Eigen::Matrix<double, 4, 3> mesh_displacements =
            point_triangle_relative_mesh_displacements(
                tangent_relative_displacement, tangent_bases[constraint_i],
                closest_points[constraint_i]);

        std::vector<long> ids = {
            { fv_constraint.vertex_index, F(fv_constraint.face_index, 0),
              F(fv_constraint.face_index, 1), F(fv_constraint.face_index, 2) }
        };
        for (int i = 0; i < mesh_displacements.rows(); i++) {
            grad.segment(dim * ids[i], dim) += mesh_displacements.row(i);
        }

        constraint_i++;
    }

    return grad;
}

Eigen::SparseMatrix<double> compute_friction_potential_hessian(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    const Eigen::VectorXd& normal_force_magnitudes,
    double epsv_times_h,
    double mu)
{
    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    Eigen::MatrixXd U = V1 - V0; // absolute linear dislacement of each point
    int dim = U.cols();
    int dim_sq = dim * dim;

    std::vector<Eigen::Triplet<double>> hess_triplets;
    hess_triplets.reserve(
        friction_constraint_set.ev_constraints.size() * /*3*3=*/9 * dim_sq
        + friction_constraint_set.ee_constraints.size() * /*4*4=*/16 * dim_sq
        + friction_constraint_set.fv_constraints.size() * /*4*4=*/16 * dim_sq);

    int constraint_i = 0;

    auto compute_common = [&](const Eigen::Vector3d& relative_displacement,
                              const Eigen::MatrixXd& TT,
                              const std::vector<long>& ids) {
        Eigen::Vector2d tangent_relative_displacement =
            tangent_bases[constraint_i].transpose() * relative_displacement;

        double tangent_relative_displacement_sqnorm =
            tangent_relative_displacement.squaredNorm();

        double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
            tangent_relative_displacement_sqnorm, epsv_times_h);
        double f2_term =
            f2_SF(tangent_relative_displacement_sqnorm, epsv_times_h);

        Eigen::MatrixXd local_hess;

        if (tangent_relative_displacement_sqnorm >= epsv_times_h_squared) {
            // no SPD projection needed
            Eigen::Vector2d ubar(
                -tangent_relative_displacement[1],
                tangent_relative_displacement[0]);
            local_hess = (TT.transpose()
                          * ((mu * normal_force_magnitudes[constraint_i]
                              * f1_div_rel_disp_norm
                              / tangent_relative_displacement_sqnorm)
                             * ubar))
                * (ubar.transpose() * TT);
        } else {
            double tangent_relative_displacement_norm =
                sqrt(tangent_relative_displacement_sqnorm);
            if (tangent_relative_displacement_norm == 0) {
                // no SPD projection needed
                local_hess = ((mu * normal_force_magnitudes[constraint_i]
                               * f1_div_rel_disp_norm)
                              * TT.transpose())
                    * TT;
            } else {
                // only need to project the inner 2x2 matrix to SPD
                Eigen::Matrix2d inner_hess =
                    ((f2_term / tangent_relative_displacement_norm)
                     * tangent_relative_displacement)
                    * tangent_relative_displacement.transpose();
                inner_hess.diagonal().array() += f1_div_rel_disp_norm;
                inner_hess =
                    project_to_psd(inner_hess); // This is not PD it is PSD
                inner_hess *= mu * normal_force_magnitudes[constraint_i];

                // tensor product:
                local_hess = TT.transpose() * inner_hess * TT;
            }
        }

        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    };

    // TODO: 2D

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        const auto& dp = U.row(ev_constraint.vertex_index);
        const auto& de0 = U.row(E(ev_constraint.edge_index, 0));
        const auto& de1 = U.row(E(ev_constraint.edge_index, 1));

        Eigen::Vector3d relative_displacement =
            point_edge_relative_displacement(
                dp, de0, de1, closest_points[constraint_i][0]);

        Eigen::Matrix<double, 2, 12> TT;
        point_edge_TT(
            tangent_bases[constraint_i], closest_points[constraint_i][0], TT);

        std::vector<long> ids = { { ev_constraint.vertex_index,
                                    E(ev_constraint.edge_index, 0),
                                    E(ev_constraint.edge_index, 1) } };
        compute_common(relative_displacement, TT, ids);

        constraint_i++;
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        const auto& dea0 = U.row(E(ee_constraint.edge0_index, 0));
        const auto& dea1 = U.row(E(ee_constraint.edge0_index, 1));
        const auto& deb0 = U.row(E(ee_constraint.edge1_index, 0));
        const auto& deb1 = U.row(E(ee_constraint.edge1_index, 1));

        Eigen::Vector3d relative_displacement = edge_edge_relative_displacement(
            dea0, dea1, deb0, deb1, closest_points[constraint_i]);

        Eigen::Matrix<double, 2, 12> TT;
        edge_edge_TT(
            tangent_bases[constraint_i], closest_points[constraint_i], TT);

        std::vector<long> ids = {
            { E(ee_constraint.edge0_index, 0), E(ee_constraint.edge0_index, 1),
              E(ee_constraint.edge1_index, 0), E(ee_constraint.edge1_index, 1) }
        };
        compute_common(relative_displacement, TT, ids);

        constraint_i++;
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        const auto& dp = U.row(fv_constraint.vertex_index);
        const auto& dt0 = U.row(F(fv_constraint.face_index, 0));
        const auto& dt1 = U.row(F(fv_constraint.face_index, 1));
        const auto& dt2 = U.row(F(fv_constraint.face_index, 2));

        Eigen::Vector3d relative_displacement =
            point_triangle_relative_displacement(
                dp, dt0, dt1, dt2, closest_points[constraint_i]);

        Eigen::Matrix<double, 2, 12> TT;
        point_triangle_TT(
            tangent_bases[constraint_i], closest_points[constraint_i], TT);

        std::vector<long> ids = {
            { fv_constraint.vertex_index, F(fv_constraint.face_index, 0),
              F(fv_constraint.face_index, 1), F(fv_constraint.face_index, 2) }
        };
        compute_common(relative_displacement, TT, ids);

        constraint_i++;
    }

    Eigen::SparseMatrix<double> hess(U.size(), U.size());
    hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
    return hess;
}

} // namespace ipc
