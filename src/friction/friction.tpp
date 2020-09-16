#include "friction.hpp"

#include <Eigen/Sparse>

#include <barrier/barrier.hpp>
#include <distance/edge_edge.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_triangle.hpp>
#include <friction/closest_point.hpp>
#include <friction/relative_displacement.hpp>
#include <friction/tangent_basis.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/local_hessian_to_global_triplets.hpp>

namespace ipc {

template <typename T>
T compute_friction_potential(
    const Eigen::MatrixXd& V0,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    const Eigen::VectorXd& normal_force_magnitudes,
    double epsv_times_h,
    double mu)
{
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
    typedef Eigen::Matrix<T, 3, 1> Vector3T;

    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    // absolute linear dislacement of each point
    MatrixXT U = V1 - V0.cast<T>();

    T friction_potential = 0;

    int constraint_i = 0;

    auto constraint_friction_potential = [&](const Vector3T& rel_ui) {
        const int& ci = constraint_i;
        return normal_force_magnitudes[ci]
            * f0_SF(
                   (rel_ui.transpose() * tangent_bases[ci].cast<T>())
                       .squaredNorm(),
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

} // namespace ipc
