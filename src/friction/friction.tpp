#include "friction.hpp"

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

template <typename T>
T compute_friction_potential(
    const Eigen::MatrixXd& V0,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h)
{
    typedef Eigen::Matrix<T, 3, 1> Vector3T;

    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    // absolute linear dislacement of each point
    auto U = V1 - V0.cast<T>();

    T friction_potential(0);

    auto compute_common = [&](const FrictionConstraint& constraint,
                              const Vector3T& rel_ui) {
        return constraint.mu * constraint.normal_force_magnitude
            * f0_SF(
                   (rel_ui.transpose() * constraint.tangent_basis.cast<T>())
                       .squaredNorm(),
                   epsv_times_h);
    };

    // TODO: 2D

    for (const auto& vv_constraint : friction_constraint_set.vv_constraints) {
        const auto& dp0 = U.row(vv_constraint.vertex0_index);
        const auto& dp1 = U.row(vv_constraint.vertex1_index);

        auto rel_disp = point_point_relative_displacement(dp0, dp1);
        friction_potential += vv_constraint.multiplicity
            * compute_common(vv_constraint, rel_disp);
    }

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        const auto& dp = U.row(ev_constraint.vertex_index);
        const auto& de0 = U.row(E(ev_constraint.edge_index, 0));
        const auto& de1 = U.row(E(ev_constraint.edge_index, 1));

        auto rel_disp = point_edge_relative_displacement(
            dp, de0, de1, T(ev_constraint.closest_point[0]));
        friction_potential += ev_constraint.multiplicity
            * compute_common(ev_constraint, rel_disp);
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        const auto& dea0 = U.row(E(ee_constraint.edge0_index, 0));
        const auto& dea1 = U.row(E(ee_constraint.edge0_index, 1));
        const auto& deb0 = U.row(E(ee_constraint.edge1_index, 0));
        const auto& deb1 = U.row(E(ee_constraint.edge1_index, 1));

        auto rel_disp = edge_edge_relative_displacement(
            dea0, dea1, deb0, deb1, ee_constraint.closest_point.cast<T>());
        friction_potential += compute_common(ee_constraint, rel_disp);
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        const auto& dp = U.row(fv_constraint.vertex_index);
        const auto& dt0 = U.row(F(fv_constraint.face_index, 0));
        const auto& dt1 = U.row(F(fv_constraint.face_index, 1));
        const auto& dt2 = U.row(F(fv_constraint.face_index, 2));

        auto rel_disp = point_triangle_relative_displacement(
            dp, dt0, dt1, dt2, fv_constraint.closest_point.cast<T>());
        friction_potential += compute_common(fv_constraint, rel_disp);
    }

    return friction_potential;
}

} // namespace ipc
