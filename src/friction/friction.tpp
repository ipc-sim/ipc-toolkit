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

template <typename DerivedRelUi, typename T = typename DerivedRelUi::Scalar>
inline T compute_friction_potential_common(
    const FrictionConstraint& constraint,
    const Eigen::MatrixBase<DerivedRelUi>& rel_ui,
    double epsv_times_h)
{
    return constraint.mu * constraint.normal_force_magnitude
        * f0_SF(
               (rel_ui.transpose() * constraint.tangent_basis.cast<T>())
                   .squaredNorm(),
               epsv_times_h);
}

template <typename DerivedDP0, typename DerivedDP1, typename T>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1,
    const VertexVertexFrictionConstraint& vv_constraint,
    double epsv_times_h)
{
    Eigen::Vector3<T> rel_disp = point_point_relative_displacement(dp0, dp1);
    return vv_constraint.multiplicity
        * compute_friction_potential_common(
               vv_constraint, rel_disp, epsv_times_h);
}

template <
    typename DerivedDP,
    typename DerivedDE0,
    typename DerivedDE1,
    typename T>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const EdgeVertexFrictionConstraint& ev_constraint,
    double epsv_times_h)
{
    Eigen::Vector3<T> rel_disp = point_edge_relative_displacement(
        dp, de0, de1, T(ev_constraint.closest_point[0]));
    return ev_constraint.multiplicity
        * compute_friction_potential_common(
               ev_constraint, rel_disp, epsv_times_h);
}

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1,
    typename T>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const EdgeEdgeFrictionConstraint& ee_constraint,
    double epsv_times_h)
{
    Eigen::Vector3<T> rel_disp = edge_edge_relative_displacement(
        dea0, dea1, deb0, deb1, ee_constraint.closest_point.cast<T>());
    return compute_friction_potential_common(
        ee_constraint, rel_disp, epsv_times_h);
}

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2,
    typename T>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const FaceVertexFrictionConstraint& fv_constraint,
    double epsv_times_h)
{
    Eigen::Vector3<T> rel_disp = point_triangle_relative_displacement(
        dp, dt0, dt1, dt2, fv_constraint.closest_point.cast<T>());
    return compute_friction_potential_common(
        fv_constraint, rel_disp, epsv_times_h);
}

template <typename T>
T compute_friction_potential(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixX<T>& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h)
{
    // absolute linear dislacement of each point
    auto U = V1 - V0.cast<T>();

    T friction_potential(0);

    // TODO: 2D

    for (const auto& vv_constraint : friction_constraint_set.vv_constraints) {
        friction_potential += compute_friction_potential(
            U.row(vv_constraint.vertex0_index),
            U.row(vv_constraint.vertex1_index), vv_constraint, epsv_times_h);
    }

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        friction_potential += compute_friction_potential(
            U.row(ev_constraint.vertex_index),
            U.row(E(ev_constraint.edge_index, 0)),
            U.row(E(ev_constraint.edge_index, 1)), ev_constraint, epsv_times_h);
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        friction_potential += compute_friction_potential(
            U.row(E(ee_constraint.edge0_index, 0)),
            U.row(E(ee_constraint.edge0_index, 1)),
            U.row(E(ee_constraint.edge1_index, 0)),
            U.row(E(ee_constraint.edge1_index, 1)), ee_constraint,
            epsv_times_h);
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        friction_potential += compute_friction_potential(
            U.row(fv_constraint.vertex_index),
            U.row(F(fv_constraint.face_index, 0)),
            U.row(F(fv_constraint.face_index, 1)),
            U.row(F(fv_constraint.face_index, 2)), fv_constraint, epsv_times_h);
    }

    return friction_potential;
}

///////////////////////////////////////////////////////////////////////////////

Eigen::Vector2d compute_friction_potential_gradient_common(
    const FrictionConstraint& constraint,
    const Eigen::Vector3d& relative_displacement,
    double epsv_times_h);

template <typename DerivedDP0, typename DerivedDP1>
inline Eigen::Matrix<double, 2, 3> compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1,
    const VertexVertexFrictionConstraint& vv_constraint,
    double epsv_times_h)
{
    Eigen::Vector2d tangent_relative_displacement =
        compute_friction_potential_gradient_common(
            vv_constraint, point_point_relative_displacement(dp0, dp1),
            epsv_times_h);
    tangent_relative_displacement *= vv_constraint.multiplicity;

    return point_point_relative_mesh_displacement(
        tangent_relative_displacement, vv_constraint.tangent_basis);
}

template <typename DerivedDP, typename DerivedDE0, typename DerivedDE1>
inline Eigen::Matrix3d compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const EdgeVertexFrictionConstraint& ev_constraint,
    double epsv_times_h)
{
    Eigen::Vector3d rel_disp = point_edge_relative_displacement(
        dp, de0, de1, ev_constraint.closest_point[0]);
    Eigen::Vector2d tangent_relative_displacement =
        compute_friction_potential_gradient_common(
            ev_constraint, rel_disp, epsv_times_h);
    tangent_relative_displacement *= ev_constraint.multiplicity;

    return point_edge_relative_mesh_displacement(
        tangent_relative_displacement, ev_constraint.tangent_basis,
        ev_constraint.closest_point[0]);
}

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1>
inline Eigen::Matrix<double, 4, 3> compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const EdgeEdgeFrictionConstraint& ee_constraint,
    double epsv_times_h)
{
    Eigen::Vector3d rel_disp = edge_edge_relative_displacement(
        dea0, dea1, deb0, deb1, ee_constraint.closest_point);
    Eigen::Vector2d tangent_relative_displacement =
        compute_friction_potential_gradient_common(
            ee_constraint, rel_disp, epsv_times_h);

    return edge_edge_relative_mesh_displacements(
        tangent_relative_displacement, ee_constraint.tangent_basis,
        ee_constraint.closest_point);
}

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2>
inline Eigen::Matrix<double, 4, 3> compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const FaceVertexFrictionConstraint& fv_constraint,
    double epsv_times_h)
{
    Eigen::Vector3d rel_disp = point_triangle_relative_displacement(
        dp, dt0, dt1, dt2, fv_constraint.closest_point);
    Eigen::Vector2d tangent_relative_displacement =
        compute_friction_potential_gradient_common(
            fv_constraint, rel_disp, epsv_times_h);

    return point_triangle_relative_mesh_displacements(
        tangent_relative_displacement, fv_constraint.tangent_basis,
        fv_constraint.closest_point);
}

///////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXd compute_friction_potential_hessian_common(
    const FrictionConstraint& constraint,
    const Eigen::Vector3d& relative_displacement,
    const Eigen::MatrixXd& TT,
    const double epsv_times_h,
    bool project_to_psd,
    const int multiplicity = 1);

template <typename DerivedDP0, typename DerivedDP1>
inline Eigen::MatrixXd compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1,
    const VertexVertexFrictionConstraint& vv_constraint,
    double epsv_times_h,
    bool project_to_psd)
{
    Eigen::Vector3d relative_displacement =
        point_point_relative_displacement(dp0, dp1);

    Eigen::Matrix<double, 2, 12> TT;
    point_point_TT(vv_constraint.tangent_basis, TT);

    return compute_friction_potential_hessian_common(
        vv_constraint, relative_displacement, TT, epsv_times_h, project_to_psd,
        vv_constraint.multiplicity);
}

template <typename DerivedDP, typename DerivedDE0, typename DerivedDE1>
inline Eigen::MatrixXd compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const EdgeVertexFrictionConstraint& ev_constraint,
    double epsv_times_h,
    bool project_to_psd)
{
    Eigen::Vector3d relative_displacement = point_edge_relative_displacement(
        dp, de0, de1, ev_constraint.closest_point[0]);

    Eigen::Matrix<double, 2, 12> TT;
    point_edge_TT(
        ev_constraint.tangent_basis, ev_constraint.closest_point[0], TT);

    return compute_friction_potential_hessian_common(
        ev_constraint, relative_displacement, TT, epsv_times_h, project_to_psd,
        ev_constraint.multiplicity);
}

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1>
inline Eigen::MatrixXd compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const EdgeEdgeFrictionConstraint& ee_constraint,
    double epsv_times_h,
    bool project_to_psd)
{
    Eigen::Vector3d relative_displacement = edge_edge_relative_displacement(
        dea0, dea1, deb0, deb1, ee_constraint.closest_point);

    Eigen::Matrix<double, 2, 12> TT;
    edge_edge_TT(ee_constraint.tangent_basis, ee_constraint.closest_point, TT);

    return compute_friction_potential_hessian_common(
        ee_constraint, relative_displacement, TT, epsv_times_h, project_to_psd);
}

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2>
inline Eigen::MatrixXd compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const FaceVertexFrictionConstraint& fv_constraint,
    double epsv_times_h,
    bool project_to_psd)
{
    Eigen::Vector3d relative_displacement =
        point_triangle_relative_displacement(
            dp, dt0, dt1, dt2, fv_constraint.closest_point);

    Eigen::Matrix<double, 2, 12> TT;
    point_triangle_TT(
        fv_constraint.tangent_basis, fv_constraint.closest_point, TT);

    return compute_friction_potential_hessian_common(
        fv_constraint, relative_displacement, TT, epsv_times_h, project_to_psd);
}

} // namespace ipc
