#include "friction.hpp"

#include <barrier/barrier.hpp>
#include <distance/edge_edge.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_triangle.hpp>
#include <friction/closest_point.hpp>
#include <friction/relative_displacement.hpp>
#include <friction/tangent_basis.hpp>
#include <utils/eigen_ext.hpp>

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
    const ccd::Candidates& contact_constraint_set,
    double dhat_squared,
    double barrier_stiffness,
    ccd::Candidates& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    Eigen::VectorXd& normal_force_magnitudes)
{
    // TODO: ignore EE constraints that are close to parallel
    friction_constraint_set = contact_constraint_set;

    closest_points.reserve(friction_constraint_set.size());
    tangent_bases.reserve(friction_constraint_set.size());
    normal_force_magnitudes.resize(friction_constraint_set.size());

    int constraint_i = 0;

    // TODO: Point-point constraints
    // for (const auto& vv_candidate : friction_constraint_set.vv_candidates) {
    //     const auto& p0 = V.row(ev_candidate.vertex0_index);
    //     const auto& p1 = V.row(ev_candidate.vertex1_index);
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

    for (const auto& ev_candidate : friction_constraint_set.ev_candidates) {
        const auto& p = V.row(ev_candidate.vertex_index);
        const auto& e0 = V.row(E(ev_candidate.edge_index, 0));
        const auto& e1 = V.row(E(ev_candidate.edge_index, 1));

        double alpha = point_edge_closest_point(p, e0, e1);
        Eigen::Vector1d alpha_vec;
        alpha_vec << alpha;
        closest_points.push_back(alpha_vec);

        tangent_bases.push_back(point_edge_tangent_basis(p, e0, e1));

        normal_force_magnitudes[constraint_i] = compute_normal_force_magnitude(
            point_edge_distance(p, e0, e1), dhat_squared, barrier_stiffness);

        constraint_i++;
    }

    for (const auto& ee_candidate : friction_constraint_set.ee_candidates) {
        const auto& ea0 = V.row(E(ee_candidate.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_candidate.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_candidate.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_candidate.edge1_index, 1));

        closest_points.push_back(edge_edge_closest_point(ea0, ea1, eb0, eb1));
        tangent_bases.push_back(edge_edge_tangent_basis(ea0, ea1, eb0, eb1));
        normal_force_magnitudes[constraint_i] = compute_normal_force_magnitude(
            edge_edge_distance(ea0, ea1, eb0, eb1), dhat_squared,
            barrier_stiffness);

        constraint_i++;
    }

    for (const auto& fv_candidate : friction_constraint_set.fv_candidates) {
        const auto& p = V.row(fv_candidate.vertex_index);
        const auto& t0 = V.row(F(fv_candidate.face_index, 0));
        const auto& t1 = V.row(F(fv_candidate.face_index, 1));
        const auto& t2 = V.row(F(fv_candidate.face_index, 2));

        closest_points.push_back(edge_edge_closest_point(p, t0, t1, t2));
        tangent_bases.push_back(edge_edge_tangent_basis(p, t0, t1, t2));
        normal_force_magnitudes[constraint_i] = compute_normal_force_magnitude(
            point_triangle_distance(p, t0, t1, t2), dhat_squared,
            barrier_stiffness);

        constraint_i++;
    }
}

double compute_friction_potential(
    const Eigen::MatrixXd& V0, // TODO: What is this
    const Eigen::MatrixXd& V1, // This is the current position
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    const Eigen::VectorXd& normal_force_magnitudes,
    double epsv_times_h_squared,
    double mu)
{
    double epsv_times_h = sqrt(epsv_times_h_squared);

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
    for (const auto& ev_candidate : friction_constraint_set.ev_candidates) {
        const auto& dp = U.row(ev_candidate.vertex_index);
        const auto& de0 = U.row(E(ev_candidate.edge_index, 0));
        const auto& de1 = U.row(E(ev_candidate.edge_index, 1));

        friction_potential +=
            constraint_friction_potential(point_edge_relative_displacement(
                dp, de0, de1, closest_points[constraint_i][0]));

        constraint_i++;
    }

    for (const auto& ee_candidate : friction_constraint_set.ee_candidates) {
        const auto& dea0 = U.row(E(ee_candidate.edge0_index, 0));
        const auto& dea1 = U.row(E(ee_candidate.edge0_index, 1));
        const auto& deb0 = U.row(E(ee_candidate.edge1_index, 0));
        const auto& deb1 = U.row(E(ee_candidate.edge1_index, 1));

        friction_potential +=
            constraint_friction_potential(edge_edge_relative_displacement(
                dea0, dea1, deb0, deb1, closest_points[constraint_i]));

        constraint_i++;
    }

    for (const auto& fv_candidate : friction_constraint_set.fv_candidates) {
        const auto& dp = U.row(fv_candidate.vertex_index);
        const auto& dt0 = U.row(F(fv_candidate.face_index, 0));
        const auto& dt1 = U.row(F(fv_candidate.face_index, 1));
        const auto& dt2 = U.row(F(fv_candidate.face_index, 2));

        friction_potential +=
            constraint_friction_potential(point_triangle_relative_displacement(
                dp, dt0, dt1, dt2, closest_points[constraint_i]));

        constraint_i++;
    }

    // TODO: μ per constraint
    return mu * friction_potential;
}

Eigen::VectorXd compute_friction_potential_gradient(
    const Eigen::MatrixXd& V0, // TODO: What is this
    const Eigen::MatrixXd& V1, // This is the current position
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& friction_constraint_set,
    std::vector<Eigen::VectorXd>& closest_points,
    std::vector<Eigen::MatrixXd>& tangent_bases,
    const Eigen::VectorXd& normal_force_magnitudes,
    double epsv_times_h_squared,
    double mu)
{
    double epsv_times_h = sqrt(epsv_times_h_squared);

    Eigen::MatrixXd U = V1 - V0; // absolute linear dislacement of each point
    int dim = U.cols();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(U.size());

    int constraint_i = 0;

    auto foo = [&](const Eigen::Vector3d& rel_ui) {
        Eigen::Vector2d tangent_relative_displacement =
            tangent_bases[constraint_i] * rel_ui;

        double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
            tangent_relative_displacement.squaredNorm(), epsv_times_h);

        tangent_relative_displacement *=
            f1_div_rel_disp_norm * mu * normal_force_magnitudes[constraint_i];

        return tangent_relative_displacement;
    };

    // TODO: 2D
    for (const auto& ev_candidate : friction_constraint_set.ev_candidates) {
        const auto& dp = U.row(ev_candidate.vertex_index);
        const auto& de0 = U.row(E(ev_candidate.edge_index, 0));
        const auto& de1 = U.row(E(ev_candidate.edge_index, 1));

        Eigen::Vector3d relative_displacement =
            point_edge_relative_displacement(
                dp, de0, de1, closest_points[constraint_i]);

        Eigen::Vector2d tangent_relative_displacement =
            tangent_bases[constraint_i] * relative_displacement;

        double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
            tangent_relative_displacement.squaredNorm(), epsv_times_h);

        tangent_relative_displacement *=
            f1_div_rel_disp_norm * mu * normal_force_magnitudes[constraint_i];

        Eigen::Matrix3d mesh_displacements =
            point_edge_relative_mesh_displacement(
                tangent_relative_displacement, tangent_bases[constraint_i],
                closest_points[constraint_i][0]);

        std::vector<long> ids = { { ev_candidate.vertex_index,
                                    E(ev_candidate.edge_index, 0),
                                    E(ev_candidate.edge_index, 1) } };
        for (int i = 0; i < mesh_displacements.rows(); i++) {
            grad.segment(dim * ids[i], dim) += mesh_displacements.row(i);
        }

        constraint_i++;
    }

    for (const auto& ee_candidate : friction_constraint_set.ee_candidates) {
        const auto& dea0 = U.row(E(ee_candidate.edge0_index, 0));
        const auto& dea1 = U.row(E(ee_candidate.edge0_index, 1));
        const auto& deb0 = U.row(E(ee_candidate.edge1_index, 0));
        const auto& deb1 = U.row(E(ee_candidate.edge1_index, 1));

        Eigen::Vector3d relative_displacement = edge_edge_relative_displacement(
            dea0, dea1, deb0, deb1, closest_points[constraint_i]);

        Eigen::Vector2d tangent_relative_displacement =
            tangent_bases[constraint_i] * relative_displacement;

        double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
            tangent_relative_displacement.squaredNorm(), epsv_times_h);

        tangent_relative_displacement *=
            f1_div_rel_disp_norm * mu * normal_force_magnitudes[constraint_i];

        Eigen::Matrix<double, 4, 3> mesh_displacements =
            edge_edge_relative_mesh_displacements(
                tangent_relative_displacement, tangent_bases[constraint_i],
                closest_points[constraint_i]);

        std::vector<long> ids = {
            { E(ee_candidate.edge0_index, 0), E(ee_candidate.edge0_index, 1),
              E(ee_candidate.edge1_index, 0), E(ee_candidate.edge1_index, 1) }
        };
        for (int i = 0; i < mesh_displacements.rows(); i++) {
            grad.segment(dim * ids[i], dim) += mesh_displacements.row(i);
        }

        constraint_i++;
    }

    for (const auto& fv_candidate : friction_constraint_set.fv_candidates) {
        const auto& dp = U.row(fv_candidate.vertex_index);
        const auto& dt0 = U.row(F(fv_candidate.face_index, 0));
        const auto& dt1 = U.row(F(fv_candidate.face_index, 1));
        const auto& dt2 = U.row(F(fv_candidate.face_index, 2));

        Eigen::Vector3d relative_displacement =
            point_triangle_relative_displacement(
                dp, dt0, dt1, dt2, closest_points[constraint_i]);

        Eigen::Vector2d tangent_relative_displacement =
            tangent_bases[constraint_i] * relative_displacement;

        double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
            tangent_relative_displacement.squaredNorm(), epsv_times_h);

        tangent_relative_displacement *=
            f1_div_rel_disp_norm * mu * normal_force_magnitudes[constraint_i];

        Eigen::Matrix<double, 4, 3> mesh_displacements =
            point_triangle_relative_mesh_displacements(
                tangent_relative_displacement, tangent_bases[constraint_i],
                closest_points[constraint_i]);

        std::vector<long> ids = {
            { fv_candidate.vertex_index, F(fv_candidate.face_index, 0),
              F(fv_candidate.face_index, 1), F(fv_candidate.face_index, 2) }
        };
        for (int i = 0; i < mesh_displacements.rows(); i++) {
            grad.segment(dim * ids[i], dim) += mesh_displacements.row(i);
        }

        constraint_i++;
    }

    return grad;
}

/*
template <class T, int dim = 3>
void Compute_Friction_Gradient(
    MESH_NODE<T, dim>& X,
    MESH_NODE<T, dim>& Xn,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    const std::vector<Eigen::Matrix<T, dim - 1, 1>>& closestPoint,
    const std::vector<Eigen::Matrix<T, dim, dim - 1>>& tanBasis,
    const std::vector<T>& normalForce,
    T epsvh2,
    T mu,
    MESH_NODE_ATTR<T, dim>& nodeAttr)
{
    TIMER_FLAG("Compute_Friction_Gradient");

    T epsvh = std::sqrt(epsvh2);

    MESH_NODE<T, dim> dX(X.size);
    X.deep_copy_to(dX);
    dX.Join(Xn).Par_Each([&](int id, auto data) {
        auto& [dx, xn] = data;
        dx -= xn;
    });

    if constexpr (dim == 3) {
        // TODO: parallelize
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            VECTOR<int, dim + 1> cIVInd =
                constraintSet[cI]; // NOTE: copy to be able to modify in the
                                   // loop if needed
            Eigen::Matrix<T, dim, 1> relDX3D;
            if (cIVInd[0] >= 0) {
                // ++++ edge-edge, no mollified stencil for friction
                const VECTOR<T, 3>& dXea0 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 3>& dXea1 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 3>& dXeb0 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                const VECTOR<T, 3>& dXeb1 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                Eigen::Matrix<T, 3, 1> dea0(dXea0.data), dea1(dXea1.data),
                    deb0(dXeb0.data), deb1(dXeb1.data);

                Edge_Edge_RelDX(
                    dea0, dea1, deb0, deb1, closestPoint[cI][0],
                    closestPoint[cI][1], relDX3D);

                Eigen::Matrix<T, dim - 1, 1> relDX =
                    tanBasis[cI].transpose() * relDX3D;
                T f1_div_relDXNorm;
                f1_SF_Div_RelDXNorm(
                    relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                relDX *= f1_div_relDXNorm * mu * normalForce[cI];

                Eigen::Matrix<T, 12, 1> TTTDX;
                Edge_Edge_RelDXTan_To_Mesh(
                    relDX, tanBasis[cI], closestPoint[cI][0],
                    closestPoint[cI][1], TTTDX);

                for (int i = 0; i < 4; ++i) {
                    VECTOR<T, dim>& g =
                        std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[i]));
                    g += TTTDX.data() + i * dim;
                }
            } else {
                // point-triangle and degenerate edge-edge
                assert(cIVInd[1] >= 0);

                cIVInd[0] = -cIVInd[0] - 1;
                const VECTOR<T, 3>& dXp =
                    std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                Eigen::Matrix<T, 3, 1> dp(dXp.data);
                if (cIVInd[2] < 0) {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& dXp1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> dp1(dXp1.data);

                    Point_Point_RelDX(dp, dp1, relDX3D);

                    Eigen::Matrix<T, dim - 1, 1> relDX =
                        tanBasis[cI].transpose() * relDX3D;
                    T f1_div_relDXNorm;
                    f1_SF_Div_RelDXNorm(
                        relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                    relDX *=
                        f1_div_relDXNorm * -cIVInd[3] * mu * normalForce[cI];

                    Eigen::Matrix<T, 6, 1> TTTDX;
                    Point_Point_RelDXTan_To_Mesh(relDX, tanBasis[cI], TTTDX);

                    for (int i = 0; i < 2; ++i) {
                        VECTOR<T, dim>& g =
                            std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(
                                nodeAttr.Get_Unchecked(cIVInd[i]));
                        g += TTTDX.data() + i * dim;
                    }
                } else if (cIVInd[3] < 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& dXe0 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXe1 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> de0(dXe0.data), de1(dXe1.data);

                    Point_Edge_RelDX(
                        dp, de0, de1, closestPoint[cI][0], relDX3D);

                    Eigen::Matrix<T, dim - 1, 1> relDX =
                        tanBasis[cI].transpose() * relDX3D;
                    T f1_div_relDXNorm;
                    f1_SF_Div_RelDXNorm(
                        relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                    relDX *=
                        f1_div_relDXNorm * -cIVInd[3] * mu * normalForce[cI];

                    Eigen::Matrix<T, 9, 1> TTTDX;
                    Point_Edge_RelDXTan_To_Mesh(
                        relDX, tanBasis[cI], closestPoint[cI][0], TTTDX);

                    for (int i = 0; i < 3; ++i) {
                        VECTOR<T, dim>& g =
                            std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(
                                nodeAttr.Get_Unchecked(cIVInd[i]));
                        g += TTTDX.data() + i * dim;
                    }
                } else {
                    // -+++ PT
                    const VECTOR<T, 3>& dXt0 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXt1 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& dXt2 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> dt0(dXt0.data), dt1(dXt1.data),
                        dt2(dXt2.data);

                    Point_Triangle_RelDX(
                        dp, dt0, dt1, dt2, closestPoint[cI][0],
                        closestPoint[cI][1], relDX3D);

                    Eigen::Matrix<T, dim - 1, 1> relDX =
                        tanBasis[cI].transpose() * relDX3D;
                    T f1_div_relDXNorm;
                    f1_SF_Div_RelDXNorm(
                        relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm * mu * normalForce[cI];

                    Eigen::Matrix<T, 12, 1> TTTDX;
                    Point_Triangle_RelDXTan_To_Mesh(
                        relDX, tanBasis[cI], closestPoint[cI][0],
                        closestPoint[cI][1], TTTDX);

                    for (int i = 0; i < 4; ++i) {
                        VECTOR<T, dim>& g =
                            std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(
                                nodeAttr.Get_Unchecked(cIVInd[i]));
                        g += TTTDX.data() + i * dim;
                    }
                }
            }
        }
    } else {
        // TODO
    }
}

template <class T, int dim = 3>
void Compute_Friction_Hessian(
    MESH_NODE<T, dim>& X,
    MESH_NODE<T, dim>& Xn,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    const std::vector<Eigen::Matrix<T, dim - 1, 1>>& closestPoint,
    const std::vector<Eigen::Matrix<T, dim, dim - 1>>& tanBasis,
    const std::vector<T>& normalForce,
    T epsvh2,
    T mu,
    std::vector<Eigen::Triplet<T>>& triplets)
{
    TIMER_FLAG("Compute_Friction_Hessian");

    T epsvh = std::sqrt(epsvh2);

    MESH_NODE<T, dim> dX(X.size);
    X.deep_copy_to(dX);
    dX.Join(Xn).Par_Each([&](int id, auto data) {
        auto& [dx, xn] = data;
        dx -= xn;
    });

    if constexpr (dim == 3) {
        BASE_STORAGE<int> threads(constraintSet.size());
        int curStartInd = triplets.size();
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            threads.Append(curStartInd);
            const VECTOR<int, 4>& cIVInd = constraintSet[cI];
            if (cIVInd[0] >= 0 || cIVInd[3] >= 0) {
                // EE or PT, 12x12
                curStartInd += 144;
            } else if (cIVInd[2] >= 0) {
                // PE, 9x9
                curStartInd += 81;
            } else {
                // PP, 6x6
                curStartInd += 36;
            }
        }
        triplets.resize(curStartInd);

        threads.Par_Each([&](int cI, auto data) {
            const auto& [tripletStart] = data;
            VECTOR<int, 4> cIVInd =
                constraintSet[cI]; // NOTE: copy to be able to modify in the
                                   // loop if needed
            assert(cIVInd[1] >= 0);

            Eigen::Matrix<T, dim, 1> relDX3D;
            if (cIVInd[0] >= 0) {
                // ++++ edge-edge, no mollified stencil for friction
                const VECTOR<T, 3>& dXea0 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 3>& dXea1 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 3>& dXeb0 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                const VECTOR<T, 3>& dXeb1 =
                    std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                Eigen::Matrix<T, 3, 1> dea0(dXea0.data), dea1(dXea1.data),
                    deb0(dXeb0.data), deb1(dXeb1.data);

                Edge_Edge_RelDX(
                    dea0, dea1, deb0, deb1, closestPoint[cI][0],
                    closestPoint[cI][1], relDX3D);
                Eigen::Matrix<T, dim - 1, 1> relDX =
                    tanBasis[cI].transpose() * relDX3D;
                T relDXSqNorm = relDX.squaredNorm();
                T relDXNorm = std::sqrt(relDXSqNorm);

                Eigen::Matrix<T, 2, 12> TT;
                Edge_Edge_TT(
                    tanBasis[cI], closestPoint[cI][0], closestPoint[cI][1], TT);

                T f1_div_relDXNorm, f2_term;
                f1_SF_Div_RelDXNorm(relDXSqNorm, epsvh, f1_div_relDXNorm);
                f2_SF_Term(relDXSqNorm, epsvh, f2_term);

                Eigen::Matrix<T, 12, 12> HessianI;
                if (relDXSqNorm >= epsvh2) {
                    // no SPD projection needed
                    Eigen::Matrix<T, 2, 1> ubar(-relDX[1], relDX[0]);
                    HessianI = (TT.transpose()
                                * ((mu * normalForce[cI] * f1_div_relDXNorm
                                    / relDXSqNorm)
                                   * ubar))
                        * (ubar.transpose() * TT);
                } else {
                    if (relDXNorm == 0) {
                        // no SPD projection needed
                        HessianI = ((mu * normalForce[cI] * f1_div_relDXNorm)
                                    * TT.transpose())
                            * TT;
                    } else {
                        // only need to project the inner 2x2 matrix to SPD
                        Eigen::Matrix<T, 2, 2> innerMtr =
                            ((f2_term / relDXNorm) * relDX) * relDX.transpose();
                        innerMtr.diagonal().array() += f1_div_relDXNorm;
                        makePD(innerMtr);
                        innerMtr *= mu * normalForce[cI];

                        // tensor product:
                        HessianI = TT.transpose() * innerMtr * TT;
                    }
                }

                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        for (int idI = 0; idI < 3; ++idI) {
                            for (int jdI = 0; jdI < 3; ++jdI) {
                                triplets
                                    [tripletStart + (i * 3 + idI) * 12 + j * 3
                                     + jdI] =
                                        Eigen::Triplet<T>(
                                            cIVInd[i] * 3 + idI,
                                            cIVInd[j] * 3 + jdI,
                                            HessianI(i * 3 + idI, j * 3 + jdI));
                            }
                        }
                    }
                }
            } else {
                // point-triangle and degenerate edge-edge
                assert(cIVInd[1] >= 0);

                cIVInd[0] = -cIVInd[0] - 1;
                const VECTOR<T, 3>& dXp =
                    std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                Eigen::Matrix<T, 3, 1> dp(dXp.data);
                if (cIVInd[2] < 0) {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& dXp1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> dp1(dXp1.data);

                    Point_Point_RelDX(dp, dp1, relDX3D);
                    Eigen::Matrix<T, dim - 1, 1> relDX =
                        tanBasis[cI].transpose() * relDX3D;
                    T relDXSqNorm = relDX.squaredNorm();
                    T relDXNorm = std::sqrt(relDXSqNorm);

                    Eigen::Matrix<T, 2, 6> TT;
                    Point_Point_TT(tanBasis[cI], TT);

                    T f1_div_relDXNorm, f2_term;
                    f1_SF_Div_RelDXNorm(relDXSqNorm, epsvh, f1_div_relDXNorm);
                    f2_SF_Term(relDXSqNorm, epsvh, f2_term);

                    Eigen::Matrix<T, 6, 6> HessianI;
                    if (relDXSqNorm >= epsvh2) {
                        // no SPD projection needed
                        Eigen::Matrix<T, 2, 1> ubar(-relDX[1], relDX[0]);
                        HessianI = (TT.transpose()
                                    * ((-cIVInd[3] * mu * normalForce[cI]
                                        * f1_div_relDXNorm / relDXSqNorm)
                                       * ubar))
                            * (ubar.transpose() * TT);
                    } else {
                        if (relDXNorm == 0) {
                            // no SPD projection needed
                            HessianI = ((-cIVInd[3] * mu * normalForce[cI]
                                         * f1_div_relDXNorm)
                                        * TT.transpose())
                                * TT;
                        } else {
                            // only need to project the inner 2x2 matrix to SPD
                            Eigen::Matrix<T, 2, 2> innerMtr =
                                ((f2_term / relDXNorm) * relDX)
                                * relDX.transpose();
                            innerMtr.diagonal().array() += f1_div_relDXNorm;
                            makePD(innerMtr);
                            innerMtr *= -cIVInd[3] * mu * normalForce[cI];

                            // tensor product:
                            HessianI = TT.transpose() * innerMtr * TT;
                        }
                    }

                    for (int i = 0; i < 2; ++i) {
                        for (int j = 0; j < 2; ++j) {
                            for (int idI = 0; idI < 3; ++idI) {
                                for (int jdI = 0; jdI < 3; ++jdI) {
                                    triplets
                                        [tripletStart + (i * 3 + idI) * 6
                                         + j * 3 + jdI] =
                                            Eigen::Triplet<T>(
                                                cIVInd[i] * 3 + idI,
                                                cIVInd[j] * 3 + jdI,
                                                HessianI(
                                                    i * 3 + idI, j * 3 + jdI));
                                }
                            }
                        }
                    }
                } else if (cIVInd[3] < 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& dXe0 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXe1 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> de0(dXe0.data), de1(dXe1.data);

                    Point_Edge_RelDX(
                        dp, de0, de1, closestPoint[cI][0], relDX3D);
                    Eigen::Matrix<T, dim - 1, 1> relDX =
                        tanBasis[cI].transpose() * relDX3D;
                    T relDXSqNorm = relDX.squaredNorm();
                    T relDXNorm = std::sqrt(relDXSqNorm);

                    Eigen::Matrix<T, 2, 9> TT;
                    Point_Edge_TT(tanBasis[cI], closestPoint[cI][0], TT);

                    T f1_div_relDXNorm, f2_term;
                    f1_SF_Div_RelDXNorm(relDXSqNorm, epsvh, f1_div_relDXNorm);
                    f2_SF_Term(relDXSqNorm, epsvh, f2_term);

                    Eigen::Matrix<T, 9, 9> HessianI;
                    if (relDXSqNorm >= epsvh2) {
                        // no SPD projection needed
                        Eigen::Matrix<T, 2, 1> ubar(-relDX[1], relDX[0]);
                        HessianI = (TT.transpose()
                                    * ((-cIVInd[3] * mu * normalForce[cI]
                                        * f1_div_relDXNorm / relDXSqNorm)
                                       * ubar))
                            * (ubar.transpose() * TT);
                    } else {
                        if (relDXNorm == 0) {
                            // no SPD projection needed
                            HessianI = ((-cIVInd[3] * mu * normalForce[cI]
                                         * f1_div_relDXNorm)
                                        * TT.transpose())
                                * TT;
                        } else {
                            // only need to project the inner 2x2 matrix to SPD
                            Eigen::Matrix<T, 2, 2> innerMtr =
                                ((f2_term / relDXNorm) * relDX)
                                * relDX.transpose();
                            innerMtr.diagonal().array() += f1_div_relDXNorm;
                            makePD(innerMtr);
                            innerMtr *= -cIVInd[3] * mu * normalForce[cI];

                            // tensor product:
                            HessianI = TT.transpose() * innerMtr * TT;
                        }
                    }

                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            for (int idI = 0; idI < 3; ++idI) {
                                for (int jdI = 0; jdI < 3; ++jdI) {
                                    triplets
                                        [tripletStart + (i * 3 + idI) * 9
                                         + j * 3 + jdI] =
                                            Eigen::Triplet<T>(
                                                cIVInd[i] * 3 + idI,
                                                cIVInd[j] * 3 + jdI,
                                                HessianI(
                                                    i * 3 + idI, j * 3 + jdI));
                                }
                            }
                        }
                    }
                } else {
                    // -+++ PT
                    const VECTOR<T, 3>& dXt0 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXt1 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& dXt2 =
                        std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> dt0(dXt0.data), dt1(dXt1.data),
                        dt2(dXt2.data);

                    Point_Triangle_RelDX(
                        dp, dt0, dt1, dt2, closestPoint[cI][0],
                        closestPoint[cI][1], relDX3D);
                    Eigen::Matrix<T, dim - 1, 1> relDX =
                        tanBasis[cI].transpose() * relDX3D;
                    T relDXSqNorm = relDX.squaredNorm();
                    T relDXNorm = std::sqrt(relDXSqNorm);

                    Eigen::Matrix<T, 2, 12> TT;
                    Point_Triangle_TT(
                        tanBasis[cI], closestPoint[cI][0], closestPoint[cI][1],
                        TT);

                    T f1_div_relDXNorm, f2_term;
                    f1_SF_Div_RelDXNorm(relDXSqNorm, epsvh, f1_div_relDXNorm);
                    f2_SF_Term(relDXSqNorm, epsvh, f2_term);

                    Eigen::Matrix<T, 12, 12> HessianI;
                    if (relDXSqNorm >= epsvh2) {
                        // no SPD projection needed
                        Eigen::Matrix<T, 2, 1> ubar(-relDX[1], relDX[0]);
                        HessianI = (TT.transpose()
                                    * ((mu * normalForce[cI] * f1_div_relDXNorm
                                        / relDXSqNorm)
                                       * ubar))
                            * (ubar.transpose() * TT);
                    } else {
                        if (relDXNorm == 0) {
                            // no SPD projection needed
                            HessianI =
                                ((mu * normalForce[cI] * f1_div_relDXNorm)
                                 * TT.transpose())
                                * TT;
                        } else {
                            // only need to project the inner 2x2 matrix to SPD
                            Eigen::Matrix<T, 2, 2> innerMtr =
                                ((f2_term / relDXNorm) * relDX)
                                * relDX.transpose();
                            innerMtr.diagonal().array() += f1_div_relDXNorm;
                            makePD(innerMtr);
                            innerMtr *= mu * normalForce[cI];

                            // tensor product:
                            HessianI = TT.transpose() * innerMtr * TT;
                        }
                    }

                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            for (int idI = 0; idI < 3; ++idI) {
                                for (int jdI = 0; jdI < 3; ++jdI) {
                                    triplets
                                        [tripletStart + (i * 3 + idI) * 12
                                         + j * 3 + jdI] =
                                            Eigen::Triplet<T>(
                                                cIVInd[i] * 3 + idI,
                                                cIVInd[j] * 3 + jdI,
                                                HessianI(
                                                    i * 3 + idI, j * 3 + jdI));
                                }
                            }
                        }
                    }
                }
            }
        });
    } else {
        // TODO
    }
}
*/
} // namespace ipc
