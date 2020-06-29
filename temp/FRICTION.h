#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <FEM/FRICTION_UTILS.h>
#include <Math/Distance/POINT_POINT.h>
#include <Math/Distance/POINT_EDGE.h>
#include <Math/Distance/POINT_TRIANGLE.h>
#include <Math/Distance/EDGE_EDGE.h>
#include <Math/BARRIER.h>

namespace py = pybind11;
namespace JGSL {

template <class T, int dim = 3>
void Compute_Friction_Basis(MESH_NODE<T, dim>& X,
    const std::vector<VECTOR<int, dim + 1>>& contactConstraintSet,
    std::vector<VECTOR<int, dim + 1>>& constraintSet,
    std::vector<Eigen::Matrix<T, dim - 1, 1>>& closestPoint,
    std::vector<Eigen::Matrix<T, dim, dim - 1>>& tanBasis,
    std::vector<T>& normalForce, 
    T dHat2, T kappa[])
{
    TIMER_FLAG("Compute_Friction_Basis");

    constraintSet.resize(0);
    constraintSet.reserve(contactConstraintSet.size());
    for (const auto& cIVInd : contactConstraintSet) {
        if (!(cIVInd[0] >= 0 && (cIVInd[2] < 0 || cIVInd[3] < 0))) {
            constraintSet.emplace_back(cIVInd);
        }
    }

    closestPoint.resize(constraintSet.size());
    tanBasis.resize(constraintSet.size());
    normalForce.resize(constraintSet.size());
    if constexpr (dim == 3) {
        //TODO: parallelize
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            const VECTOR<int, dim + 1>& cIVInd = constraintSet[cI];
            T dist2;
            if (cIVInd[0] >= 0) {
                // ++++ edge-edge, no mollified stencil for friction
                const VECTOR<T, 3>& Xea0 = std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 3>& Xea1 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 3>& Xeb0 = std::get<0>(X.Get_Unchecked(cIVInd[2]));
                const VECTOR<T, 3>& Xeb1 = std::get<0>(X.Get_Unchecked(cIVInd[3]));
                Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data), eb0(Xeb0.data), eb1(Xeb1.data);

                Edge_Edge_Closest_Point(ea0, ea1, eb0, eb1, closestPoint[cI]);
                Edge_Edge_Tangent_Basis(ea0, ea1, eb0, eb1, tanBasis[cI]);

                Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(cIVInd[1] >= 0);

                const VECTOR<T, 3>& Xp = std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                Eigen::Matrix<T, 3, 1> p(Xp.data);
                if (cIVInd[2] < 0) {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp1 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> p1(Xp1.data);

                    Point_Point_Tangent_Basis(p, p1, tanBasis[cI]);

                    Point_Point_Distance(p, p1, dist2);
                }
                else if (cIVInd[3] < 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& Xe0 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xe1 = std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> e0(Xe0.data), e1(Xe1.data);

                    Point_Edge_Closest_Point(p, e0, e1, closestPoint[cI][0]);
                    Point_Edge_Tangent_Basis(p, e0, e1, tanBasis[cI]);

                    Point_Edge_Distance(p, e0, e1, dist2);
                }
                else {
                    // -+++ PT 
                    const VECTOR<T, 3>& Xt0 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xt1 = std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xt2 = std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> t0(Xt0.data), t1(Xt1.data), t2(Xt2.data);

                    Point_Triangle_Closest_Point(p, t0, t1, t2, closestPoint[cI]);
                    Point_Triangle_Tangent_Basis(p, t0, t1, t2, tanBasis[cI]);

                    Point_Triangle_Distance(p, t0, t1, t2, dist2);
                }
            }
            T bGrad;
            Barrier_Gradient(dist2, dHat2, kappa, bGrad);
            normalForce[cI] = -bGrad * 2 * std::sqrt(dist2);
            // / (h * h) eliminated here
            // multiplicity enforced when computing friction force not here
        }
    }
    else {
        //TODO
    }
}

template <class T, int dim = 3>
void Compute_Friction_Potential(
    MESH_NODE<T, dim>& X, MESH_NODE<T, dim>& Xn,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    const std::vector<Eigen::Matrix<T, dim - 1, 1>>& closestPoint,
    const std::vector<Eigen::Matrix<T, dim, dim - 1>>& tanBasis,
    const std::vector<T>& normalForce,
    T epsvh2, T mu, T& E)
{
    TIMER_FLAG("Compute_Friction_Potential");

    T epsvh = std::sqrt(epsvh2);

    MESH_NODE<T, dim> dX(X.size);
    X.deep_copy_to(dX);
    dX.Join(Xn).Par_Each([&](int id, auto data) {
        auto &[dx, xn] = data;
        dx -= xn;
    });

    std::vector<T> EI(constraintSet.size());
    if constexpr (dim == 3) {
        //TODO: parallelize
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            const VECTOR<int, dim + 1>& cIVInd = constraintSet[cI];
            Eigen::Matrix<T, dim, 1> relDX3D;
            if (cIVInd[0] >= 0) {
                // ++++ edge-edge, no mollified stencil for friction
                const VECTOR<T, 3>& dXea0 = std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 3>& dXea1 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 3>& dXeb0 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                const VECTOR<T, 3>& dXeb1 = std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                Eigen::Matrix<T, 3, 1> dea0(dXea0.data), dea1(dXea1.data), deb0(dXeb0.data), deb1(dXeb1.data);

                Edge_Edge_RelDX(dea0, dea1, deb0, deb1, closestPoint[cI][0], closestPoint[cI][1], relDX3D);
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(cIVInd[1] >= 0);

                const VECTOR<T, 3>& dXp = std::get<0>(dX.Get_Unchecked(-cIVInd[0] - 1));
                Eigen::Matrix<T, 3, 1> dp(dXp.data);
                if (cIVInd[2] < 0) {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& dXp1 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> dp1(dXp1.data);

                    Point_Point_RelDX(dp, dp1, relDX3D);
                }
                else if (cIVInd[3] < 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& dXe0 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXe1 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> de0(dXe0.data), de1(dXe1.data);

                    Point_Edge_RelDX(dp, de0, de1, closestPoint[cI][0], relDX3D);
                }
                else {
                    // -+++ PT 
                    const VECTOR<T, 3>& dXt0 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXt1 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& dXt2 = std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> dt0(dXt0.data), dt1(dXt1.data), dt2(dXt2.data);

                    Point_Triangle_RelDX(dp, dt0, dt1, dt2, closestPoint[cI][0], closestPoint[cI][1], relDX3D);
                }
            }

            f0_SF((relDX3D.transpose() * tanBasis[cI]).squaredNorm(), epsvh, EI[cI]);
            EI[cI] *= normalForce[cI];

            if (cIVInd[3] < -1) {
                EI[cI] *= -cIVInd[3];
            }
        }
    }
    else {
        //TODO:
    }
    E += mu * std::accumulate(EI.begin(), EI.end(), T(0));
}

template <class T, int dim = 3>
void Compute_Friction_Gradient(
    MESH_NODE<T, dim>& X, MESH_NODE<T, dim>& Xn,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    const std::vector<Eigen::Matrix<T, dim - 1, 1>>& closestPoint,
    const std::vector<Eigen::Matrix<T, dim, dim - 1>>& tanBasis,
    const std::vector<T>& normalForce,
    T epsvh2, T mu, MESH_NODE_ATTR<T, dim>& nodeAttr)
{
    TIMER_FLAG("Compute_Friction_Gradient");

    T epsvh = std::sqrt(epsvh2);

    MESH_NODE<T, dim> dX(X.size);
    X.deep_copy_to(dX);
    dX.Join(Xn).Par_Each([&](int id, auto data) {
        auto &[dx, xn] = data;
        dx -= xn;
    });

    if constexpr (dim == 3) {
        //TODO: parallelize
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            VECTOR<int, dim + 1> cIVInd = constraintSet[cI]; //NOTE: copy to be able to modify in the loop if needed
            Eigen::Matrix<T, dim, 1> relDX3D;
            if (cIVInd[0] >= 0) {
                // ++++ edge-edge, no mollified stencil for friction
                const VECTOR<T, 3>& dXea0 = std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 3>& dXea1 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 3>& dXeb0 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                const VECTOR<T, 3>& dXeb1 = std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                Eigen::Matrix<T, 3, 1> dea0(dXea0.data), dea1(dXea1.data), deb0(dXeb0.data), deb1(dXeb1.data);

                Edge_Edge_RelDX(dea0, dea1, deb0, deb1, closestPoint[cI][0], closestPoint[cI][1], relDX3D);

                Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
                T f1_div_relDXNorm;
                f1_SF_Div_RelDXNorm(relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                relDX *= f1_div_relDXNorm * mu * normalForce[cI];

                Eigen::Matrix<T, 12, 1> TTTDX;
                Edge_Edge_RelDXTan_To_Mesh(relDX, tanBasis[cI], closestPoint[cI][0], closestPoint[cI][1], TTTDX);

                for (int i = 0; i < 4; ++i) {
                    VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(cIVInd[i]));
                    g += TTTDX.data() + i * dim;
                }
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(cIVInd[1] >= 0);

                cIVInd[0] = -cIVInd[0] - 1;
                const VECTOR<T, 3>& dXp = std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                Eigen::Matrix<T, 3, 1> dp(dXp.data);
                if (cIVInd[2] < 0) {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& dXp1 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> dp1(dXp1.data);

                    Point_Point_RelDX(dp, dp1, relDX3D);

                    Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
                    T f1_div_relDXNorm;
                    f1_SF_Div_RelDXNorm(relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm * -cIVInd[3] * mu * normalForce[cI];

                    Eigen::Matrix<T, 6, 1> TTTDX;
                    Point_Point_RelDXTan_To_Mesh(relDX, tanBasis[cI], TTTDX);

                    for (int i = 0; i < 2; ++i) {
                        VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(cIVInd[i]));
                        g += TTTDX.data() + i * dim;
                    }
                }
                else if (cIVInd[3] < 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& dXe0 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXe1 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> de0(dXe0.data), de1(dXe1.data);

                    Point_Edge_RelDX(dp, de0, de1, closestPoint[cI][0], relDX3D);

                    Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
                    T f1_div_relDXNorm;
                    f1_SF_Div_RelDXNorm(relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm * -cIVInd[3] * mu * normalForce[cI];

                    Eigen::Matrix<T, 9, 1> TTTDX;
                    Point_Edge_RelDXTan_To_Mesh(relDX, tanBasis[cI], closestPoint[cI][0], TTTDX);

                    for (int i = 0; i < 3; ++i) {
                        VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(cIVInd[i]));
                        g += TTTDX.data() + i * dim;
                    }
                }
                else {
                    // -+++ PT 
                    const VECTOR<T, 3>& dXt0 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXt1 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& dXt2 = std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> dt0(dXt0.data), dt1(dXt1.data), dt2(dXt2.data);

                    Point_Triangle_RelDX(dp, dt0, dt1, dt2, closestPoint[cI][0], closestPoint[cI][1], relDX3D);

                    Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
                    T f1_div_relDXNorm;
                    f1_SF_Div_RelDXNorm(relDX.squaredNorm(), epsvh, f1_div_relDXNorm);
                    relDX *= f1_div_relDXNorm * mu * normalForce[cI];

                    Eigen::Matrix<T, 12, 1> TTTDX;
                    Point_Triangle_RelDXTan_To_Mesh(relDX, tanBasis[cI], 
                        closestPoint[cI][0], closestPoint[cI][1], TTTDX);

                    for (int i = 0; i < 4; ++i) {
                        VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(cIVInd[i]));
                        g += TTTDX.data() + i * dim;
                    }
                }
            }
        }
    }
    else {
        //TODO
    }
}

template <class T, int dim = 3>
void Compute_Friction_Hessian(
    MESH_NODE<T, dim>& X, MESH_NODE<T, dim>& Xn,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    const std::vector<Eigen::Matrix<T, dim - 1, 1>>& closestPoint,
    const std::vector<Eigen::Matrix<T, dim, dim - 1>>& tanBasis,
    const std::vector<T>& normalForce,
    T epsvh2, T mu, std::vector<Eigen::Triplet<T>>& triplets)
{
    TIMER_FLAG("Compute_Friction_Hessian");

    T epsvh = std::sqrt(epsvh2);

    MESH_NODE<T, dim> dX(X.size);
    X.deep_copy_to(dX);
    dX.Join(Xn).Par_Each([&](int id, auto data) {
        auto &[dx, xn] = data;
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
            }
            else if (cIVInd[2] >= 0) {
                // PE, 9x9
                curStartInd += 81;    
            }
            else {
                // PP, 6x6
                curStartInd += 36;
            }
        }
        triplets.resize(curStartInd);

        threads.Par_Each([&](int cI, auto data) {
            const auto &[tripletStart] = data;
            VECTOR<int, 4> cIVInd = constraintSet[cI]; //NOTE: copy to be able to modify in the loop if needed
            assert(cIVInd[1] >= 0);

            Eigen::Matrix<T, dim, 1> relDX3D;
            if (cIVInd[0] >= 0) {
                // ++++ edge-edge, no mollified stencil for friction
                const VECTOR<T, 3>& dXea0 = std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 3>& dXea1 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 3>& dXeb0 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                const VECTOR<T, 3>& dXeb1 = std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                Eigen::Matrix<T, 3, 1> dea0(dXea0.data), dea1(dXea1.data), deb0(dXeb0.data), deb1(dXeb1.data);

                Edge_Edge_RelDX(dea0, dea1, deb0, deb1, closestPoint[cI][0], closestPoint[cI][1], relDX3D);
                Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
                T relDXSqNorm = relDX.squaredNorm();
                T relDXNorm = std::sqrt(relDXSqNorm);

                Eigen::Matrix<T, 2, 12> TT;
                Edge_Edge_TT(tanBasis[cI], closestPoint[cI][0], closestPoint[cI][1], TT);

                T f1_div_relDXNorm, f2_term;
                f1_SF_Div_RelDXNorm(relDXSqNorm, epsvh, f1_div_relDXNorm);
                f2_SF_Term(relDXSqNorm, epsvh, f2_term);

                Eigen::Matrix<T, 12, 12> HessianI;
                if (relDXSqNorm >= epsvh2) {
                    // no SPD projection needed
                    Eigen::Matrix<T, 2, 1> ubar(-relDX[1], relDX[0]);
                    HessianI = (TT.transpose() * ((mu * normalForce[cI] * f1_div_relDXNorm / relDXSqNorm) * ubar)) * (ubar.transpose() * TT);
                }
                else {
                    if (relDXNorm == 0) {
                        // no SPD projection needed
                        HessianI = ((mu * normalForce[cI] * f1_div_relDXNorm) * TT.transpose()) * TT;
                    }
                    else {
                        // only need to project the inner 2x2 matrix to SPD
                        Eigen::Matrix<T, 2, 2> innerMtr = ((f2_term / relDXNorm) * relDX) * relDX.transpose();
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
                                triplets[tripletStart + (i * 3 + idI) * 12 + j * 3 + jdI] = Eigen::Triplet<T>(
                                    cIVInd[i] * 3 + idI, cIVInd[j] * 3 + jdI,
                                    HessianI(i * 3 + idI, j * 3 + jdI));
                            }
                        }
                    }
                }
            }
            else {
                // point-triangle and degenerate edge-edge
                assert(cIVInd[1] >= 0);

                cIVInd[0] = -cIVInd[0] - 1;
                const VECTOR<T, 3>& dXp = std::get<0>(dX.Get_Unchecked(cIVInd[0]));
                Eigen::Matrix<T, 3, 1> dp(dXp.data);
                if (cIVInd[2] < 0) {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& dXp1 = std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> dp1(dXp1.data);

                    Point_Point_RelDX(dp, dp1, relDX3D);
                    Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
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
                        HessianI = (TT.transpose() * ((-cIVInd[3] * mu * normalForce[cI] * f1_div_relDXNorm / relDXSqNorm) * ubar)) * (ubar.transpose() * TT);
                    }
                    else {
                        if (relDXNorm == 0) {
                            // no SPD projection needed
                            HessianI = ((-cIVInd[3] * mu * normalForce[cI] * f1_div_relDXNorm) * TT.transpose()) * TT;
                        }
                        else {
                            // only need to project the inner 2x2 matrix to SPD
                            Eigen::Matrix<T, 2, 2> innerMtr = ((f2_term / relDXNorm) * relDX) * relDX.transpose();
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
                                    triplets[tripletStart + (i * 3 + idI) * 6 + j * 3 + jdI] = Eigen::Triplet<T>(cIVInd[i] * 3 + idI, cIVInd[j] * 3 + jdI,
                                        HessianI(i * 3 + idI, j * 3 + jdI));
                                }
                            }
                        }
                    }
                }
                else if (cIVInd[3] < 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& dXe0 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXe1 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> de0(dXe0.data), de1(dXe1.data);

                    Point_Edge_RelDX(dp, de0, de1, closestPoint[cI][0], relDX3D);
                    Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
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
                        HessianI = (TT.transpose() * ((-cIVInd[3] * mu * normalForce[cI] * f1_div_relDXNorm / relDXSqNorm) * ubar)) * (ubar.transpose() * TT);
                    }
                    else {
                        if (relDXNorm == 0) {
                            // no SPD projection needed
                            HessianI = ((-cIVInd[3] * mu * normalForce[cI] * f1_div_relDXNorm) * TT.transpose()) * TT;
                        }
                        else {
                            // only need to project the inner 2x2 matrix to SPD
                            Eigen::Matrix<T, 2, 2> innerMtr = ((f2_term / relDXNorm) * relDX) * relDX.transpose();
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
                                    triplets[tripletStart + (i * 3 + idI) * 9 + j * 3 + jdI] = Eigen::Triplet<T>(cIVInd[i] * 3 + idI, cIVInd[j] * 3 + jdI,
                                        HessianI(i * 3 + idI, j * 3 + jdI));
                                }
                            }
                        }
                    }
                }
                else {
                    // -+++ PT 
                    const VECTOR<T, 3>& dXt0 = std::get<0>(dX.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& dXt1 = std::get<0>(dX.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& dXt2 = std::get<0>(dX.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> dt0(dXt0.data), dt1(dXt1.data), dt2(dXt2.data);

                    Point_Triangle_RelDX(dp, dt0, dt1, dt2, closestPoint[cI][0], closestPoint[cI][1], relDX3D);
                    Eigen::Matrix<T, dim - 1, 1> relDX = tanBasis[cI].transpose() * relDX3D;
                    T relDXSqNorm = relDX.squaredNorm();
                    T relDXNorm = std::sqrt(relDXSqNorm);

                    Eigen::Matrix<T, 2, 12> TT;
                    Point_Triangle_TT(tanBasis[cI], closestPoint[cI][0], closestPoint[cI][1], TT);

                    T f1_div_relDXNorm, f2_term;
                    f1_SF_Div_RelDXNorm(relDXSqNorm, epsvh, f1_div_relDXNorm);
                    f2_SF_Term(relDXSqNorm, epsvh, f2_term);

                    Eigen::Matrix<T, 12, 12> HessianI;
                    if (relDXSqNorm >= epsvh2) {
                        // no SPD projection needed
                        Eigen::Matrix<T, 2, 1> ubar(-relDX[1], relDX[0]);
                        HessianI = (TT.transpose() * ((mu * normalForce[cI] * f1_div_relDXNorm / relDXSqNorm) * ubar)) * (ubar.transpose() * TT);
                    }
                    else {
                        if (relDXNorm == 0) {
                            // no SPD projection needed
                            HessianI = ((mu * normalForce[cI] * f1_div_relDXNorm) * TT.transpose()) * TT;
                        }
                        else {
                            // only need to project the inner 2x2 matrix to SPD
                            Eigen::Matrix<T, 2, 2> innerMtr = ((f2_term / relDXNorm) * relDX) * relDX.transpose();
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
                                    triplets[tripletStart + (i * 3 + idI) * 12 + j * 3 + jdI] = Eigen::Triplet<T>(cIVInd[i] * 3 + idI, cIVInd[j] * 3 + jdI,
                                        HessianI(i * 3 + idI, j * 3 + jdI));
                                }
                            }
                        }
                    }
                }
            }
        });
    }
    else {
        //TODO
    }
}

}