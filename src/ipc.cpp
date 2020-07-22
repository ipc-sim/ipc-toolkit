#include "ipc.hpp"

#include <barrier/barrier.hpp>
#include <spatial_hash/hash_grid.hpp>

#include <distance/edge_edge.hpp>
#include <distance/edge_edge_mollifier.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_triangle.hpp>

namespace ipc {

void compute_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat_squared,
    ccd::Candidates& constraint_set)
{
    double dhat = std::sqrt(dhat_squared);

    ccd::HashGrid hash_grid;
    hash_grid.resize(V, V, E, /*inflation_radius=*/dhat);

    hash_grid.addVertices(V, V, /*inflation_radius=*/dhat);
    hash_grid.addEdges(V, V, E, /*inflation_radius=*/dhat);
    if (V.cols() == 3) {
        // These are not needed for 2D
        hash_grid.addFaces(V, V, F, /*inflation_radius=*/dhat);
    }

    if (V.cols() == 2) {
        // This is not needed for 3D
        hash_grid.getVertexEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), constraint_set.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), constraint_set.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, /*group_ids=*/Eigen::VectorXi(), constraint_set.fv_candidates);
    }

    // TODO: Cull the candidates by measuring the distance and dropping those
    // that are greater than dhat.
}

double compute_barrier_potential(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared,
    double barrier_stiffness)
{
    double potential = 0;

    for (const auto& ev_candidate : constraint_set.ev_candidates) {
        double distance = point_edge_distance(
            V.row(ev_candidate.vertex_index),
            V.row(E(ev_candidate.edge_index, 0)),
            V.row(E(ev_candidate.edge_index, 1)));
        potential += barrier(distance, dhat_squared);
    }

    for (const auto& ee_candidate : constraint_set.ee_candidates) {
        const Eigen::Vector3d& ea0 = V.row(E(ee_candidate.edge0_index, 0));
        const Eigen::Vector3d& ea1 = V.row(E(ee_candidate.edge0_index, 1));
        const Eigen::Vector3d& eb0 = V.row(E(ee_candidate.edge1_index, 0));
        const Eigen::Vector3d& eb1 = V.row(E(ee_candidate.edge1_index, 1));

        double distance = edge_edge_distance(ea0, ea1, eb0, eb1);
        double eps_x = edge_edge_mollifier_threshold(
            V_rest.row(E(ee_candidate.edge0_index, 0)),
            V_rest.row(E(ee_candidate.edge0_index, 1)),
            V_rest.row(E(ee_candidate.edge1_index, 0)),
            V_rest.row(E(ee_candidate.edge1_index, 1)));
        potential += edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x)
            * barrier(distance, dhat_squared);
    }

    for (const auto& fv_candidate : constraint_set.fv_candidates) {
        const Eigen::Vector3d& p = V.row(fv_candidate.vertex_index);
        const Eigen::Vector3d& t0 = V.row(E(fv_candidate.face_index, 0));
        const Eigen::Vector3d& t1 = V.row(E(fv_candidate.face_index, 1));
        const Eigen::Vector3d& t2 = V.row(E(fv_candidate.face_index, 2));

        double distance = point_triangle_distance(p, t0, t1, t2);
        potential += barrier(distance, dhat_squared);
    }

    return barrier_stiffness * potential;
}

/*
template <class T, int dim>
void Compute_Barrier_Gradient(
    MESH_NODE<T, dim>& X,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    T dHat2,
    T kappa[],
    MESH_NODE_ATTR<T, dim>& nodeAttr)
{
    TIMER_FLAG("Compute_Barrier_Gradient");

    if constexpr (dim == 2) {
        // TODO: parallelize (loop contains write conflict!)
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            const VECTOR<int, 3>& cIVInd = constraintSet[cI];
            if (cIVInd[2] < 0) {
                // PP
                const VECTOR<T, 2>& Xp0 =
                    std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 2>& Xp1 =
                    std::get<0>(X.Get_Unchecked(cIVInd[1]));
                Eigen::Matrix<T, 2, 1> p0(Xp0.data), p1(Xp1.data);

                T dist2;
                Point_Point_Distance(p0, p1, dist2);
                Eigen::Matrix<T, 4, 1> distGrad;
                Point_Point_Distance_Gradient(p0, p1, distGrad);

                T barrierGrad;
                Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);

                VECTOR<T, 2>& g0 = std::get<FIELDS<MESH_NODE_ATTR<T, 2>>::g>(
                    nodeAttr.Get_Unchecked(cIVInd[0]));
                VECTOR<T, 2>& g1 = std::get<FIELDS<MESH_NODE_ATTR<T, 2>>::g>(
                    nodeAttr.Get_Unchecked(cIVInd[1]));
                distGrad *= barrierGrad;
                g0 += distGrad.data();
                g1 += distGrad.data() + 2;
            } else {
                // PE
                const VECTOR<T, 2>& Xp =
                    std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 2>& Xe0 =
                    std::get<0>(X.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 2>& Xe1 =
                    std::get<0>(X.Get_Unchecked(cIVInd[2]));
                Eigen::Matrix<T, 2, 1> p(Xp.data), e0(Xe0.data), e1(Xe1.data);

                T dist2;
                Point_Edge_Distance(p, e0, e1, dist2);
                Eigen::Matrix<T, 6, 1> distGrad;
                Point_Edge_Distance_Gradient(p, e0, e1, distGrad);

                T barrierGrad;
                Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);

                VECTOR<T, 2>& g0 = std::get<FIELDS<MESH_NODE_ATTR<T, 2>>::g>(
                    nodeAttr.Get_Unchecked(cIVInd[0]));
                VECTOR<T, 2>& g1 = std::get<FIELDS<MESH_NODE_ATTR<T, 2>>::g>(
                    nodeAttr.Get_Unchecked(cIVInd[1]));
                VECTOR<T, 2>& g2 = std::get<FIELDS<MESH_NODE_ATTR<T, 2>>::g>(
                    nodeAttr.Get_Unchecked(cIVInd[2]));
                distGrad *= barrierGrad;
                g0 += distGrad.data();
                g1 += distGrad.data() + 2;
                g2 += distGrad.data() + 4;
            }
        }
    } else {
        // TODO: parallelize (loop contains write conflict!)
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            const VECTOR<int, 4>& cIVInd = constraintSet[cI];
            assert(cIVInd[1] >= 0);
            if (cIVInd[0] >= 0) {
                // EE
                if (cIVInd[3] >= 0 && cIVInd[2] >= 0) {
                    // ++++ EE, no mollification
                    const VECTOR<T, 3>& Xea0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[0]));
                    const VECTOR<T, 3>& Xea1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xeb0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xeb1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                        eb0(Xeb0.data), eb1(Xeb1.data);

                    T dist2;
                    Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2);
                    Eigen::Matrix<T, 12, 1> distGrad;
                    Edge_Edge_Distance_Gradient(ea0, ea1, eb0, eb1, distGrad);

                    T barrierGrad;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);

                    VECTOR<T, 3>& g0 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[0]));
                    VECTOR<T, 3>& g1 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[1]));
                    VECTOR<T, 3>& g2 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[2]));
                    VECTOR<T, 3>& g3 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[3]));
                    distGrad *= barrierGrad;
                    g0 += distGrad.data();
                    g1 += distGrad.data() + 3;
                    g2 += distGrad.data() + 6;
                    g3 += distGrad.data() + 9;
                } else {
                    // EE, PE, or PP with mollification
                    if (cIVInd[3] >= 0) {
                        // ++-+ EE with mollification
                        std::array<int, 4> edgeVInd = { cIVInd[0], cIVInd[1],
                                                        -cIVInd[2] - 1,
                                                        cIVInd[3] };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                            eb0(Xeb0.data), eb1(Xeb1.data);

                        T dist2;
                        Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2);
                        Eigen::Matrix<T, 12, 1> distGrad;
                        Edge_Edge_Distance_Gradient(
                            ea0, ea1, eb0, eb1, distGrad);

                        T b, bGrad;
                        Barrier(dist2, dHat2, kappa, b);
                        Barrier_Gradient(dist2, dHat2, kappa, bGrad);

                        const VECTOR<T, 3>& Xea0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0_rest(Xea0_rest.data),
                            ea1_rest(Xea1_rest.data), eb0_rest(Xeb0_rest.data),
                            eb1_rest(Xeb1_rest.data);
                        T eps_x, e;
                        Edge_Edge_Mollifier_Threshold(
                            ea0_rest, ea1_rest, eb0_rest, eb1_rest, eps_x);
                        Edge_Edge_Mollifier(ea0, ea1, eb0, eb1, eps_x, e);
                        Eigen::Matrix<T, 12, 1> eGrad;
                        Edge_Edge_Mollifier_Gradient(
                            ea0, ea1, eb0, eb1, eps_x, eGrad);

                        VECTOR<T, 3>& g0 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        VECTOR<T, 3>& g1 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        VECTOR<T, 3>& g2 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        VECTOR<T, 3>& g3 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        eGrad = (e * bGrad) * distGrad + (b)*eGrad;
                        g0 += eGrad.data();
                        g1 += eGrad.data() + 3;
                        g2 += eGrad.data() + 6;
                        g3 += eGrad.data() + 9;
                    } else if (cIVInd[2] >= 0) {
                        // +++- PE with mollification, multiplicity 1
                        std::array<int, 4> edgeVInd = { cIVInd[0],
                                                        -cIVInd[3] - 1,
                                                        cIVInd[1], cIVInd[2] };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                            eb0(Xeb0.data), eb1(Xeb1.data);

                        T dist2;
                        Point_Edge_Distance(ea0, eb0, eb1, dist2);
                        Eigen::Matrix<T, 9, 1> distGrad;
                        Point_Edge_Distance_Gradient(ea0, eb0, eb1, distGrad);

                        T b, bGrad;
                        Barrier(dist2, dHat2, kappa, b);
                        Barrier_Gradient(dist2, dHat2, kappa, bGrad);

                        const VECTOR<T, 3>& Xea0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0_rest(Xea0_rest.data),
                            ea1_rest(Xea1_rest.data), eb0_rest(Xeb0_rest.data),
                            eb1_rest(Xeb1_rest.data);
                        T eps_x, e;
                        Edge_Edge_Mollifier_Threshold(
                            ea0_rest, ea1_rest, eb0_rest, eb1_rest, eps_x);
                        Edge_Edge_Mollifier(ea0, ea1, eb0, eb1, eps_x, e);
                        Eigen::Matrix<T, 12, 1> eGrad;
                        Edge_Edge_Mollifier_Gradient(
                            ea0, ea1, eb0, eb1, eps_x, eGrad);

                        VECTOR<T, 3>& g0 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        VECTOR<T, 3>& g1 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        VECTOR<T, 3>& g2 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        VECTOR<T, 3>& g3 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        eGrad *= b;
                        distGrad *= e * bGrad;
                        g0 += eGrad.data();
                        g1 += eGrad.data() + 3;
                        g2 += eGrad.data() + 6;
                        g3 += eGrad.data() + 9;
                        g0 += distGrad.data();
                        g2 += distGrad.data() + 3;
                        g3 += distGrad.data() + 6;
                    } else {
                        // ++-- PP with mollification, multiplicity 1
                        std::array<int, 4> edgeVInd = {
                            cIVInd[0], -cIVInd[2] - 1, cIVInd[1], -cIVInd[3] - 1
                        };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                            eb0(Xeb0.data), eb1(Xeb1.data);

                        T dist2;
                        Point_Point_Distance(ea0, eb0, dist2);
                        Eigen::Matrix<T, 6, 1> distGrad;
                        Point_Point_Distance_Gradient(ea0, eb0, distGrad);

                        T b, bGrad;
                        Barrier(dist2, dHat2, kappa, b);
                        Barrier_Gradient(dist2, dHat2, kappa, bGrad);

                        const VECTOR<T, 3>& Xea0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0_rest(Xea0_rest.data),
                            ea1_rest(Xea1_rest.data), eb0_rest(Xeb0_rest.data),
                            eb1_rest(Xeb1_rest.data);
                        T eps_x, e;
                        Edge_Edge_Mollifier_Threshold(
                            ea0_rest, ea1_rest, eb0_rest, eb1_rest, eps_x);
                        Edge_Edge_Mollifier(ea0, ea1, eb0, eb1, eps_x, e);
                        Eigen::Matrix<T, 12, 1> eGrad;
                        Edge_Edge_Mollifier_Gradient(
                            ea0, ea1, eb0, eb1, eps_x, eGrad);

                        VECTOR<T, 3>& g0 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        VECTOR<T, 3>& g1 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        VECTOR<T, 3>& g2 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        VECTOR<T, 3>& g3 =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        eGrad *= b;
                        distGrad *= e * bGrad;
                        g0 += eGrad.data();
                        g1 += eGrad.data() + 3;
                        g2 += eGrad.data() + 6;
                        g3 += eGrad.data() + 9;
                        g0 += distGrad.data();
                        g2 += distGrad.data() + 3;
                    }
                }
            } else {
                // PT, PE, and PP
                if (cIVInd[3] >= 0) {
                    // -+++ PT
                    assert(cIVInd[2] >= 0);
                    const VECTOR<T, 3>& Xp =
                        std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                    const VECTOR<T, 3>& Xt0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xt1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xt2 =
                        std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> p(Xp.data), t0(Xt0.data),
                        t1(Xt1.data), t2(Xt2.data);

                    T dist2;
                    Point_Triangle_Distance(p, t0, t1, t2, dist2);
                    Eigen::Matrix<T, 12, 1> distGrad;
                    Point_Triangle_Distance_Gradient(p, t0, t1, t2, distGrad);

                    T barrierGrad;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);

                    VECTOR<T, 3>& g0 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(-cIVInd[0] - 1));
                    VECTOR<T, 3>& g1 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[1]));
                    VECTOR<T, 3>& g2 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[2]));
                    VECTOR<T, 3>& g3 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[3]));
                    distGrad *= barrierGrad;
                    g0 += distGrad.data();
                    g1 += distGrad.data() + 3;
                    g2 += distGrad.data() + 6;
                    g3 += distGrad.data() + 9;
                } else if (cIVInd[2] >= 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp =
                        std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                    const VECTOR<T, 3>& Xe0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xe1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> p(Xp.data), e0(Xe0.data),
                        e1(Xe1.data);

                    T dist2;
                    Point_Edge_Distance(p, e0, e1, dist2);
                    Eigen::Matrix<T, 9, 1> distGrad;
                    Point_Edge_Distance_Gradient(p, e0, e1, distGrad);

                    T barrierGrad;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);

                    VECTOR<T, 3>& g0 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(-cIVInd[0] - 1));
                    VECTOR<T, 3>& g1 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[1]));
                    VECTOR<T, 3>& g2 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[2]));
                    distGrad *= -cIVInd[3] * barrierGrad;
                    g0 += distGrad.data();
                    g1 += distGrad.data() + 3;
                    g2 += distGrad.data() + 6;
                } else {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp0 =
                        std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                    const VECTOR<T, 3>& Xp1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> p0(Xp0.data), p1(Xp1.data);

                    T dist2;
                    Point_Point_Distance(p0, p1, dist2);
                    Eigen::Matrix<T, 6, 1> distGrad;
                    Point_Point_Distance_Gradient(p0, p1, distGrad);

                    T barrierGrad;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);

                    VECTOR<T, 3>& g0 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(-cIVInd[0] - 1));
                    VECTOR<T, 3>& g1 =
                        std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::g>(
                            nodeAttr.Get_Unchecked(cIVInd[1]));
                    distGrad *= -cIVInd[3] * barrierGrad;
                    g0 += distGrad.data();
                    g1 += distGrad.data() + 3;
                }
            }
        }
    }
}

template <class T, int dim>
void Compute_Barrier_Hessian(
    MESH_NODE<T, dim>& X,
    MESH_NODE_ATTR<T, dim>& nodeAttr,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    T dHat2,
    T kappa[],
    std::vector<Eigen::Triplet<T>>& triplets)
{
    TIMER_FLAG("Compute_Barrier_Hessian");

    if constexpr (dim == 2) {
        BASE_STORAGE<int> threads(constraintSet.size());
        int curStartInd = triplets.size();
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            threads.Append(curStartInd);
            const VECTOR<int, 3>& cIVInd = constraintSet[cI];
            if (cIVInd[2] < 0) {
                // PP, 4x4
                curStartInd += 16;
            } else {
                // PE, 6x6
                curStartInd += 36;
            }
        }
        triplets.resize(curStartInd);

        threads.Par_Each([&](int cI, auto data) {
            const auto& [tripletStart] = data;
            const VECTOR<int, 3>& cIVInd = constraintSet[cI];
            if (cIVInd[2] < 0) {
                // PP
                const VECTOR<T, 2>& Xp0 =
                    std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 2>& Xp1 =
                    std::get<0>(X.Get_Unchecked(cIVInd[1]));
                Eigen::Matrix<T, 2, 1> p0(Xp0.data), p1(Xp1.data);

                T dist2;
                Point_Point_Distance(p0, p1, dist2);
                Eigen::Matrix<T, 4, 1> distGrad;
                Point_Point_Distance_Gradient(p0, p1, distGrad);
                Eigen::Matrix<T, 4, 4> distH;
                Point_Point_Distance_Hessian(p0, p1, distH);

                T barrierGrad, barrierH;
                Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);
                Barrier_Hessian(dist2, dHat2, kappa, barrierH);

                Eigen::Matrix<T, 4, 4> HessianI =
                    barrierH * distGrad * distGrad.transpose()
                    + barrierGrad * distH;
                makePD(HessianI);
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        for (int idI = 0; idI < 2; ++idI) {
                            for (int jdI = 0; jdI < 2; ++jdI) {
                                triplets
                                    [tripletStart + (i * 2 + idI) * 4 + j * 2
                                     + jdI] =
                                        Eigen::Triplet<T>(
                                            cIVInd[i] * 2 + idI,
                                            cIVInd[j] * 2 + jdI,
                                            HessianI(i * 2 + idI, j * 2 + jdI));
                            }
                        }
                    }
                }
            } else {
                // PE
                const VECTOR<T, 2>& Xp =
                    std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 2>& Xe0 =
                    std::get<0>(X.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 2>& Xe1 =
                    std::get<0>(X.Get_Unchecked(cIVInd[2]));
                Eigen::Matrix<T, 2, 1> p(Xp.data), e0(Xe0.data), e1(Xe1.data);

                T dist2;
                Point_Edge_Distance(p, e0, e1, dist2);
                Eigen::Matrix<T, 6, 1> distGrad;
                Point_Edge_Distance_Gradient(p, e0, e1, distGrad);
                Eigen::Matrix<T, 6, 6> distH;
                Point_Edge_Distance_Hessian(p, e0, e1, distH);

                T barrierGrad, barrierH;
                Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);
                Barrier_Hessian(dist2, dHat2, kappa, barrierH);

                Eigen::Matrix<T, 6, 6> HessianI =
                    barrierH * distGrad * distGrad.transpose()
                    + barrierGrad * distH;
                makePD(HessianI);
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        for (int idI = 0; idI < 2; ++idI) {
                            for (int jdI = 0; jdI < 2; ++jdI) {
                                triplets
                                    [tripletStart + (i * 2 + idI) * 6 + j * 2
                                     + jdI] =
                                        Eigen::Triplet<T>(
                                            cIVInd[i] * 2 + idI,
                                            cIVInd[j] * 2 + jdI,
                                            HessianI(i * 2 + idI, j * 2 + jdI));
                            }
                        }
                    }
                }
            }
        });
    } else {
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
            if (cIVInd[0] >= 0) {
                // EE
                if (cIVInd[3] >= 0 && cIVInd[2] >= 0) {
                    // ++++ EE, no mollification
                    const VECTOR<T, 3>& Xea0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[0]));
                    const VECTOR<T, 3>& Xea1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xeb0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xeb1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                        eb0(Xeb0.data), eb1(Xeb1.data);

                    T dist2;
                    Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2);
                    Eigen::Matrix<T, 12, 1> distGrad;
                    Edge_Edge_Distance_Gradient(ea0, ea1, eb0, eb1, distGrad);
                    Eigen::Matrix<T, 12, 12> distH;
                    Edge_Edge_Distance_Hessian(ea0, ea1, eb0, eb1, distH);

                    T barrierGrad, barrierH;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);
                    Barrier_Hessian(dist2, dHat2, kappa, barrierH);

                    Eigen::Matrix<T, 12, 12> HessianI =
                        ((barrierH)*distGrad) * distGrad.transpose()
                        + (barrierGrad)*distH;
                    makePD(HessianI);
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
                } else {
                    // EE, PE, or PP with mollification
                    if (cIVInd[3] >= 0) {
                        // ++-+ EE with mollification
                        std::array<int, 4> edgeVInd = { cIVInd[0], cIVInd[1],
                                                        -cIVInd[2] - 1,
                                                        cIVInd[3] };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                            eb0(Xeb0.data), eb1(Xeb1.data);

                        T dist2;
                        Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2);
                        Eigen::Matrix<T, 12, 1> distGrad;
                        Edge_Edge_Distance_Gradient(
                            ea0, ea1, eb0, eb1, distGrad);
                        Eigen::Matrix<T, 12, 12> distH;
                        Edge_Edge_Distance_Hessian(ea0, ea1, eb0, eb1, distH);

                        T b, bGrad, bH;
                        Barrier(dist2, dHat2, kappa, b);
                        Barrier_Gradient(dist2, dHat2, kappa, bGrad);
                        Barrier_Hessian(dist2, dHat2, kappa, bH);

                        const VECTOR<T, 3>& Xea0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0_rest(Xea0_rest.data),
                            ea1_rest(Xea1_rest.data), eb0_rest(Xeb0_rest.data),
                            eb1_rest(Xeb1_rest.data);
                        T eps_x, e;
                        Edge_Edge_Mollifier_Threshold(
                            ea0_rest, ea1_rest, eb0_rest, eb1_rest, eps_x);
                        Edge_Edge_Mollifier(ea0, ea1, eb0, eb1, eps_x, e);
                        Eigen::Matrix<T, 12, 1> eGrad;
                        Edge_Edge_Mollifier_Gradient(
                            ea0, ea1, eb0, eb1, eps_x, eGrad);
                        Eigen::Matrix<T, 12, 12> eH;
                        Edge_Edge_Mollifier_Hessian(
                            ea0, ea1, eb0, eb1, eps_x, eH);

                        Eigen::Matrix<T, 12, 12> kappa_bGrad_eGradT =
                            (bGrad * distGrad) * eGrad.transpose();
                        Eigen::Matrix<T, 12, 12> HessianI = kappa_bGrad_eGradT
                            + kappa_bGrad_eGradT.transpose() + b * eH
                            + (e * bH * distGrad) * distGrad.transpose()
                            + e * bGrad * distH;
                        makePD(HessianI);
                        for (int i = 0; i < 4; ++i) {
                            for (int j = 0; j < 4; ++j) {
                                for (int idI = 0; idI < 3; ++idI) {
                                    for (int jdI = 0; jdI < 3; ++jdI) {
                                        triplets
                                            [tripletStart + (i * 3 + idI) * 12
                                             + j * 3 + jdI] =
                                                Eigen::Triplet<T>(
                                                    edgeVInd[i] * 3 + idI,
                                                    edgeVInd[j] * 3 + jdI,
                                                    HessianI(
                                                        i * 3 + idI,
                                                        j * 3 + jdI));
                                    }
                                }
                            }
                        }
                    } else if (cIVInd[2] >= 0) {
                        // +++- PE with mollification, multiplicity 1
                        std::array<int, 4> edgeVInd = { cIVInd[0],
                                                        -cIVInd[3] - 1,
                                                        cIVInd[1], cIVInd[2] };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                            eb0(Xeb0.data), eb1(Xeb1.data);

                        T dist2;
                        Point_Edge_Distance(ea0, eb0, eb1, dist2);
                        Eigen::Matrix<T, 9, 1> distGrad;
                        Point_Edge_Distance_Gradient(ea0, eb0, eb1, distGrad);
                        Eigen::Matrix<T, 9, 9> distH;
                        Point_Edge_Distance_Hessian(ea0, eb0, eb1, distH);

                        T b, bGrad, bH;
                        Barrier(dist2, dHat2, kappa, b);
                        Barrier_Gradient(dist2, dHat2, kappa, bGrad);
                        Barrier_Hessian(dist2, dHat2, kappa, bH);

                        const VECTOR<T, 3>& Xea0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0_rest(Xea0_rest.data),
                            ea1_rest(Xea1_rest.data), eb0_rest(Xeb0_rest.data),
                            eb1_rest(Xeb1_rest.data);
                        T eps_x, e;
                        Edge_Edge_Mollifier_Threshold(
                            ea0_rest, ea1_rest, eb0_rest, eb1_rest, eps_x);
                        Edge_Edge_Mollifier(ea0, ea1, eb0, eb1, eps_x, e);
                        Eigen::Matrix<T, 12, 1> eGrad;
                        Edge_Edge_Mollifier_Gradient(
                            ea0, ea1, eb0, eb1, eps_x, eGrad);
                        Eigen::Matrix<T, 12, 12> eH;
                        Edge_Edge_Mollifier_Hessian(
                            ea0, ea1, eb0, eb1, eps_x, eH);

                        Eigen::Matrix<T, 12, 12> HessianI = (b)*eH;
                        Eigen::Matrix<T, 9, 9> kappa_e_bH =
                            ((e * bH) * distGrad) * distGrad.transpose()
                            + (e * bGrad) * distH;
                        HessianI.template block<3, 3>(0, 0) +=
                            kappa_e_bH.template block<3, 3>(0, 0);
                        HessianI.template block<3, 6>(0, 6) +=
                            kappa_e_bH.template block<3, 6>(0, 3);
                        HessianI.template block<6, 3>(6, 0) +=
                            kappa_e_bH.template block<6, 3>(3, 0);
                        HessianI.template block<6, 6>(6, 6) +=
                            kappa_e_bH.template block<6, 6>(3, 3);
                        Eigen::Matrix<T, 9, 12> kappa_bGrad_eGradT =
                            ((bGrad)*distGrad) * eGrad.transpose();
                        HessianI.template block<3, 12>(0, 0) +=
                            kappa_bGrad_eGradT.template block<3, 12>(0, 0);
                        HessianI.template block<6, 12>(6, 0) +=
                            kappa_bGrad_eGradT.template block<6, 12>(3, 0);
                        HessianI.template block<12, 3>(0, 0) +=
                            kappa_bGrad_eGradT.template block<3, 12>(0, 0)
                                .transpose();
                        HessianI.template block<12, 6>(0, 6) +=
                            kappa_bGrad_eGradT.template block<6, 12>(3, 0)
                                .transpose();
                        makePD(HessianI);
                        for (int i = 0; i < 4; ++i) {
                            for (int j = 0; j < 4; ++j) {
                                for (int idI = 0; idI < 3; ++idI) {
                                    for (int jdI = 0; jdI < 3; ++jdI) {
                                        triplets
                                            [tripletStart + (i * 3 + idI) * 12
                                             + j * 3 + jdI] =
                                                Eigen::Triplet<T>(
                                                    edgeVInd[i] * 3 + idI,
                                                    edgeVInd[j] * 3 + jdI,
                                                    HessianI(
                                                        i * 3 + idI,
                                                        j * 3 + jdI));
                                    }
                                }
                            }
                        }
                    } else {
                        // ++-- PP with mollification, multiplicity 1
                        std::array<int, 4> edgeVInd = {
                            cIVInd[0], -cIVInd[2] - 1, cIVInd[1], -cIVInd[3] - 1
                        };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                            eb0(Xeb0.data), eb1(Xeb1.data);

                        T dist2;
                        Point_Point_Distance(ea0, eb0, dist2);
                        Eigen::Matrix<T, 6, 1> distGrad;
                        Point_Point_Distance_Gradient(ea0, eb0, distGrad);
                        Eigen::Matrix<T, 6, 6> distH;
                        Point_Point_Distance_Hessian(ea0, eb0, distH);

                        T b, bGrad, bH;
                        Barrier(dist2, dHat2, kappa, b);
                        Barrier_Gradient(dist2, dHat2, kappa, bGrad);
                        Barrier_Hessian(dist2, dHat2, kappa, bH);

                        const VECTOR<T, 3>& Xea0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1_rest =
                            std::get<FIELDS<MESH_NODE_ATTR<T, 3>>::x0>(
                                nodeAttr.Get_Unchecked(edgeVInd[3]));
                        Eigen::Matrix<T, 3, 1> ea0_rest(Xea0_rest.data),
                            ea1_rest(Xea1_rest.data), eb0_rest(Xeb0_rest.data),
                            eb1_rest(Xeb1_rest.data);
                        T eps_x, e;
                        Edge_Edge_Mollifier_Threshold(
                            ea0_rest, ea1_rest, eb0_rest, eb1_rest, eps_x);
                        Edge_Edge_Mollifier(ea0, ea1, eb0, eb1, eps_x, e);
                        Eigen::Matrix<T, 12, 1> eGrad;
                        Edge_Edge_Mollifier_Gradient(
                            ea0, ea1, eb0, eb1, eps_x, eGrad);
                        Eigen::Matrix<T, 12, 12> eH;
                        Edge_Edge_Mollifier_Hessian(
                            ea0, ea1, eb0, eb1, eps_x, eH);

                        Eigen::Matrix<T, 12, 12> HessianI = (b)*eH;
                        Eigen::Matrix<T, 6, 6> kappa_e_bH =
                            ((e * bH) * distGrad) * distGrad.transpose()
                            + (e * bGrad) * distH;
                        HessianI.template block<3, 3>(0, 0) +=
                            kappa_e_bH.template block<3, 3>(0, 0);
                        HessianI.template block<3, 3>(0, 6) +=
                            kappa_e_bH.template block<3, 3>(0, 3);
                        HessianI.template block<3, 3>(6, 0) +=
                            kappa_e_bH.template block<3, 3>(3, 0);
                        HessianI.template block<3, 3>(6, 6) +=
                            kappa_e_bH.template block<3, 3>(3, 3);
                        Eigen::Matrix<T, 6, 12> kappa_bGrad_eGradT =
                            ((bGrad)*distGrad) * eGrad.transpose();
                        HessianI.template block<3, 12>(0, 0) +=
                            kappa_bGrad_eGradT.template block<3, 12>(0, 0);
                        HessianI.template block<3, 12>(6, 0) +=
                            kappa_bGrad_eGradT.template block<3, 12>(3, 0);
                        HessianI.template block<12, 3>(0, 0) +=
                            kappa_bGrad_eGradT.template block<3, 12>(0, 0)
                                .transpose();
                        HessianI.template block<12, 3>(0, 6) +=
                            kappa_bGrad_eGradT.template block<3, 12>(3, 0)
                                .transpose();
                        makePD(HessianI);
                        for (int i = 0; i < 4; ++i) {
                            for (int j = 0; j < 4; ++j) {
                                for (int idI = 0; idI < 3; ++idI) {
                                    for (int jdI = 0; jdI < 3; ++jdI) {
                                        triplets
                                            [tripletStart + (i * 3 + idI) * 12
                                             + j * 3 + jdI] =
                                                Eigen::Triplet<T>(
                                                    edgeVInd[i] * 3 + idI,
                                                    edgeVInd[j] * 3 + jdI,
                                                    HessianI(
                                                        i * 3 + idI,
                                                        j * 3 + jdI));
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                // PT, PE, and PP
                cIVInd[0] = -cIVInd[0] - 1;
                if (cIVInd[3] >= 0) {
                    // -+++ PT
                    assert(cIVInd[2] >= 0);
                    const VECTOR<T, 3>& Xp =
                        std::get<0>(X.Get_Unchecked(cIVInd[0]));
                    const VECTOR<T, 3>& Xt0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xt1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xt2 =
                        std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> p(Xp.data), t0(Xt0.data),
                        t1(Xt1.data), t2(Xt2.data);

                    T dist2;
                    Point_Triangle_Distance(p, t0, t1, t2, dist2);
                    Eigen::Matrix<T, 12, 1> distGrad;
                    Point_Triangle_Distance_Gradient(p, t0, t1, t2, distGrad);
                    Eigen::Matrix<T, 12, 12> distH;
                    Point_Triangle_Distance_Hessian(p, t0, t1, t2, distH);

                    T barrierGrad, barrierH;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);
                    Barrier_Hessian(dist2, dHat2, kappa, barrierH);

                    Eigen::Matrix<T, 12, 12> HessianI =
                        ((barrierH)*distGrad) * distGrad.transpose()
                        + (barrierGrad)*distH;
                    makePD(HessianI);
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
                } else if (cIVInd[2] >= 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp =
                        std::get<0>(X.Get_Unchecked(cIVInd[0]));
                    const VECTOR<T, 3>& Xe0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xe1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> p(Xp.data), e0(Xe0.data),
                        e1(Xe1.data);

                    T dist2;
                    Point_Edge_Distance(p, e0, e1, dist2);
                    Eigen::Matrix<T, 9, 1> distGrad;
                    Point_Edge_Distance_Gradient(p, e0, e1, distGrad);
                    Eigen::Matrix<T, 9, 9> distH;
                    Point_Edge_Distance_Hessian(p, e0, e1, distH);

                    T barrierGrad, barrierH;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);
                    Barrier_Hessian(dist2, dHat2, kappa, barrierH);

                    Eigen::Matrix<T, 9, 9> HessianI =
                        ((-cIVInd[3] * barrierH) * distGrad)
                            * distGrad.transpose()
                        + (-cIVInd[3] * barrierGrad) * distH;
                    makePD(HessianI);
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
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[0]));
                    const VECTOR<T, 3>& Xp1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> p0(Xp0.data), p1(Xp1.data);

                    T dist2;
                    Point_Point_Distance(p0, p1, dist2);
                    Eigen::Matrix<T, 6, 1> distGrad;
                    Point_Point_Distance_Gradient(p0, p1, distGrad);
                    Eigen::Matrix<T, 6, 6> distH;
                    Point_Point_Distance_Hessian(p0, p1, distH);

                    T barrierGrad, barrierH;
                    Barrier_Gradient(dist2, dHat2, kappa, barrierGrad);
                    Barrier_Hessian(dist2, dHat2, kappa, barrierH);

                    Eigen::Matrix<T, 6, 6> HessianI =
                        ((-cIVInd[3] * barrierH) * distGrad)
                            * distGrad.transpose()
                        + (-cIVInd[3] * barrierGrad) * distH;
                    makePD(HessianI);
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
                }
            }
        });
    }
}

template <class T, int dim>
void Compute_Intersection_Free_StepSize(
    MESH_NODE<T, dim>& X,
    const std::vector<int>& boundaryNode,
    const std::vector<VECTOR<int, 2>>& boundaryEdge,
    const std::vector<VECTOR<int, 3>>& boundaryTri,
    const std::vector<T>& searchDir,
    T& stepSize)
{
    TIMER_FLAG("Compute_Intersection_Free_StepSize");
    // TODO: CFL?

#define USE_SH_LFSS
#ifdef USE_SH_LFSS
    SPATIAL_HASH<T, dim> sh;
    {
        TIMER_FLAG("Compute_Intersection_Free_StepSize_Build_Hash");
        sh.Build(
            X, boundaryNode, boundaryEdge, boundaryTri, searchDir, stepSize,
            1.0, 0);
    }
#endif

    if constexpr (dim == 2) {
        BASE_STORAGE<int> threads(boundaryNode.size());
        for (int i = 0; i < boundaryNode.size(); ++i) {
            threads.Append(i);
        }

        std::vector<T> stepSizeNI(boundaryNode.size(), stepSize);
        threads.Par_Each([&](int bNI, auto data) {
            int nI = boundaryNode[bNI];
            const VECTOR<T, 2>& Xp = std::get<0>(X.Get_Unchecked(nI));
            Eigen::Matrix<T, 2, 1> p(Xp.data), dp(searchDir.data() + nI * 2);
#ifdef USE_SH_LFSS
            std::unordered_set<int> edgeInds;
            sh.Query_Point_For_Edges(bNI, edgeInds);
            for (const auto& eInd : edgeInds) {
                const VECTOR<int, 2>& eI = boundaryEdge[eInd];
#else
            for (const auto& eI : boundaryEdge) {
#endif
                if (nI == eI[0] || nI == eI[1]) {
                    continue;
                }

                const VECTOR<T, 2>& Xe0 = std::get<0>(X.Get_Unchecked(eI[0]));
                const VECTOR<T, 2>& Xe1 = std::get<0>(X.Get_Unchecked(eI[1]));
                Eigen::Matrix<T, 2, 1> e0(Xe0.data), e1(Xe1.data);
                Eigen::Matrix<T, 2, 1> de0(searchDir.data() + eI[0] * 2);
                Eigen::Matrix<T, 2, 1> de1(searchDir.data() + eI[1] * 2);

                if (!Point_Edge_CCD_Broadphase(p, e0, e1, dp, de0, de1, 0)) {
                    continue;
                }

                T stepSizeNIEI;
                if (Point_Edge_CCD(
                        p, e0, e1, dp, de0, de1, 0.1, stepSizeNIEI)) {
                    if (stepSizeNI[bNI] > stepSizeNIEI) {
                        stepSizeNI[bNI] = stepSizeNIEI;
                    }
                }
            }
        });
        stepSize = *std::min_element(stepSizeNI.begin(), stepSizeNI.end());
    } else {
        // point-triangle
        {
            TIMER_FLAG("Compute_Intersection_Free_StepSize_PT");
            BASE_STORAGE<int> threadsPT(boundaryNode.size());
            for (int i = 0; i < boundaryNode.size(); ++i) {
                threadsPT.Append(i);
            }

            std::vector<T> largestAlphasPT(boundaryNode.size());
            threadsPT.Par_Each([&](int svI, auto data) {
                int vI = boundaryNode[svI];
                const VECTOR<T, 3>& Xp = std::get<0>(X.Get_Unchecked(vI));
                Eigen::Matrix<T, 3, 1> p(Xp.data),
                    dp(searchDir.data() + vI * dim);
                largestAlphasPT[svI] = 1.0;
#ifdef USE_SH_LFSS
                std::unordered_set<int> sEdgeInds, sTriInds;
                sh.Query_Point_For_Primitives(svI, sEdgeInds, sTriInds);
                for (const auto& sfI : sTriInds)
#else
                for (int sfI = 0; sfI < boundaryTri.size(); ++sfI)
#endif
                {
                    const VECTOR<int, 3>& sfVInd = boundaryTri[sfI];
                    if (!(vI == sfVInd[0] || vI == sfVInd[1]
                          || vI == sfVInd[2])) {
                        const VECTOR<T, 3>& Xt0 =
                            std::get<0>(X.Get_Unchecked(sfVInd[0]));
                        const VECTOR<T, 3>& Xt1 =
                            std::get<0>(X.Get_Unchecked(sfVInd[1]));
                        const VECTOR<T, 3>& Xt2 =
                            std::get<0>(X.Get_Unchecked(sfVInd[2]));
                        Eigen::Matrix<T, 3, 1> t0(Xt0.data), t1(Xt1.data),
                            t2(Xt2.data);
                        Eigen::Matrix<T, 3, 1> dt0(
                            searchDir.data() + sfVInd[0] * dim),
                            dt1(searchDir.data() + sfVInd[1] * dim),
                            dt2(searchDir.data() + sfVInd[2] * dim);

                        if (!Point_Triangle_CCD_Broadphase(
                                p, t0, t1, t2, dp, dt0, dt1, dt2, 0)) {
                            continue;
                        }

                        T largestAlpha = 1.0;
                        if (Point_Triangle_CCD(
                                p, t0, t1, t2, dp, dt0, dt1, dt2, 0.2, 0,
                                largestAlpha)) {
                            if (largestAlphasPT[svI] > largestAlpha) {
                                largestAlphasPT[svI] = largestAlpha;
                            }
                        }
                    }
                }

                // PE
#ifdef USE_SH_LFSS
                for (const auto& eI : sEdgeInds)
#else
                for (int eI = 0; eI < boundaryEdge.size(); ++eI)
#endif
                {
                    const VECTOR<int, 2>& meshEI = boundaryEdge[eI];
                    if (!(vI == meshEI[0] || vI == meshEI[1])) {
                        const VECTOR<T, 3>& Xe0 =
                            std::get<0>(X.Get_Unchecked(meshEI[0]));
                        const VECTOR<T, 3>& Xe1 =
                            std::get<0>(X.Get_Unchecked(meshEI[1]));
                        Eigen::Matrix<T, 3, 1> e0(Xe0.data), e1(Xe1.data);
                        Eigen::Matrix<T, 3, 1> de0(
                            searchDir.data() + meshEI[0] * dim),
                            de1(searchDir.data() + meshEI[1] * dim);

                        if (!Point_Edge_CCD_Broadphase(
                                p, e0, e1, dp, de0, de1, 0)) {
                            continue;
                        }

                        T largestAlpha = 1.0;
                        if (Point_Edge_CCD(
                                p, e0, e1, dp, de0, de1, 0.2, 0,
                                largestAlpha)) {
                            if (largestAlphasPT[svI] > largestAlpha) {
                                largestAlphasPT[svI] = largestAlpha;
                            }
                        }
                    }
                }

                // TODO: PP
            });
            stepSize = std::min(
                stepSize,
                *std::min_element(
                    largestAlphasPT.begin(), largestAlphasPT.end()));
        } // TIMER_FLAG

        // edge-edge
        {
            TIMER_FLAG("Compute_Intersection_Free_StepSize_EE");
            BASE_STORAGE<int> threadsEE(boundaryEdge.size());
            for (int i = 0; i < boundaryEdge.size(); ++i) {
                threadsEE.Append(i);
            }

            std::vector<T> largestAlphasEE(boundaryEdge.size());
            threadsEE.Par_Each([&](int eI, auto data) {
                const VECTOR<int, 2>& meshEI = boundaryEdge[eI];
                const VECTOR<T, 3>& Xea0 =
                    std::get<0>(X.Get_Unchecked(meshEI[0]));
                const VECTOR<T, 3>& Xea1 =
                    std::get<0>(X.Get_Unchecked(meshEI[1]));
                Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data);
                Eigen::Matrix<T, 3, 1> dea0(searchDir.data() + meshEI[0] * dim),
                    dea1(searchDir.data() + meshEI[1] * dim);
                largestAlphasEE[eI] = 1.0;
#ifdef USE_SH_LFSS
                std::unordered_set<int> sEdgeInds;
                sh.Query_Edge_For_Edges(eI, sEdgeInds);
                // NOTE: results may differ when computing step size with large
                // eta as long-distance pairs are dropped
                for (const auto& eJ : sEdgeInds)
#else
                for (int eJ = eI + 1; eJ < boundaryEdge.size(); ++eJ)
#endif
                {
                    const VECTOR<int, 2>& meshEJ = boundaryEdge[eJ];
                    if (!(meshEI[0] == meshEJ[0] || meshEI[0] == meshEJ[1]
                          || meshEI[1] == meshEJ[0] || meshEI[1] == meshEJ[1]
                          || eI > eJ)) {
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(meshEJ[0]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(meshEJ[1]));
                        Eigen::Matrix<T, 3, 1> eb0(Xeb0.data), eb1(Xeb1.data);
                        Eigen::Matrix<T, 3, 1> deb0(
                            searchDir.data() + meshEJ[0] * dim),
                            deb1(searchDir.data() + meshEJ[1] * dim);

                        if (!Edge_Edge_CCD_Broadphase(
                                ea0, ea1, eb0, eb1, dea0, dea1, deb0, deb1,
                                0)) {
                            continue;
                        }

                        T largestAlpha = 1.0;
                        if (Edge_Edge_CCD(
                                ea0, ea1, eb0, eb1, dea0, dea1, deb0, deb1, 0.2,
                                0, largestAlpha)) {
                            if (largestAlphasEE[eI] > largestAlpha) {
                                largestAlphasEE[eI] = largestAlpha;
                            }
                        }
                    }
                }
            });
            stepSize = std::min(
                stepSize,
                *std::min_element(
                    largestAlphasEE.begin(), largestAlphasEE.end()));
        } // TIMER_FLAG
    }
}

template <class T, int dim>
void Compute_Min_Dist2(
    MESH_NODE<T, dim>& X,
    const std::vector<VECTOR<int, dim + 1>>& constraintSet,
    std::vector<T>& dist2,
    T& minDist2)
{
    TIMER_FLAG("Compute_Min_Dist");

    if (constraintSet.empty()) {
        return;
    }

    dist2.resize(constraintSet.size());
    if constexpr (dim == 2) {
        // TODO: parallelize
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            const VECTOR<int, 3>& cIVInd = constraintSet[cI];
            if (cIVInd[2] < 0) {
                // PP
                const VECTOR<T, 2>& Xp0 =
                    std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 2>& Xp1 =
                    std::get<0>(X.Get_Unchecked(cIVInd[1]));
                Eigen::Matrix<T, 2, 1> p0(Xp0[0], Xp0[1]), p1(Xp1[0], Xp1[1]);
                Point_Point_Distance(p0, p1, dist2[cI]);
            } else {
                // PE
                const VECTOR<T, 2>& Xp =
                    std::get<0>(X.Get_Unchecked(cIVInd[0]));
                const VECTOR<T, 2>& Xe0 =
                    std::get<0>(X.Get_Unchecked(cIVInd[1]));
                const VECTOR<T, 2>& Xe1 =
                    std::get<0>(X.Get_Unchecked(cIVInd[2]));
                Eigen::Matrix<T, 2, 1> p(Xp[0], Xp[1]), e0(Xe0[0], Xe0[1]),
                    e1(Xe1[0], Xe1[1]);
                Point_Edge_Distance(p, e0, e1, dist2[cI]);
            }

            if (dist2[cI] <= 0) {
                std::cout << "0 distance detected during barrier evaluation!"
                          << std::endl;
                exit(-1);
            }
        }
    } else {
        // TODO: parallelize
        for (int cI = 0; cI < constraintSet.size(); ++cI) {
            const VECTOR<int, 4>& cIVInd = constraintSet[cI];
            assert(cIVInd[1] >= 0);
            if (cIVInd[0] >= 0) {
                // EE
                if (cIVInd[3] >= 0 && cIVInd[2] >= 0) {
                    // ++++ EE, no mollification
                    const VECTOR<T, 3>& Xea0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[0]));
                    const VECTOR<T, 3>& Xea1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xeb0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xeb1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> ea0(Xea0.data), ea1(Xea1.data),
                        eb0(Xeb0.data), eb1(Xeb1.data);

                    Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2[cI]);
                } else {
                    // EE, PE, or PP with mollification
                    std::array<int, 4> edgeVInd;
                    Eigen::Matrix<T, 3, 1> ea0, ea1, eb0, eb1;
                    if (cIVInd[3] >= 0) {
                        // ++-+ EE with mollification
                        edgeVInd = { cIVInd[0], cIVInd[1], -cIVInd[2] - 1,
                                     cIVInd[3] };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        ea0 = std::move(Eigen::Matrix<T, 3, 1>(Xea0.data));
                        ea1 = std::move(Eigen::Matrix<T, 3, 1>(Xea1.data));
                        eb0 = std::move(Eigen::Matrix<T, 3, 1>(Xeb0.data));
                        eb1 = std::move(Eigen::Matrix<T, 3, 1>(Xeb1.data));

                        Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2[cI]);
                    } else if (cIVInd[2] >= 0) {
                        // +++- PE with mollification, multiplicity 1
                        edgeVInd = { cIVInd[0], -cIVInd[3] - 1, cIVInd[1],
                                     cIVInd[2] };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        ea0 = std::move(Eigen::Matrix<T, 3, 1>(Xea0.data));
                        ea1 = std::move(Eigen::Matrix<T, 3, 1>(Xea1.data));
                        eb0 = std::move(Eigen::Matrix<T, 3, 1>(Xeb0.data));
                        eb1 = std::move(Eigen::Matrix<T, 3, 1>(Xeb1.data));

                        Point_Edge_Distance(ea0, eb0, eb1, dist2[cI]);
                    } else {
                        // ++-- PP with mollification, multiplicity 1
                        edgeVInd = { cIVInd[0], -cIVInd[2] - 1, cIVInd[1],
                                     -cIVInd[3] - 1 };
                        const VECTOR<T, 3>& Xea0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[0]));
                        const VECTOR<T, 3>& Xea1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[1]));
                        const VECTOR<T, 3>& Xeb0 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[2]));
                        const VECTOR<T, 3>& Xeb1 =
                            std::get<0>(X.Get_Unchecked(edgeVInd[3]));
                        ea0 = std::move(Eigen::Matrix<T, 3, 1>(Xea0.data));
                        ea1 = std::move(Eigen::Matrix<T, 3, 1>(Xea1.data));
                        eb0 = std::move(Eigen::Matrix<T, 3, 1>(Xeb0.data));
                        eb1 = std::move(Eigen::Matrix<T, 3, 1>(Xeb1.data));

                        Point_Point_Distance(ea0, eb0, dist2[cI]);
                    }
                }
            } else {
                // PT, PE, and PP
                if (cIVInd[3] >= 0) {
                    // -+++ PT
                    assert(cIVInd[2] >= 0);
                    const VECTOR<T, 3>& Xp =
                        std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                    const VECTOR<T, 3>& Xt0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xt1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    const VECTOR<T, 3>& Xt2 =
                        std::get<0>(X.Get_Unchecked(cIVInd[3]));
                    Eigen::Matrix<T, 3, 1> p(Xp.data), t0(Xt0.data),
                        t1(Xt1.data), t2(Xt2.data);

                    Point_Triangle_Distance(p, t0, t1, t2, dist2[cI]);
                } else if (cIVInd[2] >= 0) {
                    // -++[-] PE, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp =
                        std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                    const VECTOR<T, 3>& Xe0 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    const VECTOR<T, 3>& Xe1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[2]));
                    Eigen::Matrix<T, 3, 1> p(Xp.data), e0(Xe0.data),
                        e1(Xe1.data);

                    Point_Edge_Distance(p, e0, e1, dist2[cI]);
                } else {
                    // -+-[-] PP, last digit stores muliplicity
                    const VECTOR<T, 3>& Xp0 =
                        std::get<0>(X.Get_Unchecked(-cIVInd[0] - 1));
                    const VECTOR<T, 3>& Xp1 =
                        std::get<0>(X.Get_Unchecked(cIVInd[1]));
                    Eigen::Matrix<T, 3, 1> p0(Xp0.data), p1(Xp1.data);

                    Point_Point_Distance(p0, p1, dist2[cI]);
                }
            }
        }
    }
    minDist2 = *std::min_element(dist2.begin(), dist2.end());
}
*/

} // namespace ipc
