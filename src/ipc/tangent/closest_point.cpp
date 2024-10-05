#include "closest_point.hpp"

#include <Eigen/Core>
#include <Eigen/Cholesky>

namespace ipc {

// ============================================================================
// Point - Edge

double point_edge_closest_point(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1)
{
    const VectorMax3d e = e1 - e0;
    return (p - e0).dot(e) / e.squaredNorm();
}

VectorMax9d point_edge_closest_point_jacobian(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1)
{
    const int dim = p.size();
    assert(dim == 2 || dim == 3);
    assert(e0.size() == dim && e1.size() == dim);

    const VectorMax3d e = e1 - e0;
    const VectorMax3d e2p = p - e0;
    const double e_sqnorm = e.squaredNorm();

    VectorMax9d J(3 * dim);
    J.head(dim) = e / e_sqnorm;
    J.segment(dim, dim) = (2 / e_sqnorm * e.dot(e2p) * e - e - e2p) / e_sqnorm;
    J.tail(dim) = (e2p - 2 / e_sqnorm * e.dot(e2p) * e) / e_sqnorm;
    return J;
}

// ============================================================================
// Edge - Edge

Vector2<double> edge_edge_closest_point(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    const Eigen::Vector3d eb_to_ea = ea0 - eb0;
    const Eigen::Vector3d ea = ea1 - ea0;
    const Eigen::Vector3d eb = eb1 - eb0;

    Eigen::Matrix<double, 2, 2> coefMtr;
    coefMtr(0, 0) = ea.squaredNorm();
    coefMtr(0, 1) = coefMtr(1, 0) = -eb.dot(ea);
    coefMtr(1, 1) = eb.squaredNorm();

    Vector2<double> rhs;
    rhs[0] = -eb_to_ea.dot(ea);
    rhs[1] = eb_to_ea.dot(eb);

    const Eigen::Vector2d x = coefMtr.ldlt().solve(rhs);
    assert((coefMtr * x - rhs).norm() < 1e-10);
    return x;
}

Eigen::Matrix<double, 2, 12> edge_edge_closest_point_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    Eigen::Matrix<double, 2, 12> J;

    autogen::edge_edge_closest_point_jacobian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], J.data());

    return J;
}

// ============================================================================
// Point - Triangle

Vector2<double> point_triangle_closest_point(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2)
{
    Eigen::Matrix<double, 2, 3> basis;
    basis.row(0) = RowVector3<double>(t1 - t0); // edge 0
    basis.row(1) = RowVector3<double>(t2 - t0); // edge 1
    const Eigen::Matrix2d A = basis * basis.transpose();
    const Eigen::Vector2d b = basis * (p - t0);
    const Eigen::Vector2d x = A.ldlt().solve(b);
    assert((A * x - b).norm() < 1e-10);
    return x;
}

Eigen::Matrix<double, 2, 12> point_triangle_closest_point_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2)
{
    Eigen::Matrix<double, 2, 12> J;
    autogen::point_triangle_closest_point_jacobian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], J.data());

    return J;
}

// ============================================================================

namespace autogen {
    // J is (2×12) flattened in column-major order
    void edge_edge_closest_point_jacobian(
        double ea0_x,
        double ea0_y,
        double ea0_z,
        double ea1_x,
        double ea1_y,
        double ea1_z,
        double eb0_x,
        double eb0_y,
        double eb0_z,
        double eb1_x,
        double eb1_y,
        double eb1_z,
        double J[24])
    {
        const auto t0 = ea0_y * eb0_y;
        const auto t1 = ea0_x * eb0_x;
        const auto t2 = 2 * t1;
        const auto t3 = ea0_y * eb1_y;
        const auto t4 = ea0_x * eb1_x;
        const auto t5 = 2 * t4;
        const auto t6 = ea0_z * eb0_z;
        const auto t7 = ea0_z * eb1_z;
        const auto t8 = eb0_y * eb1_y;
        const auto t9 = 4 * t8;
        const auto t10 = ea0_x * ea1_x;
        const auto t11 = eb0_z * eb1_z;
        const auto t12 = 4 * t11;
        const auto t13 = ea1_y * eb0_y;
        const auto t14 = ea1_y * eb1_y;
        const auto t15 = ea1_z * eb0_z;
        const auto t16 = ea1_z * eb1_z;
        const auto t17 = 2 * t0;
        const auto t18 = 2 * t3;
        const auto t19 = ea1_x * eb0_x;
        const auto t20 = ea1_x * eb1_x;
        const auto t21 = ea0_y * ea1_y;
        const auto t22 = eb0_x * eb1_x;
        const auto t23 = 4 * t22;
        const auto t24 = 2 * t6;
        const auto t25 = 2 * t7;
        const auto t26 = ea0_z * ea1_z;
        const auto t27 = 2 * t19;
        const auto t28 = 2 * t20;
        const auto t29 = 2 * t13;
        const auto t30 = 2 * t14;
        const auto t31 = ea0_x * ea0_x;
        const auto t32 = eb0_y * eb0_y;
        const auto t33 = eb0_z * eb0_z;
        const auto t34 = eb1_y * eb1_y;
        const auto t35 = eb1_z * eb1_z;
        const auto t36 = ea0_y * ea0_y;
        const auto t37 = eb0_x * eb0_x;
        const auto t38 = eb1_x * eb1_x;
        const auto t39 = ea0_z * ea0_z;
        const auto t40 = ea1_x * ea1_x;
        const auto t41 = ea1_y * ea1_y;
        const auto t42 = ea1_z * ea1_z;
        const auto t43 = 2 * ea0_x;
        const auto t44 = ea1_x * t32;
        const auto t45 = ea1_x * t33;
        const auto t46 = ea1_x * t34;
        const auto t47 = ea1_x * t35;
        const auto t48 = 2 * eb0_y;
        const auto t49 = eb1_y * t31;
        const auto t50 = 2 * eb0_z;
        const auto t51 = eb1_z * t31;
        const auto t52 = 2 * ea0_y;
        const auto t53 = ea1_y * t37;
        const auto t54 = ea1_y * t33;
        const auto t55 = ea1_y * t38;
        const auto t56 = ea1_y * t35;
        const auto t57 = 2 * eb0_x;
        const auto t58 = eb1_x * t36;
        const auto t59 = eb1_z * t36;
        const auto t60 = 2 * ea0_z;
        const auto t61 = ea1_z * t37;
        const auto t62 = ea1_z * t32;
        const auto t63 = ea1_z * t38;
        const auto t64 = ea1_z * t34;
        const auto t65 = eb1_x * t39;
        const auto t66 = eb1_y * t39;
        const auto t67 = eb1_y * t40;
        const auto t68 = eb1_z * t40;
        const auto t69 = eb1_x * t41;
        const auto t70 = eb1_z * t41;
        const auto t71 = eb1_x * t42;
        const auto t72 = eb1_y * t42;
        const auto t73 = 1.0
            / (-t0 * t2 + t0 * t5 + t10 * t12 + t10 * t9 + t12 * t21 + t13 * t2
               + t13 * t24 - t13 * t25 - t13 * t27 + t13 * t28 - t13 * t5
               - t14 * t2 - t14 * t24 + t14 * t25 + t14 * t27 - t14 * t28
               + t14 * t5 + t15 * t17 - t15 * t18 + t15 * t2 - t15 * t27
               + t15 * t28 - t15 * t29 + t15 * t30 - t15 * t5 - t16 * t17
               + t16 * t18 - t16 * t2 + t16 * t27 - t16 * t28 + t16 * t29
               - t16 * t30 + t16 * t5 + t17 * t19 - t17 * t20 - t17 * t6
               + t17 * t7 - t18 * t19 + t18 * t20 + t18 * t6 - t18 * t7
               + t19 * t24 - t19 * t25 + t2 * t3 - t2 * t6 + t2 * t7 - t20 * t24
               + t20 * t25 + t21 * t23 + t23 * t26 + t26 * t9 - t3 * t5
               + t31 * t32 + t31 * t33 + t31 * t34 + t31 * t35 + t32 * t39
               + t32 * t40 + t32 * t42 + t33 * t36 + t33 * t40 + t33 * t41
               + t34 * t39 + t34 * t40 + t34 * t42 + t35 * t36 + t35 * t40
               + t35 * t41 + t36 * t37 + t36 * t38 + t37 * t39 + t37 * t41
               + t37 * t42 + t38 * t39 + t38 * t41 + t38 * t42 - t43 * t44
               - t43 * t45 - t43 * t46 - t43 * t47 - t48 * t49 - t48 * t66
               - t48 * t67 - t48 * t72 + t5 * t6 - t5 * t7 - t50 * t51
               - t50 * t59 - t50 * t68 - t50 * t70 - t52 * t53 - t52 * t54
               - t52 * t55 - t52 * t56 - t57 * t58 - t57 * t65 - t57 * t69
               - t57 * t71 - t60 * t61 - t60 * t62 - t60 * t63 - t60 * t64);
        const auto t74 = ea1_x + eb0_x - t43;
        const auto t75 = eb1_x * t57;
        const auto t76 = eb1_y * t48;
        const auto t77 = eb1_z * t50;
        const auto t78 = t32 + t33 + t34 + t35 + t37 + t38 - t75 - t76 - t77;
        const auto t79 = eb0_x - eb1_x;
        const auto t80 = ea0_x - eb0_x;
        const auto t81 = ea0_y - eb0_y;
        const auto t82 = eb0_y - eb1_y;
        const auto t83 = ea0_z - eb0_z;
        const auto t84 = eb0_z - eb1_z;
        const auto t85 = t79 * t80 + t81 * t82 + t83 * t84;
        const auto t86 = t79 * t85;
        const auto t87 =
            t0 + t1 - t13 + t14 - t15 + t16 - t19 + t20 - t3 - t4 + t6 - t7;
        const auto t88 = t80 * (-ea0_x + ea1_x) + t81 * (-ea0_y + ea1_y)
            + t83 * (-ea0_z + ea1_z);
        const auto t89 = 2 * t73;
        const auto t90 = t89
            * (ea0_x * t32 + ea0_x * t33 + ea0_x * t34 + ea0_x * t35
               + ea1_x * t76 + ea1_x * t77 - eb0_x * t0 + eb0_x * t13
               - eb0_x * t14 + eb0_x * t15 - eb0_x * t16 + eb0_x * t3
               - eb0_x * t6 + eb0_x * t7 + eb1_x * t0 - eb1_x * t13
               + eb1_x * t14 - eb1_x * t15 + eb1_x * t16 - eb1_x * t3
               + eb1_x * t6 - eb1_x * t7 - t11 * t43 - t43 * t8 - t44 - t45
               - t46 - t47);
        const auto t91 = t88 * t90;
        const auto t92 = t78 * t91;
        const auto t93 = t85 * t90;
        const auto t94 = t87 * t93;
        const auto t95 = ea1_x * t43;
        const auto t96 = ea1_y * t52;
        const auto t97 = ea1_z * t60;
        const auto t98 = t31 + t36 + t39 + t40 + t41 + t42 - t95 - t96 - t97;
        const auto t99 = ea0_x - ea1_x;
        const auto t100 = t85 * t99;
        const auto t101 = 2 * t100;
        const auto t102 = t79 * t88;
        const auto t103 = t93 * t98;
        const auto t104 = t87 * t91;
        const auto t105 = ea1_y + eb0_y - t52;
        const auto t106 = -ea0_y * t33 - ea0_y * t35 - ea0_y * t37 - ea0_y * t38
            - ea1_y * t75 - ea1_y * t77 + eb0_y * t1 - eb0_y * t15 + eb0_y * t16
            - eb0_y * t19 + eb0_y * t20 - eb0_y * t4 + eb0_y * t6 - eb0_y * t7
            - eb1_y * t1 + eb1_y * t15 - eb1_y * t16 + eb1_y * t19 - eb1_y * t20
            + eb1_y * t4 - eb1_y * t6 + eb1_y * t7 + t11 * t52 + t22 * t52 + t53
            + t54 + t55 + t56;
        const auto t107 = t106 * t89;
        const auto t108 = t107 * t88;
        const auto t109 = t107 * t85;
        const auto t110 = t109 * t87 + t82 * t85;
        const auto t111 = t82 * t88;
        const auto t112 = ea0_y - ea1_y;
        const auto t113 = t112 * t85;
        const auto t114 = t109 * t98 + 2 * t113;
        const auto t115 = ea1_z + eb0_z - t60;
        const auto t116 = -ea0_z * t32 - ea0_z * t34 - ea0_z * t37 - ea0_z * t38
            - ea1_z * t75 - ea1_z * t76 + eb0_z * t0 + eb0_z * t1 - eb0_z * t13
            + eb0_z * t14 - eb0_z * t19 + eb0_z * t20 - eb0_z * t3 - eb0_z * t4
            - eb1_z * t0 - eb1_z * t1 + eb1_z * t13 - eb1_z * t14 + eb1_z * t19
            - eb1_z * t20 + eb1_z * t3 + eb1_z * t4 + t22 * t60 + t60 * t8 + t61
            + t62 + t63 + t64;
        const auto t117 = t116 * t89;
        const auto t118 = t117 * t88;
        const auto t119 = t117 * t85;
        const auto t120 = t119 * t87 + t84 * t85;
        const auto t121 = t84 * t88;
        const auto t122 = ea0_z - ea1_z;
        const auto t123 = t122 * t85;
        const auto t124 = t119 * t98 + 2 * t123;
        const auto t125 = t80 * t87;
        const auto t126 = t112 * t81 + t122 * t83 + t80 * t99;
        const auto t127 = 2 * t126;
        const auto t128 = t127 * t73;
        const auto t129 = t106 * t128;
        const auto t130 = t81 * t87;
        const auto t131 = t126 * t82;
        const auto t132 = t116 * t128;
        const auto t133 = t83 * t87;
        const auto t134 = t126 * t84;
        const auto t135 = ea0_x + eb1_x - t57;
        const auto t136 = ea0_x * t0 - ea0_x * t13 + ea0_x * t14 - ea0_x * t15
            + ea0_x * t16 - ea0_x * t3 + ea0_x * t6 - ea0_x * t7 - ea1_x * t0
            + ea1_x * t13 - ea1_x * t14 + ea1_x * t15 - ea1_x * t16 + ea1_x * t3
            - ea1_x * t6 + ea1_x * t7 - eb0_x * t36 - eb0_x * t39 - eb0_x * t41
            - eb0_x * t42 + eb0_x * t96 + eb0_x * t97 - eb1_x * t96
            - eb1_x * t97 + t58 + t65 + t69 + t71;
        const auto t137 = t136 * t89;
        const auto t138 = t137 * t88;
        const auto t139 = t137 * t85;
        const auto t140 = t100 + t139 * t87;
        const auto t141 = t139 * t98;
        const auto t142 = ea0_y + eb1_y - t48;
        const auto t143 = -ea0_y * t1 + ea0_y * t15 - ea0_y * t16 + ea0_y * t19
            - ea0_y * t20 + ea0_y * t4 - ea0_y * t6 + ea0_y * t7 + ea1_y * t1
            - ea1_y * t15 + ea1_y * t16 - ea1_y * t19 + ea1_y * t20 - ea1_y * t4
            + ea1_y * t6 - ea1_y * t7 + eb0_y * t31 + eb0_y * t39 + eb0_y * t40
            + eb0_y * t42 - eb0_y * t95 - eb0_y * t97 + eb1_y * t95
            + eb1_y * t97 - t49 - t66 - t67 - t72;
        const auto t144 = t143 * t89;
        const auto t145 = t144 * t88;
        const auto t146 = t144 * t85;
        const auto t147 = t146 * t87;
        const auto t148 = t146 * t98;
        const auto t149 = ea0_z + eb1_z - t50;
        const auto t150 = -ea0_z * t0 - ea0_z * t1 + ea0_z * t13 - ea0_z * t14
            + ea0_z * t19 - ea0_z * t20 + ea0_z * t3 + ea0_z * t4 + ea1_z * t0
            + ea1_z * t1 - ea1_z * t13 + ea1_z * t14 - ea1_z * t19 + ea1_z * t20
            - ea1_z * t3 - ea1_z * t4 + eb0_z * t31 + eb0_z * t36 + eb0_z * t40
            + eb0_z * t41 - eb0_z * t95 - eb0_z * t96 + eb1_z * t95
            + eb1_z * t96 - t51 - t59 - t68 - t70;
        const auto t151 = t150 * t89;
        const auto t152 = t151 * t88;
        const auto t153 = t151 * t85;
        const auto t154 = t153 * t87;
        const auto t155 = t153 * t98;
        const auto t156 = t128 * t136;
        const auto t157 = t128 * t143;
        const auto t158 = t128 * t150;
        J[0] = t73 * (-t74 * t78 - t79 * t87 - t86 + t92 + t94);
        J[1] = t73 * (-t101 - t102 + t103 + t104 - t74 * t87 - t79 * t98);
        J[2] = -t73 * (t105 * t78 + t108 * t78 + t110 + t82 * t87);
        J[3] = -t73 * (t105 * t87 + t108 * t87 + t111 + t114 + t82 * t98);
        J[4] = -t73 * (t115 * t78 + t118 * t78 + t120 + t84 * t87);
        J[5] = -t73 * (t115 * t87 + t118 * t87 + t121 + t124 + t84 * t98);
        J[6] = t73 * (-t78 * t80 + t86 - t92 - t94);
        J[7] = t73 * (t101 + t102 - t103 - t104 - t125);
        J[8] = t73 * (t110 - t129 * t78 - t78 * t81);
        J[9] = t73 * (t114 - t129 * t87 - t130 - t131);
        J[10] = t73 * (t120 - t132 * t78 - t78 * t83);
        J[11] = t73 * (t124 - t132 * t87 - t133 - t134);
        J[12] = -t73 * (2 * t102 + t135 * t87 + t138 * t78 + t140 + t78 * t99);
        J[13] = -t73 * (t135 * t98 + t138 * t87 + t141 + t87 * t99 + t88 * t99);
        J[14] = t73
            * (-2 * t111 - t112 * t78 - t113 - t142 * t87 + t145 * t78 + t147);
        J[15] =
            t73 * (-t112 * t87 - t112 * t88 - t142 * t98 + t145 * t87 + t148);
        J[16] = t73
            * (-2 * t121 - t122 * t78 - t123 - t149 * t87 + t152 * t78 + t154);
        J[17] =
            t73 * (-t122 * t87 - t122 * t88 - t149 * t98 + t152 * t87 + t155);
        J[18] = t73 * (t125 - t127 * t79 + t140 - t156 * t78);
        J[19] = t73 * (-t126 * t99 + t141 - t156 * t87 + t80 * t98);
        J[20] = t73 * (t113 + t130 - 2 * t131 - t147 + t157 * t78);
        J[21] = t73 * (-t112 * t126 - t148 + t157 * t87 + t81 * t98);
        J[22] = t73 * (t123 + t133 - 2 * t134 - t154 + t158 * t78);
        J[23] = t73 * (-t122 * t126 - t155 + t158 * t87 + t83 * t98);
    }

    // J is (2×12) flattened in column-major order
    void point_triangle_closest_point_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double J[24])
    {
        const auto t0 = t0_x - t1_x;
        const auto t1 = t2_x * t2_x;
        const auto t2 = t2_y * t2_y;
        const auto t3 = t2_z * t2_z;
        const auto t4 = t0_x * t2_x;
        const auto t5 = 2 * t4;
        const auto t6 = t0_y * t2_y;
        const auto t7 = 2 * t6;
        const auto t8 = t0_z * t2_z;
        const auto t9 = 2 * t8;
        const auto t10 = t0_x * t0_x;
        const auto t11 = t0_y * t0_y;
        const auto t12 = t0_z * t0_z;
        const auto t13 = t10 + t11 + t12;
        const auto t14 = t1 + t13 + t2 + t3 - t5 - t7 - t9;
        const auto t15 = t0_x - t2_x;
        const auto t16 = t1_x * t2_x;
        const auto t17 = t1_y * t2_y;
        const auto t18 = t1_z * t2_z;
        const auto t19 = t0_x * t1_x;
        const auto t20 = t0_y * t1_y;
        const auto t21 = t0_z * t1_z;
        const auto t22 = t13 + t16 + t17 + t18 - t19 - t20 - t21 - t4 - t6 - t8;
        const auto t23 = 2 * t19;
        const auto t24 = 2 * t20;
        const auto t25 = 2 * t21;
        const auto t26 = 2 * t16;
        const auto t27 = 2 * t17;
        const auto t28 = t1_y * t1_y;
        const auto t29 = t1_z * t1_z;
        const auto t30 = t1_x * t1_x;
        const auto t31 = 2 * t18;
        const auto t32 = 1.0
            / (t1 * t11 + t1 * t12 - t1 * t24 - t1 * t25 + t1 * t28 + t1 * t29
               + t10 * t2 - t10 * t27 + t10 * t28 + t10 * t29 + t10 * t3
               - t10 * t31 - t11 * t26 + t11 * t29 + t11 * t3 + t11 * t30
               - t11 * t31 + t12 * t2 - t12 * t26 - t12 * t27 + t12 * t28
               + t12 * t30 + t16 * t24 + t16 * t25 + t16 * t7 + t16 * t9
               + t17 * t23 + t17 * t25 - t17 * t26 + t17 * t5 + t17 * t9
               + t18 * t23 + t18 * t24 - t18 * t26 - t18 * t27 + t18 * t5
               + t18 * t7 + t19 * t7 + t19 * t9 - t2 * t23 - t2 * t25 + t2 * t29
               + t2 * t30 - t20 * t23 + t20 * t5 + t20 * t9 - t21 * t23
               - t21 * t24 + t21 * t5 + t21 * t7 - t23 * t3 - t24 * t3
               + t28 * t3 - t28 * t5 - t28 * t9 - t29 * t5 - t29 * t7 + t3 * t30
               - t30 * t7 - t30 * t9 - t5 * t6 - t5 * t8 - t7 * t8);
        const auto t33 = t13 - t23 - t24 - t25 + t28 + t29 + t30;
        const auto t34 = t0_y - t1_y;
        const auto t35 = t0_y - t2_y;
        const auto t36 = t0_z - t1_z;
        const auto t37 = t0_z - t2_z;
        const auto t38 = 2 * t0_x;
        const auto t39 = -t38;
        const auto t40 = p_x + t39;
        const auto t41 = t1_x + t40;
        const auto t42 = p_x - t0_x;
        const auto t43 = p_y - t0_y;
        const auto t44 = p_z - t0_z;
        const auto t45 = t0 * t42 + t34 * t43 + t36 * t44;
        const auto t46 = t15 * t45;
        const auto t47 = 2 * t46;
        const auto t48 = t1_x + t2_x + t39;
        const auto t49 = t15 * t42 + t35 * t43 + t37 * t44;
        const auto t50 = t2_x + t40;
        const auto t51 = t0_x * t28;
        const auto t52 = t0_x * t29;
        const auto t53 = t0_x * t2;
        const auto t54 = t0_x * t3;
        const auto t55 = t1_x * t2;
        const auto t56 = t1_x * t3;
        const auto t57 = t28 * t2_x;
        const auto t58 = t29 * t2_x;
        const auto t59 = t17 * t1_x;
        const auto t60 = t18 * t1_x;
        const auto t61 = t17 * t2_x;
        const auto t62 = t18 * t2_x;
        const auto t63 = t1_x * t20;
        const auto t64 = t2_x * t6;
        const auto t65 = t1_x * t21;
        const auto t66 = t2_x * t8;
        const auto t67 = t20 * t2_x + t21 * t2_x;
        const auto t68 = t1_x * t6 + t1_x * t8;
        const auto t69 = 2 * t32;
        const auto t70 = t69
            * (-t17 * t38 - t18 * t38 + t51 + t52 + t53 + t54 - t55 - t56 - t57
               - t58 + t59 + t60 + t61 + t62 - t63 - t64 - t65 - t66 + t67
               + t68);
        const auto t71 = t14 * t45;
        const auto t72 = t22 * t70;
        const auto t73 = t0 * t49;
        const auto t74 = 2 * t73;
        const auto t75 = t33 * t49;
        const auto t76 = 2 * t0_y;
        const auto t77 = -t76;
        const auto t78 = p_y + t77;
        const auto t79 = t1_y + t78;
        const auto t80 = t35 * t45;
        const auto t81 = 2 * t80;
        const auto t82 = t1_y + t2_y + t77;
        const auto t83 = t2_y + t78;
        const auto t84 = t0_y * t30;
        const auto t85 = t0_y * t29;
        const auto t86 = t0_y * t1;
        const auto t87 = t0_y * t3;
        const auto t88 = t1 * t1_y;
        const auto t89 = t1_y * t3;
        const auto t90 = t2_y * t30;
        const auto t91 = t29 * t2_y;
        const auto t92 = t16 * t1_y;
        const auto t93 = t16 * t2_y;
        const auto t94 = t18 * t1_y;
        const auto t95 = t18 * t2_y;
        const auto t96 = t19 * t1_y;
        const auto t97 = t2_y * t4;
        const auto t98 = t1_y * t21;
        const auto t99 = t2_y * t8;
        const auto t100 = t19 * t2_y + t21 * t2_y;
        const auto t101 = t1_y * t4 + t1_y * t8;
        const auto t102 = t69
            * (t100 + t101 - t16 * t76 - t18 * t76 + t84 + t85 + t86 + t87 - t88
               - t89 - t90 - t91 + t92 + t93 + t94 + t95 - t96 - t97 - t98
               - t99);
        const auto t103 = t102 * t22;
        const auto t104 = t34 * t49;
        const auto t105 = 2 * t104;
        const auto t106 = 2 * t0_z;
        const auto t107 = -t106;
        const auto t108 = p_z + t107;
        const auto t109 = t108 + t1_z;
        const auto t110 = t37 * t45;
        const auto t111 = 2 * t110;
        const auto t112 = t107 + t1_z + t2_z;
        const auto t113 = t108 + t2_z;
        const auto t114 = t0_z * t30;
        const auto t115 = t0_z * t28;
        const auto t116 = t0_z * t1;
        const auto t117 = t0_z * t2;
        const auto t118 = t1 * t1_z;
        const auto t119 = t1_z * t2;
        const auto t120 = t2_z * t30;
        const auto t121 = t28 * t2_z;
        const auto t122 = t16 * t1_z;
        const auto t123 = t16 * t2_z;
        const auto t124 = t17 * t1_z;
        const auto t125 = t17 * t2_z;
        const auto t126 = t19 * t1_z;
        const auto t127 = t2_z * t4;
        const auto t128 = t1_z * t20;
        const auto t129 = t2_z * t6;
        const auto t130 = t19 * t2_z + t20 * t2_z;
        const auto t131 = t1_z * t4 + t1_z * t6;
        const auto t132 = t69
            * (-t106 * t16 - t106 * t17 + t114 + t115 + t116 + t117 - t118
               - t119 - t120 - t121 + t122 + t123 + t124 + t125 - t126 - t127
               - t128 - t129 + t130 + t131);
        const auto t133 = t132 * t22;
        const auto t134 = t36 * t49;
        const auto t135 = 2 * t134;
        const auto t136 = t11 * t1_x;
        const auto t137 = t12 * t1_x;
        const auto t138 = t11 * t2_x;
        const auto t139 = t12 * t2_x;
        const auto t140 = t0_x * t6;
        const auto t141 = t0_x * t8;
        const auto t142 = t0_x * t20;
        const auto t143 = t0_x * t21;
        const auto t144 = t0_x * t17 + t0_x * t18;
        const auto t145 = t69
            * (t136 + t137 - t138 - t139 + t140 + t141 - t142 - t143 + t144
               - t1_x * t7 - t1_x * t9 - t53 - t54 + t55 + t56 - t61 - t62 + t64
               + t66 + t67);
        const auto t146 = t145 * t22;
        const auto t147 = -t22 * t42;
        const auto t148 = t10 * t1_y;
        const auto t149 = t12 * t1_y;
        const auto t150 = t10 * t2_y;
        const auto t151 = t12 * t2_y;
        const auto t152 = t0_y * t4;
        const auto t153 = t0_y * t8;
        const auto t154 = t0_y * t19;
        const auto t155 = t0_y * t21;
        const auto t156 = t0_y * t16 + t0_y * t18;
        const auto t157 = t69
            * (t100 + t148 + t149 - t150 - t151 + t152 + t153 - t154 - t155
               + t156 - t1_y * t5 - t1_y * t9 - t86 - t87 + t88 + t89 - t93
               - t95 + t97 + t99);
        const auto t158 = t157 * t22;
        const auto t159 = -t22 * t43;
        const auto t160 = t10 * t1_z;
        const auto t161 = t11 * t1_z;
        const auto t162 = t10 * t2_z;
        const auto t163 = t11 * t2_z;
        const auto t164 = t0_z * t4;
        const auto t165 = t0_z * t6;
        const auto t166 = t0_z * t19;
        const auto t167 = t0_z * t20;
        const auto t168 = t0_z * t16 + t0_z * t17;
        const auto t169 = t69
            * (-t116 - t117 + t118 + t119 - t123 - t125 + t127 + t129 + t130
               + t160 + t161 - t162 - t163 + t164 + t165 - t166 - t167 + t168
               - t1_z * t5 - t1_z * t7);
        const auto t170 = t169 * t22;
        const auto t171 = -t22 * t44;
        const auto t172 = t69
            * (-t136 - t137 + t138 + t139 - t140 - t141 + t142 + t143 + t144
               - t24 * t2_x - t25 * t2_x - t51 - t52 + t57 + t58 - t59 - t60
               + t63 + t65 + t68);
        const auto t173 = t172 * t22;
        const auto t174 = t69
            * (t101 - t148 - t149 + t150 + t151 - t152 - t153 + t154 + t155
               + t156 - t23 * t2_y - t25 * t2_y - t84 - t85 + t90 + t91 - t92
               - t94 + t96 + t98);
        const auto t175 = t174 * t22;
        const auto t176 = t69
            * (-t114 - t115 + t120 + t121 - t122 - t124 + t126 + t128 + t131
               - t160 - t161 + t162 + t163 - t164 - t165 + t166 + t167 + t168
               - t23 * t2_z - t24 * t2_z);
        const auto t177 = t176 * t22;
        J[0] = t32 * (-t0 * t14 + t15 * t22);
        J[1] = t32 * (t0 * t22 - t15 * t33);
        J[2] = t32 * (-t14 * t34 + t22 * t35);
        J[3] = t32 * (t22 * t34 - t33 * t35);
        J[4] = t32 * (-t14 * t36 + t22 * t37);
        J[5] = t32 * (t22 * t36 - t33 * t37);
        J[6] = t32
            * (-t14 * t41 + t22 * t50 - t47 - t48 * t49 - t49 * t72
               + t70 * t71);
        J[7] = t32
            * (t22 * t41 - t33 * t50 - t45 * t48 - t45 * t72 + t70 * t75 - t74);
        J[8] = t32
            * (t102 * t71 - t103 * t49 - t14 * t79 + t22 * t83 - t49 * t82
               - t81);
        J[9] = t32
            * (t102 * t75 - t103 * t45 - t105 + t22 * t79 - t33 * t83
               - t45 * t82);
        J[10] = t32
            * (-t109 * t14 - t111 - t112 * t49 + t113 * t22 + t132 * t71
               - t133 * t49);
        J[11] = t32
            * (t109 * t22 - t112 * t45 - t113 * t33 + t132 * t75 - t133 * t45
               - t135);
        J[12] = t32 * (t14 * t42 + t145 * t71 - t146 * t49 - t15 * t49);
        J[13] = t32 * (t145 * t75 - t146 * t45 + t147 - t46 + t74);
        J[14] = t32 * (t14 * t43 + t157 * t71 - t158 * t49 - t35 * t49);
        J[15] = t32 * (t105 + t157 * t75 - t158 * t45 + t159 - t80);
        J[16] = t32 * (t14 * t44 + t169 * t71 - t170 * t49 - t37 * t49);
        J[17] = t32 * (-t110 + t135 + t169 * t75 - t170 * t45 + t171);
        J[18] = t32 * (t147 + t172 * t71 - t173 * t49 + t47 - t73);
        J[19] = t32 * (-t0 * t45 + t172 * t75 - t173 * t45 + t33 * t42);
        J[20] = t32 * (-t104 + t159 + t174 * t71 - t175 * t49 + t81);
        J[21] = t32 * (t174 * t75 - t175 * t45 + t33 * t43 - t34 * t45);
        J[22] = t32 * (t111 - t134 + t171 + t176 * t71 - t177 * t49);
        J[23] = t32 * (t176 * t75 - t177 * t45 + t33 * t44 - t36 * t45);
    }

} // namespace autogen
} // namespace ipc
