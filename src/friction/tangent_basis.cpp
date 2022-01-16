#include <ipc/friction/tangent_basis.hpp>

namespace ipc {
namespace autogen {

    void point_point_tangent_basis_2D_jacobian(
        double p0_x, double p0_y, double p1_x, double p1_y, double J[8])
    {
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        t0 = p0_x - p1_x;
        t1 = p0_y - p1_y;
        t2 = std::pow(t0, 2);
        t3 = std::pow(t1, 2);
        t4 = t2 + t3;
        t5 = t0 * t1 / std::pow(t4, 3.0 / 2.0);
        t6 = -t5;
        t7 = std::pow(t4, -1.0 / 2.0);
        t8 = 1.0 / t4;
        t9 = t2 * t8;
        t10 = t3 * t8;
        J[0] = t6;
        J[1] = t7 * (t9 - 1);
        J[2] = t7 * (1 - t10);
        J[3] = t5;
        J[4] = t5;
        J[5] = t7 * (1 - t9);
        J[6] = t7 * (t10 - 1);
        J[7] = t6;
    }

    void point_point_tangent_basis_3D_jacobian(
        double p0_x,
        double p0_y,
        double p0_z,
        double p1_x,
        double p1_y,
        double p1_z,
        double J[36])
    {
        const auto t0 = p0_x - p1_x;
        const auto t1 = std::pow(t0, 2);
        const auto t2 = p0_z - p1_z;
        const auto t3 = std::pow(t2, 2);
        const auto t4 = p0_y - p1_y;
        const auto t5 = std::pow(t4, 2);
        const auto t6 = t1 + t3 < t3 + t5;
        const auto t7 = -p0_z + p1_z;
        const auto t8 = t6 ? 0 : t7;
        const auto t9 =
            ((t6) ? (0) : (t3)) + ((t6) ? (t3) : (0)) + ((t6) ? (t5) : (t1));
        const auto t10 = std::pow(t9, -3.0 / 2.0);
        const auto t11 = (1.0 / 2.0) * ((t6) ? (0) : (2 * t0));
        const auto t12 = t10 * t11;
        const auto t13 = t4 * t8;
        const auto t14 = t6 ? t2 : 0;
        const auto t15 = t0 * t14;
        const auto t16 = t13 - t15;
        const auto t17 = -p0_y + p1_y;
        const auto t18 = t1 + std::pow(t7, 2) < std::pow(t17, 2) + t3;
        const auto t19 = t18 ? 0 : t7;
        const auto t20 = t19 * t7;
        const auto t21 = -p0_x + p1_x;
        const auto t22 = t18 ? t17 : t0;
        const auto t23 = t21 * t22;
        const auto t24 = -t20 + t23;
        const auto t25 = t14 * t2;
        const auto t26 = t6 ? t17 : t0;
        const auto t27 = t26 * t4;
        const auto t28 = t25 - t27;
        const auto t29 = std::pow(t16, 2) + std::pow(t24, 2) + std::pow(t28, 2);
        const auto t30 = std::pow(t29, -1.0 / 2.0);
        const auto t31 = t6 ? 0 : 1;
        const auto t32 = t31 * t4;
        const auto t33 = 1.0 / t29;
        const auto t34 = t18 ? t2 : 0;
        const auto t35 = -t17 * t19 + t21 * t34;
        const auto t36 = t34 * t35;
        const auto t37 = t18 ? 0 : 1;
        const auto t38 = t17 * t22 - t34 * t7;
        const auto t39 = t17 * t38;
        const auto t40 = t20 - t23;
        const auto t41 = t33 * (-t36 + t37 * t39 + t40 * (t0 * t37 + t22));
        const auto t42 = -t13 + t15;
        const auto t43 = -t0 * t26 + t2 * t8;
        const auto t44 = -t25 + t27;
        const auto t45 = std::pow(t42, 2) + std::pow(t43, 2) + std::pow(t44, 2);
        const auto t46 = std::pow(t45, -1.0 / 2.0);
        const auto t47 = t14 * t16;
        const auto t48 = t0 * t31 + t26;
        const auto t49 = 1.0 / t45;
        const auto t50 = t43 * t49;
        const auto t51 = std::pow(t9, -1.0 / 2.0);
        const auto t52 = 1.0 / t9;
        const auto t53 = t26 * t52;
        const auto t54 = (1.0 / 2.0) * ((t6) ? (2 * t4) : (0));
        const auto t55 = t10 * t54;
        const auto t56 = t19 * t35;
        const auto t57 = t18 ? -1 : 0;
        const auto t58 = t0 * t40;
        const auto t59 = -t22;
        const auto t60 = t33 * (t38 * (t17 * t57 + t59) + t56 + t57 * t58);
        const auto t61 = t6 ? -1 : 0;
        const auto t62 = t26 + t4 * t61;
        const auto t63 = t0 * t61;
        const auto t64 = t16 * t8;
        const auto t65 = t6 ? 0 : -1;
        const auto t66 = (1.0 / 2.0) * t8;
        const auto t67 = 2 * t2;
        const auto t68 = t52 * (((t18) ? (0) : (t67)) + ((t18) ? (t67) : (0)));
        const auto t69 = t18 ? 1 : 0;
        const auto t70 = t18 ? 0 : -1;
        const auto t71 = t28 * t33;
        const auto t72 = t6 ? 1 : 0;
        const auto t73 = t14 + t2 * t72;
        const auto t74 = (1.0 / 2.0) * t14;
        const auto t75 = t2 * t65;
        const auto t76 = t4 * t65;
        const auto t77 = t0 * t72;
        const auto t78 = t76 - t77;
        const auto t79 = t16 * t78 + t24 * (t75 + t8) + t28 * t73;
        const auto t80 = (1.0 / 2.0) * t10 * t26;
        const auto t81 = t42 * t49;
        const auto t82 = t6 ? 0 : 2 * t21;
        const auto t83 = t10 * t82;
        const auto t84 = t36 + t39 * t70 + t40 * (t0 * t70 + t59);
        const auto t85 = t0 * t65 - t26;
        const auto t86 = (1.0 / 2.0) * t53;
        const auto t87 = t16 * t33;
        const auto t88 = t6 ? 2 * t17 : 0;
        const auto t89 = t10 * t88;
        const auto t90 = t26 - t4 * t72;
        const auto t91 = -t24 * t77 + t28 * t90 - t64;
        const auto t92 = 2 * t7;
        const auto t93 = t52 * (((t18) ? (0) : (t92)) + ((t18) ? (t92) : (0)));
        const auto t94 = -t14 + t2 * t61;
        const auto t95 = t32 - t63;
        const auto t96 = -t2 * t31 + t8;
        const auto t97 = t16 * t95 - t24 * t96 + t28 * t94;
        J[0] = -t12 * t8;
        J[1] = -t30 * (t28 * t41 + t32);
        J[2] = -t12 * t14;
        J[3] = t46 * (t48 - t50 * (t24 * t48 + t28 * t32 + t47));
        J[4] = t51 * (-t11 * t53 + t31);
        J[5] = -t30 * (t14 + t16 * t41);
        J[6] = -t55 * t8;
        J[7] = -t30 * (t28 * t60 + t62);
        J[8] = -t14 * t55;
        J[9] = t46 * (t50 * (-t24 * t63 - t28 * t62 + t64) + t63);
        J[10] = t51 * (-t53 * t54 + t61);
        J[11] = t30 * (-t16 * t60 + t8);
        J[12] = t51 * (t65 - t66 * t68);
        J[13] = t30
            * (-t71
                   * (t35 * (t21 * t69 + t4 * t70) + t38 * (t2 * t69 + t34)
                      + t40 * (-t19 + t7 * t70))
               + t73);
        J[14] = t51 * (-t68 * t74 + t72);
        J[15] = t46 * (t50 * t79 - t75 - t8);
        J[16] = -t80 * (((t6) ? (0) : (t67)) + ((t6) ? (t67) : (0)));
        J[17] = t46 * (t78 + t79 * t81);
        J[18] = -t66 * t83;
        J[19] = -t30 * (t71 * t84 + t76);
        J[20] = -t74 * t83;
        J[21] = t46 * (t50 * (-t24 * t85 - t28 * t76 + t47) + t85);
        J[22] = t51 * (t65 - t82 * t86);
        J[23] = t30 * (t14 - t84 * t87);
        J[24] = -t66 * t89;
        J[25] = t46 * (t44 * t49 * t91 + t90);
        J[26] = -t74 * t89;
        J[27] = t46 * (t50 * t91 + t77);
        J[28] = t51 * (t72 - t86 * t88);
        J[29] = -t30 * (t8 + t87 * (t38 * (t17 * t69 + t22) - t56 + t58 * t69));
        J[30] = t51 * (t31 - t66 * t93);
        J[31] = t30
            * (-t71
                   * (t35 * (t21 * t57 + t37 * t4) + t38 * (t2 * t57 - t34)
                      + t40 * (t19 + t37 * t7))
               + t94);
        J[32] = t51 * (t61 - t74 * t93);
        J[33] = t46 * (t50 * t97 + t96);
        J[34] = -t80 * (((t6) ? (0) : (t92)) + ((t6) ? (t92) : (0)));
        J[35] = t46 * (t81 * t97 + t95);
    }

    void point_edge_tangent_basis_2D_jacobian(
        double p_x,
        double p_y,
        double e0_x,
        double e0_y,
        double e1_x,
        double e1_y,
        double J[12])
    {
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        t0 = e0_x - e1_x;
        t1 = std::pow(t0, 2);
        t2 = e0_y - e1_y;
        t3 = std::pow(t2, 2);
        t4 = t1 + t3;
        t5 = std::pow(t4, -1.0 / 2.0);
        t6 = 1.0 / t4;
        t7 = t1 * t6;
        t8 = t0 * t2 / std::pow(t4, 3.0 / 2.0);
        t9 = t3 * t6;
        t10 = -t8;
        J[0] = 0;
        J[1] = 0;
        J[2] = 0;
        J[3] = 0;
        J[4] = t5 * (t7 - 1);
        J[5] = t8;
        J[6] = t8;
        J[7] = t5 * (t9 - 1);
        J[8] = t5 * (1 - t7);
        J[9] = t10;
        J[10] = t10;
        J[11] = t5 * (1 - t9);
    }

    void point_edge_tangent_basis_3D_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double J[54])
    {
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14,
            t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,
            t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40,
            t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53,
            t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66,
            t67, t68, t69, t70, t71, t72, t73, t74, t75, t76;
        t0 = -e1_y;
        t1 = e0_y + t0;
        t2 = -e1_x;
        t3 = e0_x + t2;
        t4 = -p_y;
        t5 = e0_y + t4;
        t6 = -p_x;
        t7 = e0_x + t6;
        t8 = -t1 * t7 + t3 * t5;
        t9 = -e1_z;
        t10 = e0_z + t9;
        t11 = -e0_x;
        t12 = e1_x + t11;
        t13 = -e0_z;
        t14 = p_z + t13;
        t15 = t12 * t14;
        t16 = p_x + t11;
        t17 = e1_z + t13;
        t18 = t16 * t17;
        t19 = t15 - t18;
        t20 = t1 * t8 + t10 * t19;
        t21 = -p_z;
        t22 = e0_z + t21;
        t23 = t1 * t22 - t10 * t5;
        t24 = -t10 * t7 + t22 * t3;
        t25 = std::pow(t23, 2) + std::pow(t8, 2);
        t26 = std::pow(t24, 2) + t25;
        t27 = std::pow(t26, -3.0 / 2.0);
        t28 = t23 * t27;
        t29 = std::pow(t19, 2) + t25;
        t30 = std::pow(t29, -1.0 / 2.0);
        t31 = -e0_y;
        t32 = p_y + t31;
        t33 = e1_y + t31;
        t34 = t12 * t32 - t16 * t33;
        t35 = -t15 + t18;
        t36 = 1.0 / t29;
        t37 = t19 * t36;
        t38 = std::pow(t26, -1.0 / 2.0);
        t39 = 1.0 / t26;
        t40 = t39 * t8;
        t41 = t10 * t23 - t3 * t8;
        t42 = t23 * t39;
        t43 = t24 * t27;
        t44 = t14 * t33 - t17 * t32;
        t45 = t36 * t8;
        t46 = t23 * t36;
        t47 = t1 * t23 + t19 * t3;
        t48 = t24 * t39;
        t49 = t27 * t8;
        t50 = std::pow(t3, 2);
        t51 = std::pow(t1, 2);
        t52 = std::pow(t10, 2);
        t53 = t50 + t51 + t52;
        t54 = std::pow(t53, -1.0 / 2.0);
        t55 = 1.0 / t53;
        t56 = t50 * t55;
        t57 = e1_y + t4;
        t58 = e1_z + t21;
        t59 = t19 * t58 + t57 * t8;
        t60 = std::pow(t53, -3.0 / 2.0);
        t61 = t3 * t60;
        t62 = t1 * t61;
        t63 = p_z + t9;
        t64 = t10 * t61;
        t65 = e1_x + t6;
        t66 = t23 * t58 - t65 * t8;
        t67 = t51 * t55;
        t68 = t1 * t10 * t60;
        t69 = t19 * t65 + t23 * t57;
        t70 = t52 * t55;
        t71 = t19 * t22 + t5 * t8;
        t72 = -t62;
        t73 = -t64;
        t74 = -t22 * t23 + t7 * t8;
        t75 = -t68;
        t76 = t19 * t7 + t23 * t5;
        J[0] = 0;
        J[1] = -t20 * t28;
        J[2] = 0;
        J[3] = t30 * (t17 + t37 * (t1 * t34 + t17 * t35));
        J[4] = 0;
        J[5] = t38 * (t1 - t20 * t40);
        J[6] = 0;
        J[7] = t38 * (t10 - t41 * t42);
        J[8] = 0;
        J[9] = t41 * t43;
        J[10] = 0;
        J[11] = -t30 * (t3 + t45 * (t10 * t44 + t12 * t34));
        J[12] = 0;
        J[13] = -t30 * (t1 + t46 * (t3 * t35 + t33 * t44));
        J[14] = 0;
        J[15] = t38 * (t3 - t47 * t48);
        J[16] = 0;
        J[17] = t47 * t49;
        J[18] = t54 * (t56 - 1);
        J[19] = -t28 * t59;
        J[20] = t62;
        J[21] = t30 * (t37 * (t34 * t57 + t35 * t63) + t63);
        J[22] = t64;
        J[23] = t38 * (-t40 * t59 + t57);
        J[24] = t62;
        J[25] = t38 * (-t42 * t66 + t58);
        J[26] = t54 * (t67 - 1);
        J[27] = t43 * t66;
        J[28] = t68;
        J[29] = -t30 * (t45 * (t34 * (p_x + t2) + t44 * t58) + t65);
        J[30] = t64;
        J[31] = -t30 * (t46 * (t35 * t65 + t44 * (p_y + t0)) + t57);
        J[32] = t68;
        J[33] = t38 * (-t48 * t69 + t65);
        J[34] = t54 * (t70 - 1);
        J[35] = t49 * t69;
        J[36] = t54 * (1 - t56);
        J[37] = t28 * t71;
        J[38] = t72;
        J[39] = t38 * (t22 - t48 * t71);
        J[40] = t73;
        J[41] = -t30 * (t45 * (t22 * t35 + t32 * t34) + t5);
        J[42] = t72;
        J[43] = -t30 * (t22 + t46 * (t14 * t44 + t34 * t7));
        J[44] = t54 * (1 - t67);
        J[45] = t43 * t74;
        J[46] = t75;
        J[47] = t38 * (-t40 * t74 + t7);
        J[48] = t73;
        J[49] = t38 * (-t42 * t76 + t5);
        J[50] = t75;
        J[51] = t30 * (t16 + t37 * (t16 * t35 + t44 * t5));
        J[52] = t54 * (1 - t70);
        J[53] = -t49 * t76;
    }

    void edge_edge_tangent_basis_jacobian(
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
        double J[72])
    {
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14,
            t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,
            t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40,
            t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53,
            t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66,
            t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79,
            t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92,
            t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104,
            t105, t106, t107, t108, t109, t110, t111, t112, t113, t114, t115,
            t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, t126,
            t127, t128, t129, t130, t131;
        t0 = ea0_x - ea1_x;
        t1 = std::pow(t0, 2);
        t2 = ea0_y - ea1_y;
        t3 = std::pow(t2, 2);
        t4 = ea0_z - ea1_z;
        t5 = std::pow(t4, 2);
        t6 = t3 + t5;
        t7 = t1 + t6;
        t8 = std::pow(t7, -1.0 / 2.0);
        t9 = 1.0 / t7;
        t10 = t1 * t9;
        t11 = -ea0_x + ea1_x;
        t12 = -eb0_z + eb1_z;
        t13 = t11 * t12;
        t14 = -ea0_z + ea1_z;
        t15 = -eb0_x + eb1_x;
        t16 = t14 * t15;
        t17 = t13 - t16;
        t18 = eb0_z - eb1_z;
        t19 = t18 * t2;
        t20 = eb0_y - eb1_y;
        t21 = t20 * t4;
        t22 = -t21;
        t23 = t19 + t22;
        t24 = t0 * t17 + t2 * t23;
        t25 = t0 * t20;
        t26 = eb0_x - eb1_x;
        t27 = t2 * t26;
        t28 = -t27;
        t29 = t25 + t28;
        t30 = t17 * t4 + t2 * t29;
        t31 = t0 * t29;
        t32 = t23 * t4;
        t33 = t31 - t32;
        t34 = std::pow(t24, 2) + std::pow(t30, 2) + std::pow(t33, 2);
        t35 = std::pow(t34, -1.0 / 2.0);
        t36 = 1.0 / t34;
        t37 = 2 * t25;
        t38 = t27 - t37;
        t39 = -t31 + t32;
        t40 = t2 * t20;
        t41 = t18 * t4;
        t42 = t40 + t41;
        t43 = -ea0_y + ea1_y;
        t44 = -eb0_y + eb1_y;
        t45 = t11 * t44;
        t46 = t15 * t43;
        t47 = t45 - t46;
        t48 = t43 * t47;
        t49 = -t13;
        t50 = t16 + t49;
        t51 = t14 * t50;
        t52 = t48 - t51;
        t53 = t42 * t52;
        t54 = t26 * t4;
        t55 = t0 * t18;
        t56 = 2 * t55;
        t57 = t54 - t56;
        t58 = t11 * t50;
        t59 = t12 * t43;
        t60 = t14 * t44;
        t61 = -t60;
        t62 = t59 + t61;
        t63 = t43 * t62;
        t64 = t58 - t63;
        t65 = t36 * (t38 * t39 - t53 - t57 * t64);
        t66 = std::pow(t7, -3.0 / 2.0);
        t67 = t0 * t2;
        t68 = t66 * t67;
        t69 = t0 * t4;
        t70 = t66 * t69;
        t71 = t0 * t26;
        t72 = t41 + t71;
        t73 = t39 * t72;
        t74 = 2 * t19;
        t75 = t21 - t74;
        t76 = 2 * t27;
        t77 = t25 - t76;
        t78 = t36 * (-t52 * t77 - t64 * t75 + t73);
        t79 = t3 * t9;
        t80 = t2 * t4;
        t81 = t66 * t80;
        t82 = t40 + t71;
        t83 = t64 * t82;
        t84 = 2 * t21;
        t85 = t19 - t84;
        t86 = -2 * t54 + t55;
        t87 = t36 * (t39 * t85 - t52 * t86 - t83);
        t88 = std::pow(t39, 2) + std::pow(t52, 2) + std::pow(t64, 2);
        t89 = std::pow(t88, -1.0 / 2.0);
        t90 = -t58 + t63;
        t91 = t11 * t47 - t14 * t62;
        t92 = -t48 + t51;
        t93 = 1.0 / t88;
        t94 = t39 * t93;
        t95 = t5 * t9;
        t96 = t28 + t37;
        t97 = -t54 + t56;
        t98 = t36 * (t39 * t96 + t53 - t64 * t97);
        t99 = -t68;
        t100 = -t70;
        t101 = t22 + t74;
        t102 = -t25 + t76;
        t103 = t36 * (t101 * t64 + t102 * t52 + t73);
        t104 = -t81;
        t105 = 2 * t16;
        t106 = t105 + t49;
        t107 =
            t106 * t92 + t90 * (t0 * t15 + t20 * t43) + t91 * (t21 - t59 + t60);
        t108 = t52 * t93;
        t109 = -t19 + t84;
        t110 = t24 * t36;
        t111 = t39 * t67;
        t112 = t64 * t69;
        t113 = t52 * t6;
        t114 = t111 - t112 + t113;
        t115 = t114 * t36;
        t116 = t1 + t5;
        t117 = t116 * t39 + t52 * t67 + t64 * t80;
        t118 = t30 * t36;
        t119 = t14 * t90;
        t120 = t39 * t80;
        t121 = t52 * t69;
        t122 = t1 + t3;
        t123 = t122 * t64;
        t124 = t120 - t121 + t123;
        t125 = t11 * t2;
        t126 = t0 * t119 + t125 * t91 + t92 * (std::pow(t14, 2) + t3);
        t127 = -t69;
        t128 = t64 * t93;
        t129 = t4 * t43;
        t130 = t125 * t92 + t129 * t90 + t91 * (std::pow(t11, 2) + t5);
        t131 = t0 * t14 * t92 + t129 * t91 + t90 * (t1 + std::pow(t43, 2));
        J[0] = t8 * (t10 - 1);
        J[1] = t35 * (-t30 * t65 + t42);
        J[2] = t68;
        J[3] = t35 * (t33 * t65 + t38);
        J[4] = t70;
        J[5] = t35 * (t24 * t65 + t57);
        J[6] = t68;
        J[7] = t35 * (-t30 * t78 + t77);
        J[8] = t8 * (t79 - 1);
        J[9] = t35 * (t33 * t78 + t72);
        J[10] = t81;
        J[11] = t35 * (t24 * t78 + t75);
        J[12] = t70;
        J[13] = t35 * (-t30 * t87 + t86);
        J[14] = t81;
        J[15] = t89
            * (t85
               - t94
                   * (t90 * (t43 * t44 + t71) + t91 * (t4 * t44 + t62)
                      + t92 * (t14 * t26 + t17)));
        J[16] = t8 * (t95 - 1);
        J[17] = t35 * (t24 * t87 + t82);
        J[18] = t8 * (1 - t10);
        J[19] = -t35 * (t30 * t98 + t42);
        J[20] = t99;
        J[21] = t35 * (t33 * t98 + t96);
        J[22] = t100;
        J[23] = t35 * (t24 * t98 + t97);
        J[24] = t99;
        J[25] = t35 * (t102 + t103 * t30);
        J[26] = t8 * (1 - t79);
        J[27] = -t89
            * (t72
               + t94
                   * (t90 * (2 * t59 + t61) + t91 * (t11 * t26 + t12 * t4)
                      + t92 * (t27 - t45 + t46)));
        J[28] = t104;
        J[29] = t35 * (t101 - t103 * t24);
        J[30] = t100;
        J[31] = t89 * (t106 + t107 * t108);
        J[32] = t104;
        J[33] = t89 * (-t107 * t94 + t109);
        J[34] = t8 * (1 - t95);
        J[35] =
            t35 * (t110 * (t109 * t39 + t52 * (-t105 + t13) + t83) - t40 - t71);
        J[36] = 0;
        J[37] = -t35 * (t115 * t30 + t6);
        J[38] = 0;
        J[39] = t35 * (t115 * t33 + t67);
        J[40] = 0;
        J[41] = t35 * (t110 * t114 + t69);
        J[42] = 0;
        J[43] = t35 * (t117 * t118 + t67);
        J[44] = 0;
        J[45] = -t89
            * (t116
               + t94 * (t119 * t43 + t67 * t92 + t91 * (t0 * t11 + t14 * t4)));
        J[46] = 0;
        J[47] = t35 * (-t110 * t117 + t80);
        J[48] = 0;
        J[49] = t35 * (-t118 * t124 + t69);
        J[50] = 0;
        J[51] = t35 * (t124 * t33 * t36 + t80);
        J[52] = 0;
        J[53] = t35 * (-t1 + t110 * t124 - t3);
        J[54] = 0;
        J[55] = t35 * (-t118 * (-t111 + t112 - t113) + t6);
        J[56] = 0;
        J[57] = -t89 * (t126 * t94 + t67);
        J[58] = 0;
        J[59] = t89 * (t126 * t128 + t127);
        J[60] = 0;
        J[61] = t89 * (t108 * t130 - t67);
        J[62] = 0;
        J[63] = t89 * (t116 - t130 * t94);
        J[64] = 0;
        J[65] = t89 * (t128 * t130 - t80);
        J[66] = 0;
        J[67] = t89 * (t108 * t131 + t127);
        J[68] = 0;
        J[69] = -t89 * (t131 * t94 + t80);
        J[70] = 0;
        J[71] = t35 * (t110 * (-t120 + t121 - t123) + t122);
    }

    // This function was generated by the Sympy.
    // 12-Jan-2022 14:03:52
    void point_triangle_tangent_basis_jacobian(
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
        double J[72])
    {
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14,
            t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,
            t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40,
            t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53,
            t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66,
            t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79,
            t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92,
            t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104,
            t105, t106, t107, t108, t109, t110, t111, t112, t113, t114, t115,
            t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, t126,
            t127, t128, t129, t130, t131, t132, t133, t134, t135;

        t0 = t0_x - t1_x;
        t1 = std::pow(t0, 2);
        t2 = -t1_y;
        t3 = t0_y + t2;
        t4 = std::pow(t3, 2);
        t5 = -t1_z;
        t6 = t0_z + t5;
        t7 = std::pow(t6, 2);
        t8 = t4 + t7;
        t9 = t1 + t8;
        t10 = std::pow(t9, -1.0 / 2.0);
        t11 = 1.0 / t9;
        t12 = t1 * t11;
        t13 = -t0_x;
        t14 = t13 + t1_x;
        t15 = -t0_z;
        t16 = t15 + t2_z;
        t17 = t14 * t16;
        t18 = t13 + t2_x;
        t19 = t15 + t1_z;
        t20 = t18 * t19;
        t21 = t17 - t20;
        t22 = -t2_z;
        t23 = t0_z + t22;
        t24 = t23 * t3;
        t25 = -t2_y;
        t26 = t0_y + t25;
        t27 = t26 * t6;
        t28 = -t27;
        t29 = t24 + t28;
        t30 = t0 * t21 + t29 * t3;
        t31 = t0 * t26;
        t32 = -t2_x;
        t33 = t0_x + t32;
        t34 = t3 * t33;
        t35 = -t34;
        t36 = t31 + t35;
        t37 = t21 * t6 + t3 * t36;
        t38 = t0 * t36;
        t39 = t29 * t6;
        t40 = t38 - t39;
        t41 = std::pow(t30, 2) + std::pow(t37, 2) + std::pow(t40, 2);
        t42 = std::pow(t41, -1.0 / 2.0);
        t43 = 1.0 / t41;
        t44 = t1_y + t25;
        t45 = t3 * t44;
        t46 = t1_z + t22;
        t47 = t46 * t6;
        t48 = t45 + t47;
        t49 = -t0_y;
        t50 = t1_y + t49;
        t51 = t2_y + t49;
        t52 = t14 * t51;
        t53 = t18 * t50;
        t54 = t52 - t53;
        t55 = t50 * t54;
        t56 = -t17;
        t57 = t20 + t56;
        t58 = t19 * t57;
        t59 = t55 - t58;
        t60 = -t31;
        t61 = t34 + t60;
        t62 = -t0 * t44 + t61;
        t63 = -t38 + t39;
        t64 = t33 * t6;
        t65 = t0 * t23;
        t66 = t64 - t65;
        t67 = -t0 * t46 + t66;
        t68 = t14 * t57;
        t69 = t16 * t50;
        t70 = t19 * t51;
        t71 = -t70;
        t72 = t69 + t71;
        t73 = t50 * t72;
        t74 = t68 - t73;
        t75 = t43 * (-t48 * t59 + t62 * t63 - t67 * t74);
        t76 = std::pow(t9, -3.0 / 2.0);
        t77 = t0 * t3;
        t78 = t76 * t77;
        t79 = std::pow(t59, 2) + std::pow(t63, 2) + std::pow(t74, 2);
        t80 = std::pow(t79, -1.0 / 2.0);
        t81 = t2_z + t5;
        t82 = -t55 + t58;
        t83 = -t52 + t53;
        t84 = t14 * t54 - t19 * t72;
        t85 = -t68 + t73;
        t86 = 1.0 / t79;
        t87 = t63 * t86;
        t88 = t0 * t6;
        t89 = t76 * t88;
        t90 = t1_x + t32;
        t91 = t3 * t90;
        t92 = t0 * t90;
        t93 = t47 + t92;
        t94 = -t24;
        t95 = t27 - t3 * t46 + t94;
        t96 = t43 * (-t59 * (t36 - t91) + t63 * t93 - t74 * t95);
        t97 = t11 * t4;
        t98 = t3 * t6;
        t99 = t76 * t98;
        t100 = t6 * t90;
        t101 = t45 + t92;
        t102 = t29 - t44 * t6;
        t103 = -t64;
        t104 = t43 * (-t101 * t74 + t102 * t63 - t59 * (-t100 + t103 + t65));
        t105 = t2 + t2_y;
        t106 = t11 * t7;
        t107 = t26 * t3;
        t108 = t23 * t6;
        t109 = t107 + t108;
        t110 = 2 * t31 + t35;
        t111 = t103 + 2 * t65;
        t112 = t43 * (t109 * t59 + t110 * t63 - t111 * t74);
        t113 = -t78;
        t114 = -t89;
        t115 = t0 * t33;
        t116 = t108 + t115;
        t117 = 2 * t24 + t28;
        t118 = 2 * t34 + t60;
        t119 = t43 * (t116 * t63 + t117 * t74 + t118 * t59);
        t120 = -t99;
        t121 = 2 * t20;
        t122 = t121 + t56;
        t123 =
            t122 * t82 + t84 * (t27 - t69 + t70) + t85 * (t0 * t18 + t26 * t50);
        t124 = t59 * t86;
        t125 = 2 * t27 + t94;
        t126 = t30 * t43;
        t127 = t14 * t3;
        t128 = t0 * t19;
        t129 = t127 * t84 + t128 * t85 + t82 * (std::pow(t19, 2) + t4);
        t130 = -t88;
        t131 = t74 * t86;
        t132 = t50 * t6;
        t133 = t127 * t82 + t132 * t85 + t84 * (std::pow(t14, 2) + t7);
        t134 = t128 * t82 + t132 * t84 + t85 * (t1 + std::pow(t50, 2));
        t135 = t1 + t4;

        J[0] = 0;
        J[1] = 0;
        J[2] = 0;
        J[3] = 0;
        J[4] = 0;
        J[5] = 0;
        J[6] = 0;
        J[7] = 0;
        J[8] = 0;
        J[9] = 0;
        J[10] = 0;
        J[11] = 0;
        J[12] = 0;
        J[13] = 0;
        J[14] = 0;
        J[15] = 0;
        J[16] = 0;
        J[17] = 0;
        J[18] = t10 * (t12 - 1);
        J[19] = t42 * (-t37 * t75 + t48);
        J[20] = t78;
        J[21] = t80
            * (t62
               - t87
                   * (t82 * (t19 * t81 + t45) + t84 * (t14 * t44 + t83)
                      + t85 * (t0 * t81 + t57)));
        J[22] = t89;
        J[23] = t42 * (t30 * t75 + t67);
        J[24] = t78;
        J[25] = -t42 * (t37 * t96 + t61 + t91);
        J[26] = t10 * (t97 - 1);
        J[27] = t42 * (t40 * t96 + t93);
        J[28] = t99;
        J[29] = t42 * (t30 * t96 + t95);
        J[30] = t89;
        J[31] = -t42 * (t100 + t104 * t37 + t66);
        J[32] = t99;
        J[33] = t80
            * (t102
               - t87
                   * (t82 * (t19 * t90 + t21) + t84 * (t105 * t6 + t72)
                      + t85 * (t105 * t50 + t92)));
        J[34] = t10 * (t106 - 1);
        J[35] = t42 * (t101 + t104 * t30);
        J[36] = t10 * (1 - t12);
        J[37] = -t42 * (t109 + t112 * t37);
        J[38] = t113;
        J[39] = t42 * (t110 + t112 * t40);
        J[40] = t114;
        J[41] = t42 * (t111 + t112 * t30);
        J[42] = t113;
        J[43] = t42 * (t118 + t119 * t37);
        J[44] = t10 * (1 - t97);
        J[45] = -t80
            * (t116
               + t87
                   * (t82 * (t34 + t83) + t84 * (t14 * t33 + t16 * t6)
                      + t85 * (2 * t69 + t71)));
        J[46] = t120;
        J[47] = t42 * (t117 - t119 * t30);
        J[48] = t114;
        J[49] = t80 * (t122 + t123 * t124);
        J[50] = t120;
        J[51] = t80 * (-t123 * t87 + t125);
        J[52] = t10 * (1 - t106);
        J[53] = t42
            * (-t107 - t115
               + t126
                   * (t125 * t63 + t59 * (-t121 + t17) + t74 * (t107 + t115)));
        J[54] = 0;
        J[55] = t42 * (-t37 * t43 * (-t59 * t8 - t63 * t77 + t74 * t88) + t8);
        J[56] = 0;
        J[57] = -t80 * (t129 * t87 + t77);
        J[58] = 0;
        J[59] = t80 * (t129 * t131 + t130);
        J[60] = 0;
        J[61] = t80 * (t124 * t133 - t77);
        J[62] = 0;
        J[63] = t80 * (t1 - t133 * t87 + t7);
        J[64] = 0;
        J[65] = t80 * (t131 * t133 - t98);
        J[66] = 0;
        J[67] = t80 * (t124 * t134 + t130);
        J[68] = 0;
        J[69] = -t80 * (t134 * t87 + t98);
        J[70] = 0;
        J[71] = t42 * (t126 * (-t135 * t74 + t59 * t88 - t63 * t98) + t135);
    }

} // namespace autogen
} // namespace ipc
