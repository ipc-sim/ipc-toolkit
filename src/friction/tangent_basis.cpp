#include <ipc/friction/tangent_basis.hpp>

namespace ipc {
namespace autogen {

    void point_point_tangent_basis_2D_jacobian(
        double p0_x, double p0_y, double p1_x, double p1_y, double J[8])
    {
        const auto t0 = p0_x - p1_x;
        const auto t1 = p0_y - p1_y;
        const auto t2 = std::pow(t0, 2);
        const auto t3 = std::pow(t1, 2);
        const auto t4 = t2 + t3;
        const auto t5 = t0 * t1 / std::pow(t4, 3.0 / 2.0);
        const auto t6 = -t5;
        const auto t7 = std::pow(t4, -1.0 / 2.0);
        const auto t8 = 1.0 / t4;
        const auto t9 = t2 * t8;
        const auto t10 = t3 * t8;
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
        const auto t8 = ((t6) ? (0) : (t7));
        const auto t9 =
            ((t6) ? (0) : (t3)) + ((t6) ? (t3) : (0)) + ((t6) ? (t5) : (t1));
        const auto t10 = std::pow(t9, -3.0 / 2.0);
        const auto t11 = (1.0 / 2.0) * ((t6) ? (0) : (2 * t0));
        const auto t12 = t10 * t11;
        const auto t13 = ((t6) ? (t2) : (0));
        const auto t14 = std::pow(t9, -1.0 / 2.0);
        const auto t15 = ((t6) ? (0) : (1));
        const auto t16 = -p0_y + p1_y;
        const auto t17 = ((t6) ? (t16) : (t0));
        const auto t18 = 1.0 / t9;
        const auto t19 = t17 * t18;
        const auto t20 = t4 * t8;
        const auto t21 = t0 * t13;
        const auto t22 = t20 - t21;
        const auto t23 = t1 + std::pow(t7, 2) < std::pow(t16, 2) + t3;
        const auto t24 = ((t23) ? (0) : (t7));
        const auto t25 = t24 * t7;
        const auto t26 = -p0_x + p1_x;
        const auto t27 = ((t23) ? (t16) : (t0));
        const auto t28 = t26 * t27;
        const auto t29 = -t25 + t28;
        const auto t30 = t13 * t2;
        const auto t31 = t17 * t4;
        const auto t32 = t30 - t31;
        const auto t33 = std::pow(t22, 2) + std::pow(t29, 2) + std::pow(t32, 2);
        const auto t34 = std::pow(t33, -1.0 / 2.0);
        const auto t35 = t15 * t4;
        const auto t36 = 1.0 / t33;
        const auto t37 = ((t23) ? (t2) : (0));
        const auto t38 = -t16 * t24 + t26 * t37;
        const auto t39 = t37 * t38;
        const auto t40 = ((t23) ? (0) : (1));
        const auto t41 = t16 * t27 - t37 * t7;
        const auto t42 = t16 * t41;
        const auto t43 = t25 - t28;
        const auto t44 = t36 * (-t39 + t40 * t42 + t43 * (t0 * t40 + t27));
        const auto t45 = -t20 + t21;
        const auto t46 = -t0 * t17 + t2 * t8;
        const auto t47 = -t30 + t31;
        const auto t48 = std::pow(t45, 2) + std::pow(t46, 2) + std::pow(t47, 2);
        const auto t49 = std::pow(t48, -1.0 / 2.0);
        const auto t50 = t13 * t22;
        const auto t51 = t0 * t15 + t17;
        const auto t52 = 1.0 / t48;
        const auto t53 = t46 * t52;
        const auto t54 = (1.0 / 2.0) * ((t6) ? (2 * t4) : (0));
        const auto t55 = t10 * t54;
        const auto t56 = ((t6) ? (-1) : (0));
        const auto t57 = t24 * t38;
        const auto t58 = ((t23) ? (-1) : (0));
        const auto t59 = t0 * t43;
        const auto t60 = -t27;
        const auto t61 = t36 * (t41 * (t16 * t58 + t60) + t57 + t58 * t59);
        const auto t62 = t17 + t4 * t56;
        const auto t63 = t0 * t56;
        const auto t64 = t22 * t8;
        const auto t65 = ((t6) ? (0) : (-1));
        const auto t66 = (1.0 / 2.0) * t8;
        const auto t67 = 2 * t2;
        const auto t68 = t18 * (((t23) ? (0) : (t67)) + ((t23) ? (t67) : (0)));
        const auto t69 = ((t6) ? (1) : (0));
        const auto t70 = (1.0 / 2.0) * t13;
        const auto t71 = (1.0 / 2.0) * t10 * t17;
        const auto t72 = ((t23) ? (1) : (0));
        const auto t73 = ((t23) ? (0) : (-1));
        const auto t74 = t32 * t36;
        const auto t75 = t13 + t2 * t69;
        const auto t76 = t2 * t65;
        const auto t77 = t4 * t65;
        const auto t78 = t0 * t69;
        const auto t79 = t77 - t78;
        const auto t80 = t22 * t79 + t29 * (t76 + t8) + t32 * t75;
        const auto t81 = t45 * t52;
        const auto t82 = ((t6) ? (0) : (2 * t26));
        const auto t83 = t10 * t82;
        const auto t84 = (1.0 / 2.0) * t19;
        const auto t85 = t39 + t42 * t73 + t43 * (t0 * t73 + t60);
        const auto t86 = t0 * t65 - t17;
        const auto t87 = t22 * t36;
        const auto t88 = ((t6) ? (2 * t16) : (0));
        const auto t89 = t10 * t88;
        const auto t90 = t17 - t4 * t69;
        const auto t91 = -t29 * t78 + t32 * t90 - t64;
        const auto t92 = 2 * t7;
        const auto t93 = t18 * (((t23) ? (0) : (t92)) + ((t23) ? (t92) : (0)));
        const auto t94 = -t13 + t2 * t56;
        const auto t95 = t35 - t63;
        const auto t96 = -t15 * t2 + t8;
        const auto t97 = t22 * t95 - t29 * t96 + t32 * t94;
        J[0] = -t12 * t8;
        J[1] = -t12 * t13;
        J[2] = t14 * (-t11 * t19 + t15);
        J[3] = -t34 * (t32 * t44 + t35);
        J[4] = t49 * (t51 - t53 * (t29 * t51 + t32 * t35 + t50));
        J[5] = -t34 * (t13 + t22 * t44);
        J[6] = -t55 * t8;
        J[7] = -t13 * t55;
        J[8] = t14 * (-t19 * t54 + t56);
        J[9] = -t34 * (t32 * t61 + t62);
        J[10] = t49 * (t53 * (-t29 * t63 - t32 * t62 + t64) + t63);
        J[11] = t34 * (-t22 * t61 + t8);
        J[12] = t14 * (t65 - t66 * t68);
        J[13] = t14 * (-t68 * t70 + t69);
        J[14] = -t71 * (((t6) ? (0) : (t67)) + ((t6) ? (t67) : (0)));
        J[15] = t34
            * (-t74
                   * (t38 * (t26 * t72 + t4 * t73) + t41 * (t2 * t72 + t37)
                      + t43 * (-t24 + t7 * t73))
               + t75);
        J[16] = t49 * (t53 * t80 - t76 - t8);
        J[17] = t49 * (t79 + t80 * t81);
        J[18] = -t66 * t83;
        J[19] = -t70 * t83;
        J[20] = t14 * (t65 - t82 * t84);
        J[21] = -t34 * (t74 * t85 + t77);
        J[22] = t49 * (t53 * (-t29 * t86 - t32 * t77 + t50) + t86);
        J[23] = t34 * (t13 - t85 * t87);
        J[24] = -t66 * t89;
        J[25] = -t70 * t89;
        J[26] = t14 * (t69 - t84 * t88);
        J[27] = t49 * (t47 * t52 * t91 + t90);
        J[28] = t49 * (t53 * t91 + t78);
        J[29] = -t34 * (t8 + t87 * (t41 * (t16 * t72 + t27) - t57 + t59 * t72));
        J[30] = t14 * (t15 - t66 * t93);
        J[31] = t14 * (t56 - t70 * t93);
        J[32] = -t71 * (((t6) ? (0) : (t92)) + ((t6) ? (t92) : (0)));
        J[33] = t34
            * (-t74
                   * (t38 * (t26 * t58 + t4 * t40) + t41 * (t2 * t58 - t37)
                      + t43 * (t24 + t40 * t7))
               + t94);
        J[34] = t49 * (t53 * t97 + t96);
        J[35] = t49 * (t81 * t97 + t95);
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
        const auto t0 = e0_x - e1_x;
        const auto t1 = std::pow(t0, 2);
        const auto t2 = e0_y - e1_y;
        const auto t3 = std::pow(t2, 2);
        const auto t4 = t1 + t3;
        const auto t5 = std::pow(t4, -1.0 / 2.0);
        const auto t6 = 1.0 / t4;
        const auto t7 = t1 * t6;
        const auto t8 = t0 * t2 / std::pow(t4, 3.0 / 2.0);
        const auto t9 = t3 * t6;
        const auto t10 = -t8;
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
        const auto t0 = -e1_y;
        const auto t1 = e0_y + t0;
        const auto t2 = -e1_x;
        const auto t3 = e0_x + t2;
        const auto t4 = -p_y;
        const auto t5 = e0_y + t4;
        const auto t6 = -p_x;
        const auto t7 = e0_x + t6;
        const auto t8 = -t1 * t7 + t3 * t5;
        const auto t9 = -e1_z;
        const auto t10 = e0_z + t9;
        const auto t11 = -e0_x;
        const auto t12 = e1_x + t11;
        const auto t13 = -e0_z;
        const auto t14 = p_z + t13;
        const auto t15 = t12 * t14;
        const auto t16 = p_x + t11;
        const auto t17 = e1_z + t13;
        const auto t18 = t16 * t17;
        const auto t19 = t15 - t18;
        const auto t20 = t1 * t8 + t10 * t19;
        const auto t21 = -p_z;
        const auto t22 = e0_z + t21;
        const auto t23 = t1 * t22 - t10 * t5;
        const auto t24 = -t10 * t7 + t22 * t3;
        const auto t25 = std::pow(t23, 2) + std::pow(t8, 2);
        const auto t26 = std::pow(t24, 2) + t25;
        const auto t27 = std::pow(t26, -3.0 / 2.0);
        const auto t28 = t23 * t27;
        const auto t29 = std::pow(t19, 2) + t25;
        const auto t30 = std::pow(t29, -1.0 / 2.0);
        const auto t31 = -e0_y;
        const auto t32 = p_y + t31;
        const auto t33 = e1_y + t31;
        const auto t34 = t12 * t32 - t16 * t33;
        const auto t35 = -t15 + t18;
        const auto t36 = 1.0 / t29;
        const auto t37 = t19 * t36;
        const auto t38 = std::pow(t26, -1.0 / 2.0);
        const auto t39 = 1.0 / t26;
        const auto t40 = t39 * t8;
        const auto t41 = t10 * t23 - t3 * t8;
        const auto t42 = t23 * t39;
        const auto t43 = t24 * t27;
        const auto t44 = t14 * t33 - t17 * t32;
        const auto t45 = t36 * t8;
        const auto t46 = t23 * t36;
        const auto t47 = t1 * t23 + t19 * t3;
        const auto t48 = t24 * t39;
        const auto t49 = t27 * t8;
        const auto t50 = std::pow(t3, 2);
        const auto t51 = std::pow(t1, 2);
        const auto t52 = std::pow(t10, 2);
        const auto t53 = t50 + t51 + t52;
        const auto t54 = std::pow(t53, -1.0 / 2.0);
        const auto t55 = 1.0 / t53;
        const auto t56 = t50 * t55;
        const auto t57 = std::pow(t53, -3.0 / 2.0);
        const auto t58 = t3 * t57;
        const auto t59 = t1 * t58;
        const auto t60 = t10 * t58;
        const auto t61 = e1_y + t4;
        const auto t62 = e1_z + t21;
        const auto t63 = t19 * t62 + t61 * t8;
        const auto t64 = p_z + t9;
        const auto t65 = t51 * t55;
        const auto t66 = t1 * t10 * t57;
        const auto t67 = e1_x + t6;
        const auto t68 = t23 * t62 - t67 * t8;
        const auto t69 = t52 * t55;
        const auto t70 = t19 * t67 + t23 * t61;
        const auto t71 = -t59;
        const auto t72 = -t60;
        const auto t73 = t19 * t22 + t5 * t8;
        const auto t74 = -t66;
        const auto t75 = -t22 * t23 + t7 * t8;
        const auto t76 = t19 * t7 + t23 * t5;
        J[0] = 0;
        J[1] = 0;
        J[2] = 0;
        J[3] = -t20 * t28;
        J[4] = t30 * (t17 + t37 * (t1 * t34 + t17 * t35));
        J[5] = t38 * (t1 - t20 * t40);
        J[6] = 0;
        J[7] = 0;
        J[8] = 0;
        J[9] = t38 * (t10 - t41 * t42);
        J[10] = t41 * t43;
        J[11] = -t30 * (t3 + t45 * (t10 * t44 + t12 * t34));
        J[12] = 0;
        J[13] = 0;
        J[14] = 0;
        J[15] = -t30 * (t1 + t46 * (t3 * t35 + t33 * t44));
        J[16] = t38 * (t3 - t47 * t48);
        J[17] = t47 * t49;
        J[18] = t54 * (t56 - 1);
        J[19] = t59;
        J[20] = t60;
        J[21] = -t28 * t63;
        J[22] = t30 * (t37 * (t34 * t61 + t35 * t64) + t64);
        J[23] = t38 * (-t40 * t63 + t61);
        J[24] = t59;
        J[25] = t54 * (t65 - 1);
        J[26] = t66;
        J[27] = t38 * (-t42 * t68 + t62);
        J[28] = t43 * t68;
        J[29] = -t30 * (t45 * (t34 * (p_x + t2) + t44 * t62) + t67);
        J[30] = t60;
        J[31] = t66;
        J[32] = t54 * (t69 - 1);
        J[33] = -t30 * (t46 * (t35 * t67 + t44 * (p_y + t0)) + t61);
        J[34] = t38 * (-t48 * t70 + t67);
        J[35] = t49 * t70;
        J[36] = t54 * (1 - t56);
        J[37] = t71;
        J[38] = t72;
        J[39] = t28 * t73;
        J[40] = t38 * (t22 - t48 * t73);
        J[41] = -t30 * (t45 * (t22 * t35 + t32 * t34) + t5);
        J[42] = t71;
        J[43] = t54 * (1 - t65);
        J[44] = t74;
        J[45] = -t30 * (t22 + t46 * (t14 * t44 + t34 * t7));
        J[46] = t43 * t75;
        J[47] = t38 * (-t40 * t75 + t7);
        J[48] = t72;
        J[49] = t74;
        J[50] = t54 * (1 - t69);
        J[51] = t38 * (-t42 * t76 + t5);
        J[52] = t30 * (t16 + t37 * (t16 * t35 + t44 * t5));
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
        const auto t0 = ea0_x - ea1_x;
        const auto t1 = std::pow(t0, 2);
        const auto t2 = ea0_y - ea1_y;
        const auto t3 = std::pow(t2, 2);
        const auto t4 = ea0_z - ea1_z;
        const auto t5 = std::pow(t4, 2);
        const auto t6 = t3 + t5;
        const auto t7 = t1 + t6;
        const auto t8 = std::pow(t7, -1.0 / 2.0);
        const auto t9 = 1.0 / t7;
        const auto t10 = t1 * t9;
        const auto t11 = std::pow(t7, -3.0 / 2.0);
        const auto t12 = t0 * t2;
        const auto t13 = t11 * t12;
        const auto t14 = t0 * t4;
        const auto t15 = t11 * t14;
        const auto t16 = -ea0_x + ea1_x;
        const auto t17 = -eb0_z + eb1_z;
        const auto t18 = t16 * t17;
        const auto t19 = -ea0_z + ea1_z;
        const auto t20 = -eb0_x + eb1_x;
        const auto t21 = t19 * t20;
        const auto t22 = t18 - t21;
        const auto t23 = eb0_z - eb1_z;
        const auto t24 = t2 * t23;
        const auto t25 = eb0_y - eb1_y;
        const auto t26 = t25 * t4;
        const auto t27 = -t26;
        const auto t28 = t24 + t27;
        const auto t29 = t0 * t22 + t2 * t28;
        const auto t30 = t0 * t25;
        const auto t31 = eb0_x - eb1_x;
        const auto t32 = t2 * t31;
        const auto t33 = -t32;
        const auto t34 = t30 + t33;
        const auto t35 = t2 * t34 + t22 * t4;
        const auto t36 = t0 * t34;
        const auto t37 = t28 * t4;
        const auto t38 = t36 - t37;
        const auto t39 = std::pow(t29, 2) + std::pow(t35, 2) + std::pow(t38, 2);
        const auto t40 = std::pow(t39, -1.0 / 2.0);
        const auto t41 = 1.0 / t39;
        const auto t42 = 2 * t30;
        const auto t43 = t32 - t42;
        const auto t44 = -t36 + t37;
        const auto t45 = t2 * t25;
        const auto t46 = t23 * t4;
        const auto t47 = t45 + t46;
        const auto t48 = -ea0_y + ea1_y;
        const auto t49 = -eb0_y + eb1_y;
        const auto t50 = t16 * t49;
        const auto t51 = t20 * t48;
        const auto t52 = t50 - t51;
        const auto t53 = t48 * t52;
        const auto t54 = -t18;
        const auto t55 = t21 + t54;
        const auto t56 = t19 * t55;
        const auto t57 = t53 - t56;
        const auto t58 = t47 * t57;
        const auto t59 = t31 * t4;
        const auto t60 = t0 * t23;
        const auto t61 = 2 * t60;
        const auto t62 = t59 - t61;
        const auto t63 = t16 * t55;
        const auto t64 = t17 * t48;
        const auto t65 = t19 * t49;
        const auto t66 = -t65;
        const auto t67 = t64 + t66;
        const auto t68 = t48 * t67;
        const auto t69 = t63 - t68;
        const auto t70 = t41 * (t43 * t44 - t58 - t62 * t69);
        const auto t71 = t3 * t9;
        const auto t72 = t2 * t4;
        const auto t73 = t11 * t72;
        const auto t74 = t0 * t31;
        const auto t75 = t46 + t74;
        const auto t76 = t44 * t75;
        const auto t77 = 2 * t24;
        const auto t78 = t26 - t77;
        const auto t79 = 2 * t32;
        const auto t80 = t30 - t79;
        const auto t81 = t41 * (-t57 * t80 - t69 * t78 + t76);
        const auto t82 = t5 * t9;
        const auto t83 = t45 + t74;
        const auto t84 = t69 * t83;
        const auto t85 = 2 * t26;
        const auto t86 = t24 - t85;
        const auto t87 = -2 * t59 + t60;
        const auto t88 = t41 * (t44 * t86 - t57 * t87 - t84);
        const auto t89 = std::pow(t44, 2) + std::pow(t57, 2) + std::pow(t69, 2);
        const auto t90 = std::pow(t89, -1.0 / 2.0);
        const auto t91 = -t63 + t68;
        const auto t92 = t16 * t52 - t19 * t67;
        const auto t93 = -t53 + t56;
        const auto t94 = 1.0 / t89;
        const auto t95 = t44 * t94;
        const auto t96 = -t13;
        const auto t97 = -t15;
        const auto t98 = t33 + t42;
        const auto t99 = -t59 + t61;
        const auto t100 = t41 * (t44 * t98 + t58 - t69 * t99);
        const auto t101 = -t73;
        const auto t102 = t27 + t77;
        const auto t103 = -t30 + t79;
        const auto t104 = t41 * (t102 * t69 + t103 * t57 + t76);
        const auto t105 = 2 * t21;
        const auto t106 = t105 + t54;
        const auto t107 =
            t106 * t93 + t91 * (t0 * t20 + t25 * t48) + t92 * (t26 - t64 + t65);
        const auto t108 = t57 * t94;
        const auto t109 = -t24 + t85;
        const auto t110 = t29 * t41;
        const auto t111 = t12 * t44;
        const auto t112 = t14 * t69;
        const auto t113 = t57 * t6;
        const auto t114 = t111 - t112 + t113;
        const auto t115 = t114 * t41;
        const auto t116 = t1 + t5;
        const auto t117 = t116 * t44 + t12 * t57 + t69 * t72;
        const auto t118 = t35 * t41;
        const auto t119 = t19 * t91;
        const auto t120 = t44 * t72;
        const auto t121 = t14 * t57;
        const auto t122 = t1 + t3;
        const auto t123 = t122 * t69;
        const auto t124 = t120 - t121 + t123;
        const auto t125 = t16 * t2;
        const auto t126 =
            t0 * t119 + t125 * t92 + t93 * (std::pow(t19, 2) + t3);
        const auto t127 = -t14;
        const auto t128 = t69 * t94;
        const auto t129 = t4 * t48;
        const auto t130 =
            t125 * t93 + t129 * t91 + t92 * (std::pow(t16, 2) + t5);
        const auto t131 =
            t0 * t19 * t93 + t129 * t92 + t91 * (t1 + std::pow(t48, 2));
        J[0] = t8 * (t10 - 1);
        J[1] = t13;
        J[2] = t15;
        J[3] = t40 * (-t35 * t70 + t47);
        J[4] = t40 * (t38 * t70 + t43);
        J[5] = t40 * (t29 * t70 + t62);
        J[6] = t13;
        J[7] = t8 * (t71 - 1);
        J[8] = t73;
        J[9] = t40 * (-t35 * t81 + t80);
        J[10] = t40 * (t38 * t81 + t75);
        J[11] = t40 * (t29 * t81 + t78);
        J[12] = t15;
        J[13] = t73;
        J[14] = t8 * (t82 - 1);
        J[15] = t40 * (-t35 * t88 + t87);
        J[16] = t90
            * (t86
               - t95
                   * (t91 * (t48 * t49 + t74) + t92 * (t4 * t49 + t67)
                      + t93 * (t19 * t31 + t22)));
        J[17] = t40 * (t29 * t88 + t83);
        J[18] = t8 * (1 - t10);
        J[19] = t96;
        J[20] = t97;
        J[21] = -t40 * (t100 * t35 + t47);
        J[22] = t40 * (t100 * t38 + t98);
        J[23] = t40 * (t100 * t29 + t99);
        J[24] = t96;
        J[25] = t8 * (1 - t71);
        J[26] = t101;
        J[27] = t40 * (t103 + t104 * t35);
        J[28] = -t90
            * (t75
               + t95
                   * (t91 * (2 * t64 + t66) + t92 * (t16 * t31 + t17 * t4)
                      + t93 * (t32 - t50 + t51)));
        J[29] = t40 * (t102 - t104 * t29);
        J[30] = t97;
        J[31] = t101;
        J[32] = t8 * (1 - t82);
        J[33] = t90 * (t106 + t107 * t108);
        J[34] = t90 * (-t107 * t95 + t109);
        J[35] =
            t40 * (t110 * (t109 * t44 + t57 * (-t105 + t18) + t84) - t45 - t74);
        J[36] = 0;
        J[37] = 0;
        J[38] = 0;
        J[39] = -t40 * (t115 * t35 + t6);
        J[40] = t40 * (t115 * t38 + t12);
        J[41] = t40 * (t110 * t114 + t14);
        J[42] = 0;
        J[43] = 0;
        J[44] = 0;
        J[45] = t40 * (t117 * t118 + t12);
        J[46] = -t90
            * (t116
               + t95 * (t119 * t48 + t12 * t93 + t92 * (t0 * t16 + t19 * t4)));
        J[47] = t40 * (-t110 * t117 + t72);
        J[48] = 0;
        J[49] = 0;
        J[50] = 0;
        J[51] = t40 * (-t118 * t124 + t14);
        J[52] = t40 * (t124 * t38 * t41 + t72);
        J[53] = t40 * (-t1 + t110 * t124 - t3);
        J[54] = 0;
        J[55] = 0;
        J[56] = 0;
        J[57] = t40 * (-t118 * (-t111 + t112 - t113) + t6);
        J[58] = -t90 * (t12 + t126 * t95);
        J[59] = t90 * (t126 * t128 + t127);
        J[60] = 0;
        J[61] = 0;
        J[62] = 0;
        J[63] = t90 * (t108 * t130 - t12);
        J[64] = t90 * (t116 - t130 * t95);
        J[65] = t90 * (t128 * t130 - t72);
        J[66] = 0;
        J[67] = 0;
        J[68] = 0;
        J[69] = t90 * (t108 * t131 + t127);
        J[70] = -t90 * (t131 * t95 + t72);
        J[71] = t40 * (t110 * (-t120 + t121 - t123) + t122);
    }

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
        const auto t0 = t0_x - t1_x;
        const auto t1 = std::pow(t0, 2);
        const auto t2 = -t1_y;
        const auto t3 = t0_y + t2;
        const auto t4 = std::pow(t3, 2);
        const auto t5 = -t1_z;
        const auto t6 = t0_z + t5;
        const auto t7 = std::pow(t6, 2);
        const auto t8 = t4 + t7;
        const auto t9 = t1 + t8;
        const auto t10 = std::pow(t9, -1.0 / 2.0);
        const auto t11 = 1.0 / t9;
        const auto t12 = t1 * t11;
        const auto t13 = std::pow(t9, -3.0 / 2.0);
        const auto t14 = t0 * t3;
        const auto t15 = t13 * t14;
        const auto t16 = t0 * t6;
        const auto t17 = t13 * t16;
        const auto t18 = -t0_x;
        const auto t19 = t18 + t1_x;
        const auto t20 = -t0_z;
        const auto t21 = t20 + t2_z;
        const auto t22 = t19 * t21;
        const auto t23 = t18 + t2_x;
        const auto t24 = t1_z + t20;
        const auto t25 = t23 * t24;
        const auto t26 = t22 - t25;
        const auto t27 = -t2_z;
        const auto t28 = t0_z + t27;
        const auto t29 = t28 * t3;
        const auto t30 = -t2_y;
        const auto t31 = t0_y + t30;
        const auto t32 = t31 * t6;
        const auto t33 = -t32;
        const auto t34 = t29 + t33;
        const auto t35 = t0 * t26 + t3 * t34;
        const auto t36 = t0 * t31;
        const auto t37 = -t2_x;
        const auto t38 = t0_x + t37;
        const auto t39 = t3 * t38;
        const auto t40 = -t39;
        const auto t41 = t36 + t40;
        const auto t42 = t26 * t6 + t3 * t41;
        const auto t43 = t0 * t41;
        const auto t44 = t34 * t6;
        const auto t45 = t43 - t44;
        const auto t46 = std::pow(t35, 2) + std::pow(t42, 2) + std::pow(t45, 2);
        const auto t47 = std::pow(t46, -1.0 / 2.0);
        const auto t48 = 1.0 / t46;
        const auto t49 = t1_y + t30;
        const auto t50 = t3 * t49;
        const auto t51 = t1_z + t27;
        const auto t52 = t51 * t6;
        const auto t53 = t50 + t52;
        const auto t54 = -t0_y;
        const auto t55 = t1_y + t54;
        const auto t56 = t2_y + t54;
        const auto t57 = t19 * t56;
        const auto t58 = t23 * t55;
        const auto t59 = t57 - t58;
        const auto t60 = t55 * t59;
        const auto t61 = -t22;
        const auto t62 = t25 + t61;
        const auto t63 = t24 * t62;
        const auto t64 = t60 - t63;
        const auto t65 = -t36;
        const auto t66 = t39 + t65;
        const auto t67 = -t0 * t49 + t66;
        const auto t68 = -t43 + t44;
        const auto t69 = t38 * t6;
        const auto t70 = t0 * t28;
        const auto t71 = t69 - t70;
        const auto t72 = -t0 * t51 + t71;
        const auto t73 = t19 * t62;
        const auto t74 = t21 * t55;
        const auto t75 = t24 * t56;
        const auto t76 = -t75;
        const auto t77 = t74 + t76;
        const auto t78 = t55 * t77;
        const auto t79 = t73 - t78;
        const auto t80 = t48 * (-t53 * t64 + t67 * t68 - t72 * t79);
        const auto t81 = std::pow(t64, 2) + std::pow(t68, 2) + std::pow(t79, 2);
        const auto t82 = std::pow(t81, -1.0 / 2.0);
        const auto t83 = t2_z + t5;
        const auto t84 = -t60 + t63;
        const auto t85 = -t57 + t58;
        const auto t86 = t19 * t59 - t24 * t77;
        const auto t87 = -t73 + t78;
        const auto t88 = 1.0 / t81;
        const auto t89 = t68 * t88;
        const auto t90 = t11 * t4;
        const auto t91 = t3 * t6;
        const auto t92 = t13 * t91;
        const auto t93 = t1_x + t37;
        const auto t94 = t3 * t93;
        const auto t95 = t0 * t93;
        const auto t96 = t52 + t95;
        const auto t97 = -t29;
        const auto t98 = -t3 * t51 + t32 + t97;
        const auto t99 = t48 * (-t64 * (t41 - t94) + t68 * t96 - t79 * t98);
        const auto t100 = t11 * t7;
        const auto t101 = t6 * t93;
        const auto t102 = t50 + t95;
        const auto t103 = t34 - t49 * t6;
        const auto t104 = -t69;
        const auto t105 =
            t48 * (-t102 * t79 + t103 * t68 - t64 * (-t101 + t104 + t70));
        const auto t106 = t2 + t2_y;
        const auto t107 = -t15;
        const auto t108 = -t17;
        const auto t109 = t3 * t31;
        const auto t110 = t28 * t6;
        const auto t111 = t109 + t110;
        const auto t112 = 2 * t36 + t40;
        const auto t113 = t104 + 2 * t70;
        const auto t114 = t48 * (t111 * t64 + t112 * t68 - t113 * t79);
        const auto t115 = -t92;
        const auto t116 = t0 * t38;
        const auto t117 = t110 + t116;
        const auto t118 = 2 * t29 + t33;
        const auto t119 = 2 * t39 + t65;
        const auto t120 = t48 * (t117 * t68 + t118 * t79 + t119 * t64);
        const auto t121 = 2 * t25;
        const auto t122 = t121 + t61;
        const auto t123 =
            t122 * t84 + t86 * (t32 - t74 + t75) + t87 * (t0 * t23 + t31 * t55);
        const auto t124 = t64 * t88;
        const auto t125 = 2 * t32 + t97;
        const auto t126 = t35 * t48;
        const auto t127 = t19 * t3;
        const auto t128 = t0 * t24;
        const auto t129 =
            t127 * t86 + t128 * t87 + t84 * (std::pow(t24, 2) + t4);
        const auto t130 = -t16;
        const auto t131 = t79 * t88;
        const auto t132 = t55 * t6;
        const auto t133 =
            t127 * t84 + t132 * t87 + t86 * (std::pow(t19, 2) + t7);
        const auto t134 =
            t128 * t84 + t132 * t86 + t87 * (t1 + std::pow(t55, 2));
        const auto t135 = t1 + t4;
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
        J[19] = t15;
        J[20] = t17;
        J[21] = t47 * (-t42 * t80 + t53);
        J[22] = t82
            * (t67
               - t89
                   * (t84 * (t24 * t83 + t50) + t86 * (t19 * t49 + t85)
                      + t87 * (t0 * t83 + t62)));
        J[23] = t47 * (t35 * t80 + t72);
        J[24] = t15;
        J[25] = t10 * (t90 - 1);
        J[26] = t92;
        J[27] = -t47 * (t42 * t99 + t66 + t94);
        J[28] = t47 * (t45 * t99 + t96);
        J[29] = t47 * (t35 * t99 + t98);
        J[30] = t17;
        J[31] = t92;
        J[32] = t10 * (t100 - 1);
        J[33] = -t47 * (t101 + t105 * t42 + t71);
        J[34] = t82
            * (t103
               - t89
                   * (t84 * (t24 * t93 + t26) + t86 * (t106 * t6 + t77)
                      + t87 * (t106 * t55 + t95)));
        J[35] = t47 * (t102 + t105 * t35);
        J[36] = t10 * (1 - t12);
        J[37] = t107;
        J[38] = t108;
        J[39] = -t47 * (t111 + t114 * t42);
        J[40] = t47 * (t112 + t114 * t45);
        J[41] = t47 * (t113 + t114 * t35);
        J[42] = t107;
        J[43] = t10 * (1 - t90);
        J[44] = t115;
        J[45] = t47 * (t119 + t120 * t42);
        J[46] = -t82
            * (t117
               + t89
                   * (t84 * (t39 + t85) + t86 * (t19 * t38 + t21 * t6)
                      + t87 * (2 * t74 + t76)));
        J[47] = t47 * (t118 - t120 * t35);
        J[48] = t108;
        J[49] = t115;
        J[50] = t10 * (1 - t100);
        J[51] = t82 * (t122 + t123 * t124);
        J[52] = t82 * (-t123 * t89 + t125);
        J[53] = t47
            * (-t109 - t116
               + t126
                   * (t125 * t68 + t64 * (-t121 + t22) + t79 * (t109 + t116)));
        J[54] = 0;
        J[55] = 0;
        J[56] = 0;
        J[57] = t47 * (-t42 * t48 * (-t14 * t68 + t16 * t79 - t64 * t8) + t8);
        J[58] = -t82 * (t129 * t89 + t14);
        J[59] = t82 * (t129 * t131 + t130);
        J[60] = 0;
        J[61] = 0;
        J[62] = 0;
        J[63] = t82 * (t124 * t133 - t14);
        J[64] = t82 * (t1 - t133 * t89 + t7);
        J[65] = t82 * (t131 * t133 - t91);
        J[66] = 0;
        J[67] = 0;
        J[68] = 0;
        J[69] = t82 * (t124 * t134 + t130);
        J[70] = -t82 * (t134 * t89 + t91);
        J[71] = t47 * (t126 * (-t135 * t79 + t16 * t64 - t68 * t91) + t135);
    }

} // namespace autogen
} // namespace ipc
