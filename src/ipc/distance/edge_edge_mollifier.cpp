#include "edge_edge_mollifier.hpp"

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

double edge_edge_cross_squarednorm(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    return (ea1 - ea0).cross(eb1 - eb0).squaredNorm();
}

Vector12d edge_edge_cross_squarednorm_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    Vector12d grad;
    autogen::edge_edge_cross_squarednorm_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
    return grad;
}

Matrix12d edge_edge_cross_squarednorm_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    Matrix12d hess;
    autogen::edge_edge_cross_squarednorm_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
    return hess;
}

double edge_edge_mollifier(const double x, const double eps_x)
{
    if (x < eps_x) {
        const double x_div_eps_x = x / eps_x;
        return (-x_div_eps_x + 2) * x_div_eps_x;
    } else {
        return 1;
    }
}

double edge_edge_mollifier_gradient(const double x, const double eps_x)
{
    if (x < eps_x) {
        const double one_div_eps_x = 1 / eps_x;
        return 2 * one_div_eps_x * (-one_div_eps_x * x + 1);
    } else {
        return 0;
    }
}

double edge_edge_mollifier_hessian(const double x, const double eps_x)
{
    if (x < eps_x) {
        return -2 / (eps_x * eps_x);
    } else {
        return 0;
    }
}

double
edge_edge_mollifier_derivative_wrt_eps_x(const double x, const double eps_x)
{
    return x < eps_x ? (2 * x * (-eps_x + x) / (eps_x * eps_x * eps_x)) : 0.0;
}

double edge_edge_mollifier_gradient_derivative_wrt_eps_x(
    const double x, const double eps_x)
{
    return x < eps_x ? (2 * (-eps_x + 2 * x) / (eps_x * eps_x * eps_x)) : 0.0;
}

double edge_edge_mollifier(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    const double eps_x)
{
    const double ee_cross_norm_sqr =
        edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        return edge_edge_mollifier(ee_cross_norm_sqr, eps_x);
    } else {
        return 1;
    }
}

Vector12d edge_edge_mollifier_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    const double eps_x)
{
    const double ee_cross_norm_sqr =
        edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        return edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x)
            * edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1);
    } else {
        return Vector12d::Zero();
    }
}

Matrix12d edge_edge_mollifier_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    const double eps_x)
{
    const double ee_cross_norm_sqr =
        edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        const Vector12d grad =
            edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1);

        return (edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x)
                * edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1))
            + ((edge_edge_mollifier_hessian(ee_cross_norm_sqr, eps_x) * grad)
               * grad.transpose());
    } else {
        return Matrix12d::Zero();
    }
}

Vector12d edge_edge_mollifier_gradient_wrt_x(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    const double eps_x =
        edge_edge_mollifier_threshold(ea0_rest, ea1_rest, eb0_rest, eb1_rest);
    const double ee_cross_norm_sqr =
        edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        // ∇ₓ m = ∂m/∂ε ∇ₓε
        return edge_edge_mollifier_derivative_wrt_eps_x(
                   ee_cross_norm_sqr, eps_x)
            * edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x);
    } else {
        return Vector12d::Zero();
    }
}

Matrix12d edge_edge_mollifier_gradient_jacobian_wrt_x(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    const double eps_x =
        edge_edge_mollifier_threshold(ea0_rest, ea1_rest, eb0_rest, eb1_rest);
    const double ee_cross_norm_sqr =
        edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    if (ee_cross_norm_sqr < eps_x) {
        // ∂²m/∂ε∂s (∇ₓε)(∇ᵤs(x+u))ᵀ + ∂m/∂s ∇ᵤ²s(x+u)
        return edge_edge_mollifier_gradient_derivative_wrt_eps_x(
                   ee_cross_norm_sqr, eps_x)
            * edge_edge_mollifier_threshold_gradient(
                   ea0_rest, ea1_rest, eb0_rest, eb1_rest)
            * edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1)
                  .transpose()
            + edge_edge_mollifier_gradient(ee_cross_norm_sqr, eps_x)
            * edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1);
    } else {
        return Matrix12d::Zero();
    }
}

double edge_edge_mollifier_threshold(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest)
{
    return 1e-3 * (ea0_rest - ea1_rest).squaredNorm()
        * (eb0_rest - eb1_rest).squaredNorm();
}

Vector12d edge_edge_mollifier_threshold_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& ea1_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb0_rest,
    const Eigen::Ref<const Eigen::Vector3d>& eb1_rest)
{
    Vector12d grad;
    autogen::edge_edge_mollifier_threshold_gradient(
        ea0_rest[0], ea0_rest[1], ea0_rest[2], ea1_rest[0], ea1_rest[1],
        ea1_rest[2], eb0_rest[0], eb0_rest[1], eb0_rest[2], eb1_rest[0],
        eb1_rest[1], eb1_rest[2], grad.data(), /*scale=*/1e-3);
    return grad;
}

namespace autogen {

    // This function was generated by the Symbolic Math Toolbox version 8.3.
    // 01-Nov-2019 16:54:23
    void edge_edge_cross_squarednorm_gradient(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double g[12])
    {
        double t8, t9, t10, t11, t12, t13, t23, t24, t25, t26, t27, t28, t29,
            t30, t31, t32, t33;

        t8 = -v11 + v01;
        t9 = -v12 + v02;
        t10 = -v13 + v03;
        t11 = -v31 + v21;
        t12 = -v32 + v22;
        t13 = -v33 + v23;
        t23 = t8 * t12 + -(t9 * t11);
        t24 = t8 * t13 + -(t10 * t11);
        t25 = t9 * t13 + -(t10 * t12);
        t26 = t8 * t23 * 2.0;
        t27 = t9 * t23 * 2.0;
        t28 = t8 * t24 * 2.0;
        t29 = t10 * t24 * 2.0;
        t30 = t9 * t25 * 2.0;
        t31 = t10 * t25 * 2.0;
        t32 = t11 * t23 * 2.0;
        t33 = t12 * t23 * 2.0;
        t23 = t11 * t24 * 2.0;
        t10 = t13 * t24 * 2.0;
        t9 = t12 * t25 * 2.0;
        t8 = t13 * t25 * 2.0;
        g[0] = t33 + t10;
        g[1] = -t32 + t8;
        g[2] = -t23 - t9;
        g[3] = -t33 - t10;
        g[4] = t32 - t8;
        g[5] = t23 + t9;
        g[6] = -t27 - t29;
        g[7] = t26 - t31;
        g[8] = t28 + t30;
        g[9] = t27 + t29;
        g[10] = -t26 + t31;
        g[11] = -t28 - t30;
    }

    // This function was generated by the Symbolic Math Toolbox version 8.3.
    // 01-Nov-2019 16:54:23
    void edge_edge_cross_squarednorm_hessian(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double H[144])
    {
        double t8, t9, t10, t11, t12, t13, t32, t33, t34, t35, t48, t36, t49,
            t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t50, t51,
            t52, t20, t23, t24, t25, t86, t87, t88, t74, t75, t76, t77, t78,
            t79, t89, t90, t91, t92, t93, t94, t95;

        t8 = -v11 + v01;
        t9 = -v12 + v02;
        t10 = -v13 + v03;
        t11 = -v31 + v21;
        t12 = -v32 + v22;
        t13 = -v33 + v23;
        t32 = t8 * t9 * 2.0;
        t33 = t8 * t10 * 2.0;
        t34 = t9 * t10 * 2.0;
        t35 = t8 * t11 * 2.0;
        t48 = t8 * t12;
        t36 = t48 * 2.0;
        t49 = t9 * t11;
        t37 = t49 * 2.0;
        t38 = t48 * 4.0;
        t48 = t8 * t13;
        t39 = t48 * 2.0;
        t40 = t49 * 4.0;
        t41 = t9 * t12 * 2.0;
        t49 = t10 * t11;
        t42 = t49 * 2.0;
        t43 = t48 * 4.0;
        t48 = t9 * t13;
        t44 = t48 * 2.0;
        t45 = t49 * 4.0;
        t49 = t10 * t12;
        t46 = t49 * 2.0;
        t47 = t48 * 4.0;
        t48 = t49 * 4.0;
        t49 = t10 * t13 * 2.0;
        t50 = t11 * t12 * 2.0;
        t51 = t11 * t13 * 2.0;
        t52 = t12 * t13 * 2.0;
        t20 = t8 * t8 * 2.0;
        t9 = t9 * t9 * 2.0;
        t8 = t10 * t10 * 2.0;
        t23 = t11 * t11 * 2.0;
        t24 = t12 * t12 * 2.0;
        t25 = t13 * t13 * 2.0;
        t86 = t35 + t41;
        t87 = t35 + t49;
        t88 = t41 + t49;
        t74 = t20 + t9;
        t75 = t20 + t8;
        t76 = t9 + t8;
        t77 = t23 + t24;
        t78 = t23 + t25;
        t79 = t24 + t25;
        t89 = t40 + -t36;
        t90 = t36 + -t40;
        t91 = t37 + -t38;
        t92 = t38 + -t37;
        t93 = t45 + -t39;
        t94 = t39 + -t45;
        t95 = t42 + -t43;
        t37 = t43 + -t42;
        t39 = t48 + -t44;
        t45 = t44 + -t48;
        t38 = t46 + -t47;
        t40 = t47 + -t46;
        t36 = -t35 + -t41;
        t13 = -t35 + -t49;
        t11 = -t41 + -t49;
        t12 = -t20 + -t9;
        t10 = -t20 + -t8;
        t8 = -t9 + -t8;
        t9 = -t23 + -t24;
        t49 = -t23 + -t25;
        t48 = -t24 + -t25;
        H[0] = t79;
        H[1] = -t50;
        H[2] = -t51;
        H[3] = t48;
        H[4] = t50;
        H[5] = t51;
        H[6] = t11;
        H[7] = t92;
        H[8] = t37;
        H[9] = t88;
        H[10] = t91;
        H[11] = t95;
        H[12] = -t50;
        H[13] = t78;
        H[14] = -t52;
        H[15] = t50;
        H[16] = t49;
        H[17] = t52;
        H[18] = t89;
        H[19] = t13;
        H[20] = t40;
        H[21] = t90;
        H[22] = t87;
        H[23] = t38;
        H[24] = -t51;
        H[25] = -t52;
        H[26] = t77;
        H[27] = t51;
        H[28] = t52;
        H[29] = t9;
        H[30] = t93;
        H[31] = t39;
        H[32] = t36;
        H[33] = t94;
        H[34] = t45;
        H[35] = t86;
        H[36] = t48;
        H[37] = t50;
        H[38] = t51;
        H[39] = t79;
        H[40] = -t50;
        H[41] = -t51;
        H[42] = t88;
        H[43] = t91;
        H[44] = t95;
        H[45] = t11;
        H[46] = t92;
        H[47] = t37;
        H[48] = t50;
        H[49] = t49;
        H[50] = t52;
        H[51] = -t50;
        H[52] = t78;
        H[53] = -t52;
        H[54] = t90;
        H[55] = t87;
        H[56] = t38;
        H[57] = t89;
        H[58] = t13;
        H[59] = t40;
        H[60] = t51;
        H[61] = t52;
        H[62] = t9;
        H[63] = -t51;
        H[64] = -t52;
        H[65] = t77;
        H[66] = t94;
        H[67] = t45;
        H[68] = t86;
        H[69] = t93;
        H[70] = t39;
        H[71] = t36;
        H[72] = t11;
        H[73] = t89;
        H[74] = t93;
        H[75] = t88;
        H[76] = t90;
        H[77] = t94;
        H[78] = t76;
        H[79] = -t32;
        H[80] = -t33;
        H[81] = t8;
        H[82] = t32;
        H[83] = t33;
        H[84] = t92;
        H[85] = t13;
        H[86] = t39;
        H[87] = t91;
        H[88] = t87;
        H[89] = t45;
        H[90] = -t32;
        H[91] = t75;
        H[92] = -t34;
        H[93] = t32;
        H[94] = t10;
        H[95] = t34;
        H[96] = t37;
        H[97] = t40;
        H[98] = t36;
        H[99] = t95;
        H[100] = t38;
        H[101] = t86;
        H[102] = -t33;
        H[103] = -t34;
        H[104] = t74;
        H[105] = t33;
        H[106] = t34;
        H[107] = t12;
        H[108] = t88;
        H[109] = t90;
        H[110] = t94;
        H[111] = t11;
        H[112] = t89;
        H[113] = t93;
        H[114] = t8;
        H[115] = t32;
        H[116] = t33;
        H[117] = t76;
        H[118] = -t32;
        H[119] = -t33;
        H[120] = t91;
        H[121] = t87;
        H[122] = t45;
        H[123] = t92;
        H[124] = t13;
        H[125] = t39;
        H[126] = t32;
        H[127] = t10;
        H[128] = t34;
        H[129] = -t32;
        H[130] = t75;
        H[131] = -t34;
        H[132] = t95;
        H[133] = t38;
        H[134] = t86;
        H[135] = t37;
        H[136] = t40;
        H[137] = t36;
        H[138] = t33;
        H[139] = t34;
        H[140] = t12;
        H[141] = -t33;
        H[142] = -t34;
        H[143] = t74;
    }

    void edge_edge_mollifier_threshold_gradient(
        double ea0x,
        double ea0y,
        double ea0z,
        double ea1x,
        double ea1y,
        double ea1z,
        double eb0x,
        double eb0y,
        double eb0z,
        double eb1x,
        double eb1y,
        double eb1z,
        double grad[12],
        double scale)
    {
        const auto t0 = ea0x - ea1x;
        const auto t1 = eb0x - eb1x;
        const auto t2 = eb0y - eb1y;
        const auto t3 = eb0z - eb1z;
        const auto t4 = 2 * scale;
        const auto t5 = t4 * ((t1 * t1) + (t2 * t2) + (t3 * t3));
        const auto t6 = t0 * t5;
        const auto t7 = ea0y - ea1y;
        const auto t8 = t5 * t7;
        const auto t9 = ea0z - ea1z;
        const auto t10 = t5 * t9;
        const auto t11 = t4 * ((t0 * t0) + (t7 * t7) + (t9 * t9));
        const auto t12 = t1 * t11;
        const auto t13 = t11 * t2;
        const auto t14 = t11 * t3;
        grad[0] = t6;
        grad[1] = t8;
        grad[2] = t10;
        grad[3] = -t6;
        grad[4] = -t8;
        grad[5] = -t10;
        grad[6] = t12;
        grad[7] = t13;
        grad[8] = t14;
        grad[9] = -t12;
        grad[10] = -t13;
        grad[11] = -t14;
    }

} // namespace autogen
} // namespace ipc
