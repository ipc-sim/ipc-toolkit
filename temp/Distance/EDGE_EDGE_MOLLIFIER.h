#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace JGSL {
    
template <class T>
void Edge_Edge_Cross_Norm2(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    T& result)
{
    result = (ea1 - ea0).cross(eb1 - eb0).squaredNorm();
}

template <class T>
void g_EECN2(T v01, T v02, T v03, T v11,
    T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T g[12])
{
    T t8;
    T t9;
    T t10;
    T t11;
    T t12;
    T t13;
    T t23;
    T t24;
    T t25;
    T t26;
    T t27;
    T t28;
    T t29;
    T t30;
    T t31;
    T t32;
    T t33;

    /* COMPUTEEECROSSSQNORMGRADIENT */
    /*     G = COMPUTEEECROSSSQNORMGRADIENT(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     01-Nov-2019 16:54:23 */
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

template <class T>
void Edge_Edge_Cross_Norm2_Gradient(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    Eigen::Matrix<T, 12, 1>& grad)
{
    g_EECN2(
        ea0[0], ea0[1], ea0[2],
        ea1[0], ea1[1], ea1[2],
        eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2],
        grad.data());
}

template <class T>
void H_EECN2(T v01, T v02, T v03, T v11,
    T v12, T v13, T v21, T v22, T v23, T v31, T v32, T v33, T H[144])
{
    T t8;
    T t9;
    T t10;
    T t11;
    T t12;
    T t13;
    T t32;
    T t33;
    T t34;
    T t35;
    T t48;
    T t36;
    T t49;
    T t37;
    T t38;
    T t39;
    T t40;
    T t41;
    T t42;
    T t43;
    T t44;
    T t45;
    T t46;
    T t47;
    T t50;
    T t51;
    T t52;
    T t20;
    T t23;
    T t24;
    T t25;
    T t86;
    T t87;
    T t88;
    T t74;
    T t75;
    T t76;
    T t77;
    T t78;
    T t79;
    T t89;
    T t90;
    T t91;
    T t92;
    T t93;
    T t94;
    T t95;

    /* COMPUTEEECROSSSQNORMHESSIAN */
    /*     H = COMPUTEEECROSSSQNORMHESSIAN(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     01-Nov-2019 16:54:23 */
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

template <class T>
void Edge_Edge_Cross_Norm2_Hessian(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    Eigen::Matrix<T, 12, 12>& Hessian)
{
    H_EECN2(
        ea0[0], ea0[1], ea0[2],
        ea1[0], ea1[1], ea1[2],
        eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2],
        Hessian.data());
}

template <class T>
void derivTest_EECross(void)
{
    Eigen::Matrix<T, 3, 1> v[4];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;
    v[3] << 0.0, 0.1, -1.1;

    T f0;
    Edge_Edge_Cross_Norm2(v[0], v[1], v[2], v[3], f0);

    T eps = 1.0e-6;
    Eigen::Matrix<T, 12, 1> fd, s;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::Matrix<T, 3, 1> vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            T f;
            Edge_Edge_Cross_Norm2(vv[0], vv[1], vv[2], vv[3], f);
            fd[i * 3 + j] = (f - f0) / eps;
        }
    }

    Edge_Edge_Cross_Norm2_Gradient(v[0], v[1], v[2], v[3], s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<T, 12, 12> Hfd, Hs;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::Matrix<T, 3, 1> vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            Eigen::Matrix<T, 12, 1> g;
            Edge_Edge_Cross_Norm2_Gradient(vv[0], vv[1], vv[2], vv[3], g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    Edge_Edge_Cross_Norm2_Hessian(v[0], v[1], v[2], v[3], Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

template <class T>
void EEM(T input, T eps_x, T& e)
{
    T input_div_eps_x = input / eps_x;
    e = (-input_div_eps_x + 2.0) * input_div_eps_x;
}

template <class T>
void g_EEM(T input, T eps_x, T& g)
{
    T one_div_eps_x = 1.0 / eps_x;
    g = 2.0 * one_div_eps_x * (-one_div_eps_x * input + 1.0);
}

template <class T>
void H_EEM(T input, T eps_x, T& H)
{
    H = -2.0 / (eps_x * eps_x);
}

template <class T>
void Edge_Edge_Mollifier(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    T eps_x, T& e)
{
    T EECrossSqNorm;
    Edge_Edge_Cross_Norm2(ea0, ea1, eb0, eb1, EECrossSqNorm);
    if (EECrossSqNorm < eps_x) {
        EEM(EECrossSqNorm, eps_x, e);
    }
    else {
        e = 1.0;
    }
}

template <class T>
void Edge_Edge_Mollifier_Gradient(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    T eps_x, Eigen::Matrix<T, 12, 1>& g)
{
    T EECrossSqNorm;
    Edge_Edge_Cross_Norm2(ea0, ea1, eb0, eb1, EECrossSqNorm);
    if (EECrossSqNorm < eps_x) {
        T q_g;
        g_EEM(EECrossSqNorm, eps_x, q_g);
        Edge_Edge_Cross_Norm2_Gradient(ea0, ea1, eb0, eb1, g);
        g *= q_g;
    }
    else {
        g.setZero();
    }
}

template <class T>
void Edge_Edge_Mollifier_Hessian(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    T eps_x, Eigen::Matrix<T, 12, 12>& H)
{
    T EECrossSqNorm;
    Edge_Edge_Cross_Norm2(ea0, ea1, eb0, eb1, EECrossSqNorm);
    if (EECrossSqNorm < eps_x) {
        T q_g, q_H;
        g_EEM(EECrossSqNorm, eps_x, q_g);
        H_EEM(EECrossSqNorm, eps_x, q_H);

        Eigen::Matrix<T, 12, 1> g;
        Edge_Edge_Cross_Norm2_Gradient(ea0, ea1, eb0, eb1, g);

        Edge_Edge_Cross_Norm2_Hessian(ea0, ea1, eb0, eb1, H);
        H *= q_g;
        H += (q_H * g) * g.transpose();
    }
    else {
        H.setZero();
    }
}

template <class T>
void derivTest_e(T eps_x = 10.0)
{
    Eigen::Matrix<T, 3, 1> v[4];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;
    v[3] << 0.0, 0.1, -1.1;

    T f0;
    Edge_Edge_Mollifier(v[0], v[1], v[2], v[3], eps_x, f0);

    T eps = 1.0e-6;
    Eigen::Matrix<T, 12, 1> fd, s;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::Matrix<T, 3, 1> vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            T f;
            Edge_Edge_Mollifier(vv[0], vv[1], vv[2], vv[3], eps_x, f);
            fd[i * 3 + j] = (f - f0) / eps;
        }
    }

    Edge_Edge_Mollifier_Gradient(v[0], v[1], v[2], v[3], eps_x, s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<T, 12, 12> Hfd, Hs;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::Matrix<T, 3, 1> vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            Eigen::Matrix<T, 12, 1> g;
            Edge_Edge_Mollifier_Gradient(vv[0], vv[1], vv[2], vv[3], eps_x, g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    Edge_Edge_Mollifier_Hessian(v[0], v[1], v[2], v[3], eps_x, Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

template <class T>
void Edge_Edge_Mollifier_Threshold(
    const Eigen::Matrix<T, 3, 1>& ea0_rest,
    const Eigen::Matrix<T, 3, 1>& ea1_rest,
    const Eigen::Matrix<T, 3, 1>& eb0_rest,
    const Eigen::Matrix<T, 3, 1>& eb1_rest, 
    T& eps_x)
{
    eps_x = 1.0e-3 * (ea0_rest - ea1_rest).squaredNorm() * (eb0_rest - eb1_rest).squaredNorm();
}

}