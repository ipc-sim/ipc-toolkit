#include "mollifier.hpp"
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include "point_edge.hpp"

namespace ipc {
namespace {
    double func_aux1(const double &a, const double &b, const double &c)
    {
        return (a - c) / b;
    }

    Eigen::Vector3d func_aux1_grad(const double &a, const double &b, const double &c)
    {
        return Eigen::Vector3d(1. / b, -func_aux1(a, b, c) / b, -1. / b);
    }

    Eigen::Matrix3d func_aux1_hess(const double &a, const double &b, const double &c)
    {
        const double b2 = 1. / b / b;
        Eigen::Matrix3d out;
        out << 0.,           -b2, 0.,
              -b2, 2. * func_aux1(a, b, c) * b2, b2,
               0.,            b2, 0.;
        return out;
    }

    double func_aux2(const double &a, const double &b, const double &c)
    {
        return Math<double>::mollifier(func_aux1(a, b, c));
    }

    Eigen::Vector3d func_aux2_grad(const double &a, const double &b, const double &c)
    {
        double val = func_aux1(a, b, c);
        return Math<double>::mollifier_grad(val) * func_aux1_grad(a, b, c);
    }

    Eigen::Matrix3d func_aux2_hess(const double &a, const double &b, const double &c)
    {
        double val = func_aux1(a, b, c);
        Eigen::Vector3d g = func_aux1_grad(a, b, c);
        Eigen::Matrix3d h = func_aux1_hess(a, b, c);
        return g * Math<double>::mollifier_hess(val) * g.transpose() + Math<double>::mollifier_grad(val) * h;
    }
}

/// @brief Compute the gradient of the mollifier function wrt. 4 edge points and the distance squared
std::pair<double, Vector<double, 13>> edge_edge_mollifier_grad(
    const Vector3<double> &ea0, const Vector3<double> &ea1,
    const Vector3<double> &eb0, const Vector3<double> &eb1,
    const double &dist_sqr)
{
    const double da0 = point_edge_distance(ea0, eb0, eb1);
    const double da1 = point_edge_distance(ea1, eb0, eb1);
    const double db0 = point_edge_distance(eb0, ea0, ea1);
    const double db1 = point_edge_distance(eb1, ea0, ea1);

    Vector12d da0_wrt_x, da1_wrt_x, db0_wrt_x, db1_wrt_x;
    {
        Vector9d tmp = point_edge_distance_gradient(ea0, eb0, eb1);
        da0_wrt_x << tmp.head(3), Eigen::Vector3d::Zero(), tmp.tail(6);

        tmp = point_edge_distance_gradient(ea1, eb0, eb1);
        da1_wrt_x << Eigen::Vector3d::Zero(), tmp;

        tmp = point_edge_distance_gradient(eb0, ea0, ea1);
        db0_wrt_x << tmp.tail(6), tmp.head(3), Eigen::Vector3d::Zero();

        tmp = point_edge_distance_gradient(eb1, ea0, ea1);
        db1_wrt_x << tmp.tail(6), Eigen::Vector3d::Zero(), tmp.head(3);
    }

    const double la = point_point_distance(ea0, ea1);
    const double lb = point_point_distance(eb0, eb1);

    const double db = lb * mollifier_threshold_eps;
    const double da = la * mollifier_threshold_eps;

    Vector6d da_wrt_x = point_point_distance_gradient(ea0, ea1) * mollifier_threshold_eps;
    Vector6d db_wrt_x = point_point_distance_gradient(eb0, eb1) * mollifier_threshold_eps;

    const double ma0 = func_aux1(da0, db, dist_sqr);
    const double ma1 = func_aux1(da1, db, dist_sqr);
    const double mb0 = func_aux1(db0, da, dist_sqr);
    const double mb1 = func_aux1(db1, da, dist_sqr);

    const double a = Math<double>::mollifier(ma0);
    const double b = Math<double>::mollifier(ma1);
    const double c = Math<double>::mollifier(mb0);
    const double d = Math<double>::mollifier(mb1);

    const double a_wrt_ma0 = Math<double>::mollifier_grad(ma0);
    const double b_wrt_ma1 = Math<double>::mollifier_grad(ma1);
    const double c_wrt_mb0 = Math<double>::mollifier_grad(mb0);
    const double d_wrt_mb1 = Math<double>::mollifier_grad(mb1);

    const double ab = a * b, cd = c * d;
    Vector12d ab_grad = a * b_wrt_ma1 * da1_wrt_x + b * a_wrt_ma0 * da0_wrt_x;
    ab_grad.tail(6)  -= db_wrt_x * (a * b_wrt_ma1 * ma1 + b * a_wrt_ma0 * ma0);
    ab_grad /= db;
    Vector12d cd_grad = c * d_wrt_mb1 * db1_wrt_x + d * c_wrt_mb0 * db0_wrt_x;
    cd_grad.head(6)  -= da_wrt_x * (c * d_wrt_mb1 * mb1 + d * c_wrt_mb0 * mb0);
    cd_grad /= da;

    Vector<double, 13> out;
    out.head(12) = ab_grad * cd + cd_grad * ab;
    out(12) = - (a * b_wrt_ma1 + b * a_wrt_ma0) * cd / db - (c * d_wrt_mb1 + d * c_wrt_mb0) * ab / da;

    return std::make_pair(ab * cd, out);
}

/// @brief Compute the hessian of the mollifier function wrt. 4 edge points and the distance squared
std::tuple<double, Vector<double, 13>, Eigen::Matrix<double, 13, 13>>
edge_edge_mollifier_hessian(
    const Vector3<double> &ea0, const Vector3<double> &ea1,
    const Vector3<double> &eb0, const Vector3<double> &eb1, 
    const double &dist_sqr)
{
    DiffScalarBase::setVariableCount(13);
    using T = ADHessian<13>;
    Vector<double, 13> input;
    input << ea0, ea1, eb0, eb1, dist_sqr;
    Vector<T, 13> input_ad = slice_positions<T, 13, 1>(input);

    const double da0 = point_edge_distance(ea0, eb0, eb1);
    const double da1 = point_edge_distance(ea1, eb0, eb1);
    const double db0 = point_edge_distance(eb0, ea0, ea1);
    const double db1 = point_edge_distance(eb1, ea0, ea1);

    Vector<double,13> da0_wrt_x_1, da1_wrt_x_1, db0_wrt_x_1, db1_wrt_x_1;
    {
        Vector9d tmp = point_edge_distance_gradient(ea0, eb0, eb1);
        da0_wrt_x_1 << tmp.head(3), Eigen::Vector3d::Zero(), tmp.tail(6), 0.;

        tmp = point_edge_distance_gradient(ea1, eb0, eb1);
        da1_wrt_x_1 << Eigen::Vector3d::Zero(), tmp, 0.;

        tmp = point_edge_distance_gradient(eb0, ea0, ea1);
        db0_wrt_x_1 << tmp.tail(6), tmp.head(3), Eigen::Vector3d::Zero(), 0.;

        tmp = point_edge_distance_gradient(eb1, ea0, ea1);
        db1_wrt_x_1 << tmp.tail(6), Eigen::Vector3d::Zero(), tmp.head(3), 0.;
    }

    Eigen::Matrix<double, 13, 13> da0_wrt_x_2, da1_wrt_x_2, db0_wrt_x_2, db1_wrt_x_2;
    {
        Matrix9d tmp = point_edge_distance_hessian(ea0, eb0, eb1);
        da0_wrt_x_2 << tmp.topLeftCorner(3,3),  Eigen::Matrix<double, 3,3>::Zero() , tmp.topRightCorner(3,6), Eigen::Matrix<double, 3,1>::Zero(),
                                                Eigen::Matrix<double, 3,13>::Zero(),
                     tmp.bottomLeftCorner(6,3), Eigen::Matrix<double, 6,3>::Zero() , tmp.bottomRightCorner(6,6), Eigen::Matrix<double, 6,1>::Zero(),
                                                Eigen::Matrix<double, 1,13>::Zero();

        tmp = point_edge_distance_hessian(ea1, eb0, eb1);
        da1_wrt_x_2 << Eigen::Matrix<double, 3,13>::Zero(),
                       Eigen::Matrix<double, 9,3>::Zero(), tmp, Eigen::Matrix<double, 9, 1>::Zero(),
                       Eigen::Matrix<double, 1,13>::Zero();

        tmp = point_edge_distance_hessian(eb0, ea0, ea1);
        db0_wrt_x_2 << tmp.bottomRightCorner(6, 6), tmp.bottomLeftCorner(6, 3), Eigen::Matrix<double, 6, 4>::Zero(),
                          tmp.topRightCorner(3, 6), tmp.topLeftCorner(3, 3)   , Eigen::Matrix<double, 3, 4>::Zero(), 
                                                    Eigen::Matrix<double, 4, 13>::Zero();

        tmp = point_edge_distance_hessian(eb1, ea0, ea1);
        db1_wrt_x_2 << tmp.bottomRightCorner(6, 6), Eigen::Matrix<double, 6, 3>::Zero(), tmp.bottomLeftCorner(6, 3), Eigen::Matrix<double, 6, 1>::Zero(),
                                                    Eigen::Matrix<double, 3, 13>::Zero(),
                       tmp.topRightCorner(3, 6)   , Eigen::Matrix<double, 3, 3>::Zero(), tmp.topLeftCorner(3, 3), Eigen::Matrix<double, 3, 1>::Zero(),
                                                    Eigen::Matrix<double, 1, 13>::Zero();
    }

    T laAD = (input_ad.head(3) - input_ad.segment(3,3)).squaredNorm() * mollifier_threshold_eps;
    T lbAD = (input_ad.segment(6,3) - input_ad.segment(9,3)).squaredNorm() * mollifier_threshold_eps;
    T aAD = Math<T>::mollifier((T(da0, da0_wrt_x_1, da0_wrt_x_2) - input_ad(12)) / lbAD);
    T bAD = Math<T>::mollifier((T(da1, da1_wrt_x_1, da1_wrt_x_2) - input_ad(12)) / lbAD);
    T cAD = Math<T>::mollifier((T(db0, db0_wrt_x_1, db0_wrt_x_2) - input_ad(12)) / laAD);
    T dAD = Math<T>::mollifier((T(db1, db1_wrt_x_1, db1_wrt_x_2) - input_ad(12)) / laAD);

    T outAD = aAD * bAD * cAD * dAD;
    return std::make_tuple(outAD.getValue(), outAD.getGradient(), outAD.getHessian());
}

}
