#pragma once

#include <Math/Distance/DISTANCE_TYPE.h>
#include <Math/Distance/DISTANCE_UNCLASSIFIED.h>

#include <Math/Distance/EVCTCD/CTCD.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace JGSL {

template <class T>
bool Point_Edge_CD_Broadphase(
    const Eigen::Matrix<T, 2, 1>& x0, 
    const Eigen::Matrix<T, 2, 1>& x1, 
    const Eigen::Matrix<T, 2, 1>& x2,
    T dist)
{
    const Eigen::Array<T, 2, 1> max_e = x1.array().max(x2.array());
    const Eigen::Array<T, 2, 1> min_e = x1.array().min(x2.array());
    if ((x0.array() - max_e > dist).any() || (min_e - x0.array() > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Point_Edge_CCD_Broadphase(
    const Eigen::Matrix<T, 2, 1>& p, 
    const Eigen::Matrix<T, 2, 1>& e0, 
    const Eigen::Matrix<T, 2, 1>& e1,
    const Eigen::Matrix<T, 2, 1>& dp, 
    const Eigen::Matrix<T, 2, 1>& de0, 
    const Eigen::Matrix<T, 2, 1>& de1,
    T dist)
{
    const Eigen::Array<T, 2, 1> max_p = p.array().max((p + dp).array());
    const Eigen::Array<T, 2, 1> min_p = p.array().min((p + dp).array());
    const Eigen::Array<T, 2, 1> max_e = e0.array().max(e1.array()).
        max((e0 + de0).array()).max((e1 + de1).array());
    const Eigen::Array<T, 2, 1> min_e = e0.array().min(e1.array()).
        min((e0 + de0).array()).min((e1 + de1).array());
    if ((min_p - max_e > dist).any() || (min_e - max_p > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Point_Edge_CCD(const Eigen::Matrix<T, 2, 1>& x0, 
    const Eigen::Matrix<T, 2, 1>& x1, 
    const Eigen::Matrix<T, 2, 1>& x2,
    const Eigen::Matrix<T, 2, 1>& d0, 
    const Eigen::Matrix<T, 2, 1>& d1, 
    const Eigen::Matrix<T, 2, 1>& d2,
    T eta, T& toc)
{
    T a = d0[0] * (d2[1] - d1[1]) + d0[1] * (d1[0] - d2[0]) + d2[0] * d1[1] - d2[1] * d1[0];
    T b = x0[0] * (d2[1] - d1[1]) + d0[0] * (x2[1] - x1[1]) + 
        d0[1] * (x1[0] - x2[0]) + x0[1] * (d1[0] - d2[0]) +
        d1[1] * x2[0] + d2[0] * x1[1] - d1[0] * x2[1] - d2[1] * x1[0];
    T c = x0[0] * (x2[1] - x1[1]) + x0[1] * (x1[0] - x2[0]) + x2[0] * x1[1] - x2[1] * x1[0];

    T roots[2];
    int rootAmt = 0;
    if (a == 0) {
        if (b == 0) {
            // parallel motion, only need to handle colinear case
            if (c == 0) {
                // colinear PP CCD
                if ((x0 - x1).dot(d0 - d1) < 0) {
                    roots[rootAmt] = std::sqrt((x0 - x1).squaredNorm() / (d0 - d1).squaredNorm());
                    if (roots[rootAmt] > 0 && roots[rootAmt] <= 1) {
                        ++rootAmt;
                    }
                }
                if ((x0 - x2).dot(d0 - d2) < 0) {
                    roots[rootAmt] = std::sqrt((x0 - x2).squaredNorm() / (d0 - d2).squaredNorm());
                    if (roots[rootAmt] > 0 && roots[rootAmt] <= 1) {
                        ++rootAmt;
                    }
                }

                if (rootAmt == 2) {
                    toc = std::min(roots[0], roots[1]) * (1 - eta);
                    return true;
                }
                else if (rootAmt == 1) {
                    toc = roots[0] * (1 - eta);
                    return true;
                }
                else {
                    return false;
                }
            }
        }
        else {
            rootAmt = 1;
            roots[0] = -c / b;
        }
    }
    else {
        T delta = b * b - 4 * a * c;
        if (delta == 0) {
            rootAmt = 1;
            roots[0] = -b / (2 * a);
        }
        else if (delta > 0) {
            rootAmt = 2;
            // accurate expression differs in b's sign
            if (b > 0) {
                roots[0] = (-b - std::sqrt(delta)) / (2 * a);
                roots[1] = 2 * c / (-b - std::sqrt(delta));
            }
            else {
                roots[0] = 2 * c / (-b + std::sqrt(delta));
                roots[1] = (-b + std::sqrt(delta)) / (2 * a);
            }

            if (roots[0] > roots[1]) {
                std::swap(roots[0], roots[1]);
            }
        }
    }

    for (int i = 0; i < rootAmt; ++i) {
        if (roots[i] > 0 && roots[i] <= 1) {
            // check overlap
            T ratio;
            if(Point_Edge_Distance_Type(Eigen::Matrix<T, 2, 1>(x0 + roots[i] * d0), 
                Eigen::Matrix<T, 2, 1>(x1 + roots[i] * d1), 
                Eigen::Matrix<T, 2, 1>(x2 + roots[i] * d2), ratio) == 2) {
                toc = roots[i] * (1 - eta); //TODO: distance eta
                return true;
            }
        }
    }

    return false;
}

template <class T>
bool Point_Triangle_CD_Broadphase(
    const Eigen::Matrix<T, 3, 1>& p, 
    const Eigen::Matrix<T, 3, 1>& t0, 
    const Eigen::Matrix<T, 3, 1>& t1,
    const Eigen::Matrix<T, 3, 1>& t2,
    T dist)
{
    const Eigen::Array<T, 3, 1> max_tri = t0.array().max(t1.array()).max(t2.array());
    const Eigen::Array<T, 3, 1> min_tri = t0.array().min(t1.array()).min(t2.array());
    if ((p.array() - max_tri > dist).any() || (min_tri - p.array() > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Edge_Edge_CD_Broadphase(
    const Eigen::Matrix<T, 3, 1>& ea0, 
    const Eigen::Matrix<T, 3, 1>& ea1, 
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    T dist)
{
    const Eigen::Array<T, 3, 1> max_a = ea0.array().max(ea1.array());
    const Eigen::Array<T, 3, 1> min_a = ea0.array().min(ea1.array());
    const Eigen::Array<T, 3, 1> max_b = eb0.array().max(eb1.array());
    const Eigen::Array<T, 3, 1> min_b = eb0.array().min(eb1.array());
    if ((min_a - max_b > dist).any() || (min_b - max_a > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Point_Triangle_CCD_Broadphase(
    const Eigen::Matrix<T, 3, 1>& p, 
    const Eigen::Matrix<T, 3, 1>& t0, 
    const Eigen::Matrix<T, 3, 1>& t1,
    const Eigen::Matrix<T, 3, 1>& t2,
    const Eigen::Matrix<T, 3, 1>& dp, 
    const Eigen::Matrix<T, 3, 1>& dt0, 
    const Eigen::Matrix<T, 3, 1>& dt1,
    const Eigen::Matrix<T, 3, 1>& dt2,
    T dist)
{
    const Eigen::Array<T, 3, 1> max_p = p.array().max((p + dp).array());
    const Eigen::Array<T, 3, 1> min_p = p.array().min((p + dp).array());
    const Eigen::Array<T, 3, 1> max_tri = t0.array().max(t1.array()).max(t2.array()).
        max((t0 + dt0).array()).max((t1 + dt1).array()).max((t2 + dt2).array());
    const Eigen::Array<T, 3, 1> min_tri = t0.array().min(t1.array()).min(t2.array()).
        min((t0 + dt0).array()).min((t1 + dt1).array()).min((t2 + dt2).array());
    if ((min_p - max_tri > dist).any() || (min_tri - max_p > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Edge_Edge_CCD_Broadphase(
    const Eigen::Matrix<T, 3, 1>& ea0, 
    const Eigen::Matrix<T, 3, 1>& ea1, 
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    const Eigen::Matrix<T, 3, 1>& dea0, 
    const Eigen::Matrix<T, 3, 1>& dea1, 
    const Eigen::Matrix<T, 3, 1>& deb0,
    const Eigen::Matrix<T, 3, 1>& deb1,
    T dist)
{
    const Eigen::Array<T, 3, 1> max_a = ea0.array().max(ea1.array()).max((ea0 + dea0).array()).max((ea1 + dea1).array());
    const Eigen::Array<T, 3, 1> min_a = ea0.array().min(ea1.array()).min((ea0 + dea0).array()).min((ea1 + dea1).array());
    const Eigen::Array<T, 3, 1> max_b = eb0.array().max(eb1.array()).max((eb0 + deb0).array()).max((eb1 + deb1).array());
    const Eigen::Array<T, 3, 1> min_b = eb0.array().min(eb1.array()).min((eb0 + deb0).array()).min((eb1 + deb1).array());
    if ((min_a - max_b > dist).any() || (min_b - max_a > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Point_Edge_CCD_Broadphase(
    const Eigen::Matrix<T, 3, 1>& p, 
    const Eigen::Matrix<T, 3, 1>& e0, 
    const Eigen::Matrix<T, 3, 1>& e1,
    const Eigen::Matrix<T, 3, 1>& dp, 
    const Eigen::Matrix<T, 3, 1>& de0, 
    const Eigen::Matrix<T, 3, 1>& de1,
    T dist)
{
    const Eigen::Array<T, 3, 1> max_p = p.array().max((p + dp).array());
    const Eigen::Array<T, 3, 1> min_p = p.array().min((p + dp).array());
    const Eigen::Array<T, 3, 1> max_e = e0.array().max(e1.array()).max((e0 + de0).array()).max((e1 + de1).array());
    const Eigen::Array<T, 3, 1> min_e = e0.array().min(e1.array()).min((e0 + de0).array()).min((e1 + de1).array());
    if ((min_p - max_e > dist).any() || (min_e - max_p > dist).any()) {
        return false;
    }
    else {
        return true;
    }
}

template <class T>
bool Point_Triangle_CCD(
    const Eigen::Matrix<T, 3, 1>& p, 
    const Eigen::Matrix<T, 3, 1>& t0, 
    const Eigen::Matrix<T, 3, 1>& t1,
    const Eigen::Matrix<T, 3, 1>& t2,
    const Eigen::Matrix<T, 3, 1>& dp, 
    const Eigen::Matrix<T, 3, 1>& dt0, 
    const Eigen::Matrix<T, 3, 1>& dt1,
    const Eigen::Matrix<T, 3, 1>& dt2,
    T eta, T& toc)
{
    T dist2_cur;
    Point_Triangle_Distance_Unclassified(p, t0, t1, t2, dist2_cur);
    if (CTCD::vertexFaceCTCD(p, t0, t1, t2,
        p + dp, t0 + dt0, t1 + dt1, t2 + dt2,
        eta * std::sqrt(dist2_cur), toc))
    {
        if (toc < 1.0e-6) {
            std::cout << "PT CCD tiny!" << std::endl;
            if (CTCD::vertexFaceCTCD(p, t0, t1, t2,
                    p + dp, t0 + dt0, t1 + dt1, t2 + dt2,
                    0, toc))
            {
                toc *= (1.0 - eta);
                return true;
            }
            else {
                return false;
            }
        }
        return true;
    }
    return false;
}

template <class T>
bool Edge_Edge_CCD(
    const Eigen::Matrix<T, 3, 1>& ea0, 
    const Eigen::Matrix<T, 3, 1>& ea1, 
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    const Eigen::Matrix<T, 3, 1>& dea0, 
    const Eigen::Matrix<T, 3, 1>& dea1, 
    const Eigen::Matrix<T, 3, 1>& deb0,
    const Eigen::Matrix<T, 3, 1>& deb1,
    T eta, T& toc)
{
    T dist2_cur;
    Edge_Edge_Distance_Unclassified(ea0, ea1, eb0, eb1, dist2_cur);
    if (CTCD::edgeEdgeCTCD(ea0, ea1, eb0, eb1,
        ea0 + dea0, ea1 + dea1, eb0 + deb0, eb1 + deb1,
        eta * std::sqrt(dist2_cur), toc))
    {
        if (toc < 1.0e-6) {
            std::cout << "EE CCD tiny!" << std::endl;
            if (CTCD::edgeEdgeCTCD(ea0, ea1, eb0, eb1,
                    ea0 + dea0, ea1 + dea1, eb0 + deb0, eb1 + deb1,
                    0, toc))
            {
                toc *= (1.0 - eta);
                return true;
            }
            else {
                return false;
            }
        }
        return true;
    }
    return false;
}

template <class T>
bool Point_Edge_CCD(
    const Eigen::Matrix<T, 3, 1>& p, 
    const Eigen::Matrix<T, 3, 1>& e0,
    const Eigen::Matrix<T, 3, 1>& e1,
    const Eigen::Matrix<T, 3, 1>& dp, 
    const Eigen::Matrix<T, 3, 1>& de0,
    const Eigen::Matrix<T, 3, 1>& de1,
    T eta, T& toc)
{
    T dist2_cur;
    Point_Edge_Distance_Unclassified(p, e0, e1, dist2_cur);
    if (CTCD::vertexEdgeCTCD(p, e0, e1, p + dp, e0 + de0, e1 + de1,
        eta * std::sqrt(dist2_cur), toc))
    {
        if (toc < 1.0e-6) {
            std::cout << "PE CCD tiny!" << std::endl;
            if (CTCD::vertexEdgeCTCD(p, e0, e1, 
                    p + dp, e0 + de0, e1 + de1,
                    0, toc))
            {
                toc *= (1.0 - eta);
                return true;
            }
            else {
                return false;
            }
        }
        return true;
    }
    return false;
}

}