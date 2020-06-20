#pragma once

#include <Math/Distance/DISTANCE_TYPE.h>
#include <Math/Distance/POINT_TRIANGLE.h>
#include <Math/Distance/POINT_EDGE.h>
#include <Math/Distance/POINT_POINT.h>
#include <Math/Distance/EDGE_EDGE.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
namespace JGSL {
    
template <class T>
void Point_Triangle_Distance_Unclassified(
    const Eigen::Matrix<T, 3, 1>& p,
    const Eigen::Matrix<T, 3, 1>& t0,
    const Eigen::Matrix<T, 3, 1>& t1,
    const Eigen::Matrix<T, 3, 1>& t2,
    T& dist2)
{
    switch (Point_Triangle_Distance_Type(p, t0, t1, t2)) {
    case 0: {
        Point_Point_Distance(p, t0, dist2);
        break;
    }

    case 1: {
        Point_Point_Distance(p, t1, dist2);
        break;
    }

    case 2: {
        Point_Point_Distance(p, t2, dist2);
        break;
    }

    case 3: {
        Point_Edge_Distance(p, t0, t1, dist2);
        break;
    }

    case 4: {
        Point_Edge_Distance(p, t1, t2, dist2);
        break;
    }

    case 5: {
        Point_Edge_Distance(p, t2, t0, dist2);
        break;
    }

    case 6: {
        Point_Triangle_Distance(p, t0, t1, t2, dist2);
        break;
    }

    default:
        break;
    }
}

template <class T>
void Edge_Edge_Distance_Unclassified(
    const Eigen::Matrix<T, 3, 1>& ea0,
    const Eigen::Matrix<T, 3, 1>& ea1,
    const Eigen::Matrix<T, 3, 1>& eb0,
    const Eigen::Matrix<T, 3, 1>& eb1,
    T& dist2)
{
    switch (Edge_Edge_Distance_Type(ea0, ea1, eb0, eb1)) {
    case 0: {
        Point_Point_Distance(ea0, eb0, dist2);
        break;
    }

    case 1: {
        Point_Point_Distance(ea0, eb1, dist2);
        break;
    }

    case 2: {
        Point_Edge_Distance(ea0, eb0, eb1, dist2);
        break;
    }

    case 3: {
        Point_Point_Distance(ea1, eb0, dist2);
        break;
    }

    case 4: {
        Point_Point_Distance(ea1, eb1, dist2);
        break;
    }

    case 5: {
        Point_Edge_Distance(ea1, eb0, eb1, dist2);
        break;
    }

    case 6: {
        Point_Edge_Distance(eb0, ea0, ea1, dist2);
        break;
    }

    case 7: {
        Point_Edge_Distance(eb1, ea0, ea1, dist2);
        break;
    }

    case 8: {
        Edge_Edge_Distance(ea0, ea1, eb0, eb1, dist2);
        break;
    }

    default:
        break;
    }
}

// http://geomalgorithms.com/a02-_lines.html
template <class T>
void Point_Edge_Distance_Unclassified(
    const Eigen::Matrix<T, 3, 1>& p,
    const Eigen::Matrix<T, 3, 1>& e0,
    const Eigen::Matrix<T, 3, 1>& e1,
    T& dist2)
{
    Eigen::Matrix<T, 3, 1> v = e1 - e0;
    Eigen::Matrix<T, 3, 1> w = p - e0;

    double c1 = w.dot(v);
    if (c1 <= 0.0) {
        Point_Point_Distance(p, e0, dist2);
    }
    else {
        double c2 = v.squaredNorm();
        if (c2 <= c1) {
            Point_Point_Distance(p, e1, dist2);
        }
        else {
            double b = c1 / c2;
            Point_Point_Distance(p, std::move(Eigen::Matrix<T, 3, 1>(e0 + b * v)), dist2);
        }
    }
}

}