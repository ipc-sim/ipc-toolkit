#pragma once

#include <ipc/distance/distance_type.hpp>
#include "math.hpp"

namespace ipc {
    template <typename scalar>
    scalar point_point_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& a,
        const Eigen::Ref<const Vector3<scalar>>& b)
    {
        return (a - b).squaredNorm();
    }

    template <typename scalar>
    scalar line_line_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1)
    {
        const Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
        const scalar line_to_line = (eb0 - ea0).dot(normal);
        return line_to_line * line_to_line / normal.squaredNorm();
    }

    template <typename scalar>
    scalar point_line_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1)
    {
        return (e0 - p).cross(e1 - p).squaredNorm() / (e1 - e0).squaredNorm();
    }

    template <typename scalar>
    scalar point_plane_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2)
    {
        const Vector3<scalar> normal = (f2 - f0).cross(f1 - f0);
        return intpow<2>(normal.dot(p - f0)) / normal.squaredNorm();
    }

    template <typename scalar>
    scalar point_triangle_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& t0,
        const Eigen::Ref<const Vector3<scalar>>& t1,
        const Eigen::Ref<const Vector3<scalar>>& t2,
        PointTriangleDistanceType dtype)
    {
        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            return point_point_sqr_distance<scalar>(p, t0);

        case PointTriangleDistanceType::P_T1:
            return point_point_sqr_distance<scalar>(p, t1);

        case PointTriangleDistanceType::P_T2:
            return point_point_sqr_distance<scalar>(p, t2);

        case PointTriangleDistanceType::P_E0:
            return point_line_sqr_distance<scalar>(p, t0, t1);

        case PointTriangleDistanceType::P_E1:
            return point_line_sqr_distance<scalar>(p, t1, t2);

        case PointTriangleDistanceType::P_E2:
            return point_line_sqr_distance<scalar>(p, t2, t0);

        case PointTriangleDistanceType::P_T:
            return point_plane_sqr_distance<scalar>(p, t0, t1, t2);

        default:
            throw std::invalid_argument(
                "Invalid distance type for point-triangle distance!");
        }
    }

    template <typename scalar>
    scalar point_edge_sqr_distance(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType &dtype)
    {
        switch (dtype)
        {
        case PointEdgeDistanceType::P_E:
            return point_line_sqr_distance<scalar>(p, e0, e1);
        case PointEdgeDistanceType::P_E0:
            return point_point_sqr_distance<scalar>(p, e0);
        case PointEdgeDistanceType::P_E1:
            return point_point_sqr_distance<scalar>(p, e1);
        case PointEdgeDistanceType::AUTO:
        default:
            VectorMax3<scalar> t = e1 - e0;
            const scalar len = t.norm();
            t = t / len;

            const VectorMax3<scalar> pos = p - e0;
            const scalar s = pos.dot(t) / len;
            return (pos - (L_ns(s) * len) * t).squaredNorm();
        }
    }

    template <typename scalar>
    scalar edge_edge_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        EdgeEdgeDistanceType dtype)
    {
        switch (dtype) {
        case EdgeEdgeDistanceType::EA0_EB0:
            return point_point_sqr_distance<scalar>(ea0, eb0);

        case EdgeEdgeDistanceType::EA0_EB1:
            return point_point_sqr_distance<scalar>(ea0, eb1);

        case EdgeEdgeDistanceType::EA1_EB0:
            return point_point_sqr_distance<scalar>(ea1, eb0);

        case EdgeEdgeDistanceType::EA1_EB1:
            return point_point_sqr_distance<scalar>(ea1, eb1);

        case EdgeEdgeDistanceType::EA_EB0:
            return point_line_sqr_distance<scalar>(eb0, ea0, ea1);

        case EdgeEdgeDistanceType::EA_EB1:
            return point_line_sqr_distance<scalar>(eb1, ea0, ea1);

        case EdgeEdgeDistanceType::EA0_EB:
            return point_line_sqr_distance<scalar>(ea0, eb0, eb1);

        case EdgeEdgeDistanceType::EA1_EB:
            return point_line_sqr_distance<scalar>(ea1, eb0, eb1);

        case EdgeEdgeDistanceType::EA_EB:
            return line_line_sqr_distance<scalar>(ea0, ea1, eb0, eb1);

        default:
            throw std::invalid_argument(
                "Invalid distance type for edge-edge distance!");
        }
    }

    template <typename scalar>
    VectorMax3<scalar> point_line_closest_point_direction(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1)
    {
        const VectorMax3<scalar> d = p - e0;
        const VectorMax3<scalar> t = e1 - e0;
        return d - (d.dot(t) / t.squaredNorm()) * t;
    }

    template <typename scalar>
    VectorMax3<scalar> point_edge_closest_point_direction(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType &dtype = PointEdgeDistanceType::AUTO)
    {
        switch (dtype)
        {
        case PointEdgeDistanceType::P_E:
            return point_line_closest_point_direction<scalar>(p, e0, e1);
        case PointEdgeDistanceType::P_E0:
            return p - e0;
        case PointEdgeDistanceType::P_E1:
            return p - e1;
        case PointEdgeDistanceType::AUTO:
        default:
            VectorMax3<scalar> t = e1 - e0;
            const scalar len = t.norm();
            t = t / len;

            const VectorMax3<scalar> pos = p - e0;
            const scalar s = pos.dot(t) / len;
            return pos - (L_ns(s) * len) * t;
        }
    }

    template <typename scalar>
    Vector3<scalar> line_line_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1)
    {
        const Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
        return ((eb0 - ea0).dot(normal) / normal.squaredNorm()) * normal;
    }

    /// @brief Computes the direction of the closest point pair
    /// @param ea0 Vertex 0 of edge 0
    /// @param ea1 Vertex 1 of edge 0
    /// @param eb0 Vertex 0 of edge 1
    /// @param eb1 Vertex 1 of edge 1
    /// @param dtype Edge-edge distance type
    /// @return Difference of the pair of closest point, pointing from edge 0 to edge 1
    template <typename scalar>
    Vector3<scalar> edge_edge_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        EdgeEdgeDistanceType dtype)
    {
        switch (dtype) {
        case EdgeEdgeDistanceType::EA0_EB0:
            return (eb0 - ea0);

        case EdgeEdgeDistanceType::EA0_EB1:
            return (eb1 - ea0);

        case EdgeEdgeDistanceType::EA1_EB0:
            return (eb0 - ea1);

        case EdgeEdgeDistanceType::EA1_EB1:
            return (eb1 - ea1);

        case EdgeEdgeDistanceType::EA_EB0:
            return point_line_closest_point_direction<scalar>(eb0, ea0, ea1);

        case EdgeEdgeDistanceType::EA_EB1:
            return point_line_closest_point_direction<scalar>(eb1, ea0, ea1);

        case EdgeEdgeDistanceType::EA0_EB:
            return -point_line_closest_point_direction<scalar>(ea0, eb0, eb1);

        case EdgeEdgeDistanceType::EA1_EB:
            return -point_line_closest_point_direction<scalar>(ea1, eb0, eb1);

        case EdgeEdgeDistanceType::EA_EB:
            return line_line_closest_point_direction<scalar>(ea0, ea1, eb0, eb1);

        default:
            throw std::invalid_argument(
                "Invalid distance type for edge-edge distance!");
        }
    }

    template <typename scalar>
    Vector3<scalar> point_plane_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2)
    {
        const Vector3<scalar> normal = (f2 - f0).cross(f1 - f0);
        return (normal.dot(p - f0) / normal.squaredNorm()) * normal;
    }

    template <typename scalar>
    Vector3<scalar> point_triangle_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& t0,
        const Eigen::Ref<const Vector3<scalar>>& t1,
        const Eigen::Ref<const Vector3<scalar>>& t2,
        PointTriangleDistanceType dtype)
    {
        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            return p - t0;

        case PointTriangleDistanceType::P_T1:
            return p - t1;

        case PointTriangleDistanceType::P_T2:
            return p - t2;

        case PointTriangleDistanceType::P_E0:
            return point_line_closest_point_direction<scalar>(p, t0, t1);

        case PointTriangleDistanceType::P_E1:
            return point_line_closest_point_direction<scalar>(p, t1, t2);

        case PointTriangleDistanceType::P_E2:
            return point_line_closest_point_direction<scalar>(p, t2, t0);

        case PointTriangleDistanceType::P_T:
            return point_plane_closest_point_direction<scalar>(p, t0, t1, t2);

        default:
            throw std::invalid_argument(
                "Invalid distance type for point-triangle distance!");
        }
    }

    template <typename scalar>
    scalar edge_mollifier(const VectorMax3<scalar> &p, const VectorMax3<scalar> &e0, const VectorMax3<scalar> &e1, const scalar &dist_sqr)
    {
        const scalar denominator = dist_sqr * mollifier_threshold_eps;
        return mollifier<scalar>(((p - e0).squaredNorm() - dist_sqr) / denominator) *
            mollifier<scalar>(((p - e1).squaredNorm() - dist_sqr) / denominator);
    }

    template <typename scalar>
    scalar edge_edge_mollifier(
        const Vector3<scalar> &ea0, const Vector3<scalar> &ea1,
        const Vector3<scalar> &eb0, const Vector3<scalar> &eb1, 
        const scalar &dist_sqr)
    {
        const scalar denominator = dist_sqr * mollifier_threshold_eps;
        scalar a = mollifier<scalar>((point_edge_sqr_distance<scalar>(ea0, eb0, eb1, PointEdgeDistanceType::AUTO) - dist_sqr) / denominator);
        scalar b = mollifier<scalar>((point_edge_sqr_distance<scalar>(ea1, eb0, eb1, PointEdgeDistanceType::AUTO) - dist_sqr) / denominator);
        scalar c = mollifier<scalar>((point_edge_sqr_distance<scalar>(eb0, ea0, ea1, PointEdgeDistanceType::AUTO) - dist_sqr) / denominator);
        scalar d = mollifier<scalar>((point_edge_sqr_distance<scalar>(eb1, ea0, ea1, PointEdgeDistanceType::AUTO) - dist_sqr) / denominator);
        
        return a * b * c * d;
    }

    template <typename scalar>
    scalar triangle_mollifier(
        const VectorMax3<scalar> &p, 
        const VectorMax3<scalar> &e0, 
        const VectorMax3<scalar> &e1,
        const VectorMax3<scalar> &e2,
        const scalar &dist_sqr)
    {
        const scalar denominator = dist_sqr * mollifier_threshold_eps;
        return mollifier<scalar>(point_edge_sqr_distance<scalar>(p, e0, e1, PointEdgeDistanceType::AUTO) / denominator) *
            mollifier<scalar>(point_edge_sqr_distance<scalar>(p, e2, e1, PointEdgeDistanceType::AUTO) / denominator) *
            mollifier<scalar>(point_edge_sqr_distance<scalar>(p, e0, e2, PointEdgeDistanceType::AUTO) / denominator);
    }
}