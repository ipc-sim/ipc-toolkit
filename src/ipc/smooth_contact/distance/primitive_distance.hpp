#pragma once
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/distance/distance_type.hpp>

#include <ipc/smooth_contact/primitives/point.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>

namespace ipc {
    template <typename PrimitiveA, typename PrimitiveB>
    struct PrimitiveDistType {};

    template <>
    struct PrimitiveDistType<Point3, Point3>
    { using type = PointPointDistanceType; };

    template <>
    struct PrimitiveDistType<Edge3, Point3>
    { using type = PointEdgeDistanceType; };

    template <>
    struct PrimitiveDistType<Face, Point3>
    { using type = PointTriangleDistanceType; };

    template <>
    struct PrimitiveDistType<Edge3, Edge3>
    { using type = EdgeEdgeDistanceType; };

    template <typename PrimitiveA, typename PrimitiveB>
    class PrimitiveDistance
    {
    public:
        static typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type compute_distance_type(const Vector<double, -1, 12>& x);
        
        static double compute_distance(const CollisionMesh &mesh, const Eigen::MatrixXd &V, const long &a, const long &b, typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);

        // points from primitiveA to primitiveB
        static VectorMax3d compute_closest_direction(const CollisionMesh &mesh, const Eigen::MatrixXd &V, const long &a, const long &b, typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
    };

    template <typename PrimitiveA, typename PrimitiveB, typename T>
    class PrimitiveDistanceTemplate
    {
    public:
        static VectorMax3<T> compute_closest_direction(const Vector<T, -1, 12>& x, typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
        // mollifier only depends on at most 4 points
        static T mollifier(const Vector<T, -1, 12>& x, const T& dist_sqr);
    };
    
}