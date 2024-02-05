#include "primitive_distance.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <ipc/smooth_contact/distance/edge_edge.hpp>
#include <ipc/smooth_contact/distance/point_edge.hpp>
#include <ipc/smooth_contact/distance/point_face.hpp>

namespace ipc {

template <>
typename PrimitiveDistType<Face, Point3>::type
PrimitiveDistance<Face, Point3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return point_triangle_distance_type(
        x.tail(3) /* point */, x.head(3), x.segment(3, 3),
        x.segment(6, 3) /* face */);
}

template <>
typename PrimitiveDistType<Edge3, Edge3>::type
PrimitiveDistance<Edge3, Edge3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return edge_edge_distance_type(
        x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
        x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */);
}

template <>
typename PrimitiveDistType<Edge3, Point3>::type
PrimitiveDistance<Edge3, Point3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return point_edge_distance_type(
        x.tail(3) /* point */, x.head(3) /* edge */,
        x.segment(3, 3) /* edge */);
}

template <>
typename PrimitiveDistType<Point3, Point3>::type
PrimitiveDistance<Point3, Point3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return PointPointDistanceType::AUTO;
}

template <>
double PrimitiveDistance<Face, Point3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    return point_triangle_distance(
        V.row(b), V.row(mesh.faces()(a, 0)), V.row(mesh.faces()(a, 1)),
        V.row(mesh.faces()(a, 2)), dtype);
}

template <>
double PrimitiveDistance<Edge3, Edge3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return edge_edge_distance(
        V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)),
        V.row(mesh.edges()(b, 0)), V.row(mesh.edges()(b, 1)), dtype);
}

template <>
double PrimitiveDistance<Edge3, Point3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    return point_edge_distance(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
double PrimitiveDistance<Point3, Point3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    return point_point_distance(V.row(a), V.row(b));
}

template <>
Vector<double, PrimitiveDistance<Face, Point3>::dim>
PrimitiveDistance<Face, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    return point_triangle_closest_point_direction<double>(
        V.row(b), V.row(mesh.faces()(a, 0)), V.row(mesh.faces()(a, 1)),
        V.row(mesh.faces()(a, 2)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge3, Edge3>::dim>
PrimitiveDistance<Edge3, Edge3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return edge_edge_closest_point_direction<double>(
        V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)),
        V.row(mesh.edges()(b, 0)), V.row(mesh.edges()(b, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge3, Point3>::dim>
PrimitiveDistance<Edge3, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    return PointEdgeDistance<double, 3>::point_edge_closest_point_direction(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Point3, Point3>::dim>
PrimitiveDistance<Point3, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    return V.row(b) - V.row(a);
}

template <typename T> class PrimitiveDistanceTemplate<Face, Point3, T> {
    static_assert(
        Face::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Face::dim;
    constexpr static int n_core_dofs =
        Face::n_core_points * Face::dim + Point3::n_core_points * Point3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Face, Point3>::type dtype)
    {
        return point_triangle_closest_point_direction<T>(
            x.tail(3) /* point */, x.head(3), x.segment(3, 3),
            x.segment(6, 3) /* face */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_face_mollifier<T>(
            x.tail(3) /* point */, x.head(3), x.segment(3, 3),
            x.segment(6, 3) /* face */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Edge3, T> {
    static_assert(
        Edge3::dim == Edge3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Edge3::dim;
    constexpr static int n_core_dofs =
        Edge3::n_core_points * Edge3::dim + Edge3::n_core_points * Edge3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Edge3>::type dtype)
    {
        return edge_edge_closest_point_direction<T>(
            x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
            x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        std::array<HEAVISIDE_TYPE, 4> types;
        types.fill(HEAVISIDE_TYPE::VARIANT);
        return edge_edge_mollifier<T>(
            x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
            x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */, types,
            dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Point3, T> {
    static_assert(
        Edge3::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Edge3::dim;
    constexpr static int n_core_dofs =
        Edge3::n_core_points * Edge3::dim + Point3::n_core_points * Point3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Point3>::type dtype)
    {
        return PointEdgeDistance<T, dim>::point_edge_closest_point_direction(
            x.tail(3) /* point */, x.head(3) /* edge */,
            x.segment(3, 3) /* edge */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_edge_mollifier<T>(
            x.tail(3) /* point */, x.segment(3, 3) /* edge */,
            x.head(3) /* edge */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Point3, Point3, T> {
    static_assert(
        Point3::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Point3::dim;
    constexpr static int n_core_dofs = Point3::n_core_points * Point3::dim
        + Point3::n_core_points * Point3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Point3, Point3>::type dtype)
    {
        return x.tail(3) - x.head(3);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return T(1.);
    }
};

template class PrimitiveDistance<Edge3, Point3>;
template class PrimitiveDistance<Edge3, Edge3>;
template class PrimitiveDistance<Point3, Point3>;
template class PrimitiveDistance<Face, Point3>;

template class PrimitiveDistanceTemplate<Edge3, Point3, double>;
template class PrimitiveDistanceTemplate<Edge3, Edge3, double>;
template class PrimitiveDistanceTemplate<Point3, Point3, double>;
template class PrimitiveDistanceTemplate<Face, Point3, double>;

template class PrimitiveDistanceTemplate<Edge3, Point3, ADGrad<9>>;
template class PrimitiveDistanceTemplate<Edge3, Edge3, ADGrad<12>>;
template class PrimitiveDistanceTemplate<Point3, Point3, ADGrad<6>>;
template class PrimitiveDistanceTemplate<Face, Point3, ADGrad<12>>;

template class PrimitiveDistanceTemplate<Edge3, Point3, ADHessian<9>>;
template class PrimitiveDistanceTemplate<Edge3, Edge3, ADHessian<12>>;
template class PrimitiveDistanceTemplate<Point3, Point3, ADHessian<6>>;
template class PrimitiveDistanceTemplate<Face, Point3, ADHessian<12>>;
} // namespace ipc
