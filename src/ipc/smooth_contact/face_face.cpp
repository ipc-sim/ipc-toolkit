#include "face_face.hpp"
#include "smooth_point_face.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    namespace {
        template <class T>
        std::array<Vector3<T>, 6> slice_positions(const Vector<double, 18> &positions)
        {
            std::array<Vector3<T>, 6> points;
            points.fill(Vector3<T>::Zero(3));
            
            for (int i = 0, id = 0; i < 6; i++)
                for (int d = 0; d < 3; d++, id++)
                    if constexpr (std::is_same<T, double>::value)
                        points[i](d) = positions(id);
                    else
                        points[i](d) = T(id, positions(id));

            return points;
        }
    }

    template <typename scalar> 
    scalar SmoothFaceFaceCollision::evaluate_quadrature(const Vector<double, -1, 18>& positions, const ParameterType &params) const
    {
        std::array<Vector3<scalar>, 6> points = slice_positions<scalar>(positions);
        scalar out = scalar(0.);

        for (const int t : {0, 1})
        {
            const int tt = 1 - t;
            const scalar area = (points[tt * 3 + 2] - points[tt * 3 + 0]).cross(points[tt * 3 + 1] - points[tt * 3 + 0]).norm() / scalar(2.);
            
            // face - vertex potential
            for (const int i : {0, 1, 2})
            {
                if (vertices[t * 3 + i] == vertices[tt * 3 + 0] ||
                    vertices[t * 3 + i] == vertices[tt * 3 + 1] ||
                    vertices[t * 3 + i] == vertices[tt * 3 + 2])
                    continue;

                const PointTriangleDistanceType dtype = point_triangle_distance_type(
                    positions.segment<3>(3 * (t * 3 + i)), positions.segment<3>(3 * (tt * 3 + 0)), 
                    positions.segment<3>(3 * (tt * 3 + 1)), positions.segment<3>(3 * (tt * 3 + 2)));
                out += (area / scalar(3.)) * smooth_point_face_potential_single_point<scalar>(
                    points[t * 3 + i], points[tt * 3 + 0], points[tt * 3 + 1], points[tt * 3 + 2], params, dtype);
            }

            // TODO: edge - edge potential
        }

        return out;
    }

    double SmoothFaceFaceCollision::compute_distance(const Vector<double, -1, 18>& positions) const
    {
        assert(false);
        return 0.;
    }

    Vector<double, -1, 18>
    SmoothFaceFaceCollision::compute_distance_gradient(const Vector<double, -1, 18>& positions) const
    {
        assert(false);
        return Vector<double, -1, 18>::Zero(18);
    }

    MatrixMax<double, 18, 18>
    SmoothFaceFaceCollision::compute_distance_hessian(const Vector<double, -1, 18>& positions) const
    {
        assert(false);
        return MatrixMax<double, 18, 18>::Zero(18, 18);
    }

    double SmoothFaceFaceCollision::operator()(const Vector<double, -1, 18>& positions, 
        const ParameterType &params) const
    {
        return evaluate_quadrature<double>(positions, params);
    }

    Vector<double, -1, 18> SmoothFaceFaceCollision::gradient(
        const Vector<double, -1, 18>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(18);
        using Diff=AutodiffScalarGrad<18>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(18);
    }

    MatrixMax<double, 18, 18> SmoothFaceFaceCollision::hessian(
        const Vector<double, -1, 18>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(18);
        using Diff=AutodiffScalarHessian<18>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(18, 18);
    }

}