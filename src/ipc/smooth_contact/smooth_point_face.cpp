#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include "smooth_point_face.hpp"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <stdexcept>

namespace ipc {
    template <typename scalar>
    scalar smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params)
    {
        Vector3<scalar> tangent1 = v1 - v0;
        Vector3<scalar> tangent2 = v2 - v0;
        Vector3<scalar> sample = v0 + tangent1 * scalar(uv(0)) + tangent2 * scalar(uv(1));
        Vector3<scalar> normal = tangent1.cross(tangent2).normalized();

        const scalar dist_sqr = (p - sample).squaredNorm();
        const scalar Phi = normal.cross(p - sample).squaredNorm() / dist_sqr;

        return inv_barrier(dist_sqr, params.eps, params.r) * cubic_spline(Phi * (2. / params.alpha));
    }

    template double smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Ref<const Vector3<double>>& v0,
        const Eigen::Ref<const Vector3<double>>& v1,
        const Eigen::Ref<const Vector3<double>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params);
    template AutodiffScalarGrad<12> smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params);
    template AutodiffScalarHessian<12> smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params);

    template <typename scalar>
    scalar smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const ParameterType &params)
    {
        Eigen::Matrix<double, -1, 2> pts;
        Eigen::VectorXd weights;
        ipc::triangle_quadrature(params.n_quadrature, pts, weights);

        scalar val(0.);
        for (int i = 0; i < pts.size(); i++)
        {
            val += smooth_point_face_potential_pointwise<scalar>(p, v0, v1, v2, pts.row(i), params) * weights(i);
        }

        return val;
    }

    template double smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Ref<const Vector3<double>>& v0,
        const Eigen::Ref<const Vector3<double>>& v1,
        const Eigen::Ref<const Vector3<double>>& v2,
        const ParameterType &params);
    template AutodiffScalarGrad<12> smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v2,
        const ParameterType &params);
    template AutodiffScalarHessian<12> smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v2,
        const ParameterType &params);

    template <typename scalar>
    scalar smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype)
    {
        const Vector3<scalar> tangent1 = v1 - v0;
        const Vector3<scalar> tangent2 = v2 - v0;
        const Vector3<scalar> tangent3 = v2 - v1;
        const Vector3<scalar> normal = tangent1.cross(tangent2);
        const scalar area_sqr = normal.squaredNorm();
        const Vector3<scalar> pos = p - v0;

        Vector3<scalar> diff;
        switch (dtype)
        {
        case PointTriangleDistanceType::P_E0:
            diff = pos - (pos.dot(tangent1) / tangent1.squaredNorm()) * tangent1;
            break;
        case PointTriangleDistanceType::P_E1:
            diff = p - (v1 + (pos.dot(tangent3) / tangent3.squaredNorm()) * tangent3);
            break;
        case PointTriangleDistanceType::P_E2:
            diff = pos - (pos.dot(tangent2) / tangent2.squaredNorm()) * tangent2;
            break;
        case PointTriangleDistanceType::P_T0:
            diff = p - v0;
            break;
        case PointTriangleDistanceType::P_T1:
            diff = p - v1;
            break;
        case PointTriangleDistanceType::P_T2:
            diff = p - v2;
            break;
        case PointTriangleDistanceType::P_T:
            diff = (normal.dot(pos) / area_sqr) * normal;
            break;
        default:
            assert(false);
        }

        const scalar dist_sqr = diff.squaredNorm();
        const scalar Phi = cross2_sqr<scalar>(diff, normal) / dist_sqr / area_sqr;
        return sqrt(area_sqr) * inv_barrier(dist_sqr, params.eps, params.r) * cubic_spline(Phi * (2. / params.alpha));
    }

    template double smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Ref<const Vector3<double>>& v0,
        const Eigen::Ref<const Vector3<double>>& v1,
        const Eigen::Ref<const Vector3<double>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype);
    template AutodiffScalarGrad<12> smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<12>>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype);
    template AutodiffScalarHessian<12> smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<12>>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype);
}