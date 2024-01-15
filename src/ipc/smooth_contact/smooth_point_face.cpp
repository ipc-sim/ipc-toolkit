#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/distance_autodiff.hpp>
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

        if (Phi > params.alpha)
            return scalar(0.);
        if (dist_sqr > params.eps)
            return scalar(0.);

        return inv_barrier(dist_sqr / params.eps, params.r) * cubic_spline(Phi * (2. / params.alpha));
    }

    template double smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Ref<const Vector3<double>>& v0,
        const Eigen::Ref<const Vector3<double>>& v1,
        const Eigen::Ref<const Vector3<double>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params);
    template AutodiffScalarGrad<18> smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& v2,
        const Vector2<double> &uv,
        const ParameterType &params);
    template AutodiffScalarHessian<18> smooth_point_face_potential_pointwise(
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& v2,
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
    template AutodiffScalarGrad<18> smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<18>>>& v2,
        const ParameterType &params);
    template AutodiffScalarHessian<18> smooth_point_face_potential_quadrature(
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<18>>>& v2,
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
        const Vector3<scalar> normal = (v1 - v0).cross(v2 - v0).normalized();
        const scalar dist_sqr = point_triangle_distance(p, v0, v1, v2, dtype);
        const scalar Phi = 1 - (p - v0).dot(normal) / sqrt(dist_sqr); // cross2_sqr<scalar>(diff, normal) / dist_sqr / normal_len_sqr;

        if (Phi > params.alpha)
            return scalar(0.);

        return inv_barrier(dist_sqr / params.eps, params.r) * cubic_spline(Phi * (2. / params.alpha));
    }

    template double smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Ref<const Vector3<double>>& v0,
        const Eigen::Ref<const Vector3<double>>& v1,
        const Eigen::Ref<const Vector3<double>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype);
    template AutodiffScalarGrad<24> smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<24>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<24>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<24>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarGrad<24>>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype);
    template AutodiffScalarHessian<24> smooth_point_face_potential_single_point(
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<24>>>& p,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<24>>>& v0,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<24>>>& v1,
        const Eigen::Ref<const Vector3<AutodiffScalarHessian<24>>>& v2,
        const ParameterType &params,
        const PointTriangleDistanceType &dtype);
}