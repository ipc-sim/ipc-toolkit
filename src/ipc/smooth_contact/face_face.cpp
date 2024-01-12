#include "face_face.hpp"
#include "smooth_point_face.hpp"
#include "smooth_edge_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <ipc/utils/logger.hpp>
// #include <mutex>
// std::mutex mut;

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

        // enum class FD_RULE { central, left, right };
        
        // void finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule)
        // {
        //     const double eps = 1e-7;
        //     grad.setZero(x.size());
        //     switch (rule)
        //     {
        //     case FD_RULE::central:
        //         for (int i = 0; i < x.size(); i++)
        //             for (int d : {-1, 1})
        //             {
        //                 auto y = x;
        //                 y(i) += d * eps;
        //                 grad(i) += d * f(y) / (2*eps);
        //             }
        //         break;
        //     case FD_RULE::left:
        //         for (int i = 0; i < x.size(); i++)
        //         {
        //                 auto y = x;
        //                 grad(i) += f(y) / eps;
        //                 y(i) -= eps;
        //                 grad(i) -= f(y) / eps;
        //         }
        //         break;
        //     case FD_RULE::right:
        //         for (int i = 0; i < x.size(); i++)
        //         {
        //                 auto y = x;
        //                 grad(i) -= f(y) / eps;
        //                 y(i) += eps;
        //                 grad(i) += f(y) / eps;
        //         }
        //         break;
        //     default:
        //     assert(false);
        //     }
        // }
    }

    std::array<long, 6> SmoothFaceFaceCollision::vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const
    {
        return {{faces(face0_id, 0), faces(face0_id, 1), faces(face0_id, 2),
                faces(face1_id, 0), faces(face1_id, 1), faces(face1_id, 2)}};
    }

    template <typename scalar> 
    scalar SmoothFaceFaceCollision::evaluate_quadrature(const Vector<double, -1, 18>& positions, const ParameterType &params) const
    {
        std::array<Vector3<scalar>, 6> points = slice_positions<scalar>(positions);
        std::array<Vector3<double>, 6> points_double = slice_positions<double>(positions);
        scalar out = scalar(0.);

        // face - vertex potential
        for (const int t : {0, 1})
        {
            const int tt = 1 - t;
            const scalar area = (points[t * 3 + 2] - points[t * 3 + 0]).cross(points[t * 3 + 1] - points[t * 3 + 0]).norm() / scalar(2.);
            
            const std::array<long, 3> ttv = {{vertices[tt * 3 + 0], vertices[tt * 3 + 1], vertices[tt * 3 + 2]}};
            
            for (const int i : {0, 1, 2})
            {
                const int p_id = t * 3 + i;
                if (std::find(ttv.begin(), ttv.end(), vertices[p_id]) != std::end(ttv))
                    continue;

                const PointTriangleDistanceType dtype = point_triangle_distance_type(
                    points_double[p_id], points_double[tt * 3 + 0], points_double[tt * 3 + 1], points_double[tt * 3 + 2]);
                out += (area / scalar(3.)) * smooth_point_face_potential_single_point<scalar>(
                    points[p_id], points[tt * 3 + 0], points[tt * 3 + 1], points[tt * 3 + 2], params, dtype);
                
                // for debugging derivatives using finite difference
                // if constexpr (std::is_same<scalar, AutodiffScalarGrad<18>>::value)
                // {
                //     mut.lock();
                //     Eigen::VectorXd fgrad, fgrad1, fgrad2;
                //     auto f = [&](const Eigen::VectorXd& x) {
                //         std::array<Vector3<double>, 6> points_ = slice_positions<double>(x);
                //         const PointTriangleDistanceType dtype_ = point_triangle_distance_type(
                //             x.segment<3>(3 * p_id), x.segment<3>(3 * (tt * 3 + 0)),
                //             x.segment<3>(3 * (tt * 3 + 1)), x.segment<3>(3 * (tt * 3 + 2)));
                //         return smooth_point_face_potential_single_point<double>(
                //         points_[p_id], points_[tt * 3 + 0], points_[tt * 3 + 1], points_[tt * 3 + 2], params, dtype_);
                //     };
                //     finite_gradient(positions, f, fgrad, FD_RULE::central);

                //     if (tmp.getGradient().norm() > 1e-8)
                //     {
                //         double err = (tmp.getGradient() - fgrad).norm() / tmp.getGradient().norm();
                //         if (err > 1e-3)
                //         {
                //             finite_gradient(positions, f, fgrad1, FD_RULE::left);
                //             finite_gradient(positions, f, fgrad2, FD_RULE::right);

                //             logger().error("p {}, t0 {}, t1 {}, t2 {}", p_id, tt * 3 + 0, tt * 3 + 1, tt * 3 + 2);
                //             logger().error("distance type {}", static_cast<int>(dtype));
                //             logger().error("err {}, norm {}", err, tmp.getGradient().norm());
                //             logger().error("positions {}", positions.transpose());
                //             logger().error("grad {}", tmp.getGradient().transpose());
                //             logger().error("fgrad {}", fgrad.transpose());
                //             logger().error("fgrad1 {}", fgrad1.transpose());
                //             logger().error("fgrad2 {}", fgrad2.transpose());
                //         }
                //     }
                //     mut.unlock();
                // }
            }
        }

        // edge - edge potential
        {
            const int t = 0;
            const int tt = 1;
            const scalar area = ((points[t * 3 + 2] - points[t * 3 + 0]).cross(points[t * 3 + 1] - points[t * 3 + 0]).norm() / scalar(2.)) * ((points[tt * 3 + 2] - points[tt * 3 + 0]).cross(points[tt * 3 + 1] - points[tt * 3 + 0]).norm() / scalar(2.));
            for (const int e0 : {0, 1, 2})
            {
                const std::array<long, 2> e0v = {{vertices[t * 3 + e0], vertices[t * 3 + (e0 + 1) % 3]}};
                for (const int e1 : {0, 1, 2})
                {
                    const std::array<long, 2> e1v = {{vertices[tt * 3 + e1], vertices[tt * 3 + (e1 + 1) % 3]}};
                    
                    // skip if two edges share at least one end point
                    if (std::find(e0v.begin(), e0v.end(), e1v[0]) != std::end(e0v) ||
                        std::find(e0v.begin(), e0v.end(), e1v[1]) != std::end(e0v))
                        continue;

                    const EdgeEdgeDistanceType dtype = edge_edge_distance_type(points_double[t * 3 + e0], points_double[t * 3 + (e0 + 1) % 3], points_double[tt * 3 + e1], points_double[tt * 3 + (e1 + 1) % 3]);

                    // Use original edge-edge for now
                    out += (area / scalar(9.)) * smooth_edge_edge_potential_pointwise<scalar>(points[t * 3 + e0], points[t * 3 + (e0 + 1) % 3], points[tt * 3 + e1], points[tt * 3 + (e1 + 1) % 3], params, dtype);
                }
            }
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

        // Eigen::VectorXd fgrad;
        // {
        //     auto f = [&](const Eigen::VectorXd& x) {
        //         return evaluate_quadrature<double>(x, params);
        //     };
        //     finite_gradient(positions, f, fgrad);
        // }

        Vector<double, -1, 18> grad = evaluate_quadrature<Diff>(positions, params).getGradient().head(18);

        // if (grad.norm() > 1e-8)
        // {
        //     double err = (grad - fgrad).norm() / grad.norm();
        //     if (err > 1e-3)
        //     {
        //         // std::cout << "err " << err << " grad norm " << grad.norm() << "\n";
        //         logger().error("err {}, norm {}", err, grad.norm());
        //         logger().error("positions {}", positions.transpose());
        //         logger().error("grad {}", grad.transpose());
        //         logger().error("fgrad {}", fgrad.transpose());
        //         exit(0);
        //     }
        //     else
        //         logger().warn("err {}, norm {}", err, grad.norm());
        // }

        return grad;
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