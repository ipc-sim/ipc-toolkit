#include "edge_edge_face.hpp"
#include "smooth_edge_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <iterator>
#include <ipc/utils/logger.hpp>
// #include <ipc/utils/finitediff.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/line_line.hpp>
// #include <mutex>
// std::mutex mut_edge;

namespace ipc {
    namespace {
        enum class FD_RULE { CENTRAL, LEFT, RIGHT };
        
        void my_finite_gradient(const Eigen::VectorXd& x, const std::function<double(const Eigen::VectorXd&)> &f, Eigen::VectorXd &grad, FD_RULE rule = FD_RULE::CENTRAL, const double eps = 1e-7)
        {
            grad.setZero(x.size());
            switch (rule)
            {
            case FD_RULE::CENTRAL:
                for (int i = 0; i < x.size(); i++)
                    for (int d : {-1, 1})
                    {
                        auto y = x;
                        y(i) += d * eps;
                        grad(i) += d * f(y) / (2*eps);
                    }
                break;
            case FD_RULE::LEFT:
                for (int i = 0; i < x.size(); i++)
                {
                        auto y = x;
                        grad(i) += f(y) / eps;
                        y(i) -= eps;
                        grad(i) -= f(y) / eps;
                }
                break;
            case FD_RULE::RIGHT:
                for (int i = 0; i < x.size(); i++)
                {
                        auto y = x;
                        grad(i) -= f(y) / eps;
                        y(i) += eps;
                        grad(i) += f(y) / eps;
                }
                break;
            default:
            assert(false);
            }
        }
        template <typename Iter>
        size_t index_of(Iter first, Iter last, const typename std::iterator_traits<Iter>::value_type& x)
        {
            size_t i = 0;
            while (first != last && *first != x)
            ++first, ++i;
            return i;
        }

        bool not_in(const long &a, const std::array<long, 2> &b)
        {
            return a != b[0] && a != b[1];
        }
    }

    SmoothEdgeEdge3Collision::SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh)
    : SmoothCollision(primitive0_, primitive1_, mesh)
    {
        std::array<long, 4> faces = {{mesh.edges_to_faces()(primitive0, 0), mesh.edges_to_faces()(primitive0, 1),
                  mesh.edges_to_faces()(primitive1, 0), mesh.edges_to_faces()(primitive1, 1)}};
        std::array<std::array<long, 2>, 2> edges = {{ {{mesh.edges()(primitive0, 0), mesh.edges()(primitive0, 1)}}, {{mesh.edges()(primitive1, 0), mesh.edges()(primitive1, 1)}} }};
        vertices.fill(-1);

        face_to_vertex.setConstant(-1);

        for (int j : {0, 1, 2, 3})
        {
            const int le = j / 2;
            int i;
            for (i = 0; i < 3; i++)
            {
                const auto va = mesh.faces()(faces[j], i);

                if (not_in(va, edges[le]))
                {
                    face_to_vertex(j, 0) = 4 + j;
                    vertices[4 + j] = va;
                    break;
                }
            }
            const auto vb = mesh.faces()(faces[j], (i+1) % 3);
            const auto vc = mesh.faces()(faces[j], (i+2) % 3);

            if (j%2 == 0)
            {
                vertices[j + 0] = vb;
                vertices[j + 1] = vc;
            }

            face_to_vertex(j, 1) = index_of(vertices.begin(), vertices.begin() + 4, vb);
            face_to_vertex(j, 2) = index_of(vertices.begin(), vertices.begin() + 4, vc);
            assert(face_to_vertex(j, 1) < 4 && face_to_vertex(j, 2) < 4);
        }
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_quadrature(const Vector<double, 24>& positions, const ParameterType &params) const
    {
        std::array<Vector3<scalar>, 8> points = slice_positions<scalar, 8, 3>(positions);
        std::array<Vector3<double>, 8> points_double = slice_positions<double, 8, 3>(positions);

        const EdgeEdgeDistanceType dtype = edge_edge_distance_type(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)],
            points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)]);
        const Eigen::Vector3d direc = edge_edge_closest_point_direction<double>(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)], points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)], dtype).normalized();
        
        std::array<PointEdgeDistanceType, 4> edge_dtypes;
        {
            edge_dtypes.fill(PointEdgeDistanceType::AUTO);
            int id = 0;
            for (int e : {0, 1})
                for (int v : {0, 1})
                    edge_dtypes[id++] = point_edge_distance_type(points_double[face_to_vertex(2*e, v+1)], points_double[face_to_vertex(2*(1-e), 1)], points_double[face_to_vertex(2*(1-e), 2)]);
        }
        
        std::array<HEAVISIDE_TYPE, 4> normal_types;
        {
            for (int f : {0, 1, 2, 3})
            {
                const Eigen::Vector3d normal = (points_double[face_to_vertex(f, 1)] - points_double[face_to_vertex(f, 0)]).cross(points_double[face_to_vertex(f, 2)] - points_double[face_to_vertex(f, 0)]).normalized();
                const double val = (f >= 2 ? -1. : 1.) * direc.dot(normal);
                if (val < -params.alpha)
                    return scalar(0.);
                    // normal_types[f] = HEAVISIDE_TYPE::ZERO;
                else if (val > 0)
                    normal_types[f] = HEAVISIDE_TYPE::ONE;
                else
                    normal_types[f] = HEAVISIDE_TYPE::VARIANT;
            }
        }

        std::array<HEAVISIDE_TYPE, 4> tangent_types;
        {
            for (int f : {0, 1, 2, 3})
            {
                const double val = (f >= 2 ? 1. : -1.) * direc.dot(
                    point_line_closest_point_direction<double>(
                    points_double[face_to_vertex(f, 0)],
                    points_double[face_to_vertex(f, 1)],
                    points_double[face_to_vertex(f, 2)]).normalized());
                if (val < -params.alpha)
                    tangent_types[f] = HEAVISIDE_TYPE::ZERO;
                else if (val > 0)
                    tangent_types[f] = HEAVISIDE_TYPE::ONE;
                else
                    tangent_types[f] = HEAVISIDE_TYPE::VARIANT;
            }
        }
        
        const scalar out = smooth_edge_edge_potential_single_point<scalar>(
            points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
            points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)],
            points[face_to_vertex(0, 0)], points[face_to_vertex(1, 0)],
            points[face_to_vertex(2, 0)], points[face_to_vertex(3, 0)], 
            params, dtype, edge_dtypes, normal_types, tangent_types);
        
        // if constexpr (std::is_same<scalar, AutodiffScalarGrad<24>>::value)
        // {
        //     std::unordered_set<long> vert_set;
        //     for (auto v : vertices)
        //         vert_set.insert(v);
        //     if (vert_set.size() < vertices.size())
        //         return out;

        //     Eigen::VectorXd fgrad, fgrad1, fgrad2;
        //     auto f = [&](const Eigen::VectorXd& x) {
        //         auto points_ = slice_positions<double, 8, 3>(x);
        //         const EdgeEdgeDistanceType dtype_ = edge_edge_distance_type(
        //             points_[face_to_vertex(0, 1)], points_[face_to_vertex(0, 2)],
        //             points_[face_to_vertex(2, 1)], points_[face_to_vertex(2, 2)]);

        //         std::array<Vector3<double>, 4> normals_;
        //         for (int i = 0; i < 4; i++)
        //             normals_[i] = (points_[face_to_vertex(i, 1)] - points_[face_to_vertex(i, 0)]).cross(points_[face_to_vertex(i, 2)] - points_[face_to_vertex(i, 0)]).normalized();
                
        //         return smooth_edge_edge_potential_single_point<double>(
        //             points_[face_to_vertex(0, 1)], points_[face_to_vertex(0, 2)],
        //             points_[face_to_vertex(2, 1)], points_[face_to_vertex(2, 2)],
        //             normals_[0], normals_[1], normals_[2], normals_[3], params, dtype_);
        //     };
        //     my_finite_gradient(positions, f, fgrad, FD_RULE::CENTRAL, 1e-8);

        //     if (out.getGradient().norm() > 1e-8)
        //     {
        //         double err = (out.getGradient() - fgrad).norm() / out.getGradient().norm();
        //         if (err > 1e-4)
        //         {
        //             mut_edge.lock();
        //             my_finite_gradient(positions, f, fgrad1, FD_RULE::LEFT, 1e-8);
        //             my_finite_gradient(positions, f, fgrad2, FD_RULE::RIGHT, 1e-8);

        //             logger().error("fa0 {} {} {}", face_to_vertex(0, 0), face_to_vertex(0, 1), face_to_vertex(0, 2));
        //             logger().error("fa1 {} {} {}", face_to_vertex(1, 0), face_to_vertex(1, 1), face_to_vertex(1, 2));
        //             logger().error("fb0 {} {} {}", face_to_vertex(2, 0), face_to_vertex(2, 1), face_to_vertex(2, 2));
        //             logger().error("fb1 {} {} {}", face_to_vertex(3, 0), face_to_vertex(3, 1), face_to_vertex(3, 2));

        //             double line_line_distance_ = line_line_distance(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)], points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)]);
        //             double vert_line_distance_ = std::numeric_limits<double>::max();
        //             {
        //                 for (int i : {2, 1})
        //                     vert_line_distance_ = std::min(vert_line_distance_, 
        //                     std::min(point_line_distance(points_double[face_to_vertex(0, i)], points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)]),
        //                              point_line_distance(points_double[face_to_vertex(2, i)], points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)])));
        //             }

        //             Vector3<double> u = points_double[face_to_vertex(0, 1)] - points_double[face_to_vertex(0, 2)];
        //             Vector3<double> v = points_double[face_to_vertex(2, 1)] - points_double[face_to_vertex(2, 2)];
        //             logger().error("distance type {}, parallel threshold {}", static_cast<int>(dtype), u.cross(v).squaredNorm() / u.squaredNorm() / v.squaredNorm());
        //             logger().error("line line distance {}", line_line_distance_);
        //             logger().error("vertex line distance {}", vert_line_distance_);
        //             logger().error("err {}, norm {}", err, out.getGradient().norm());
        //             logger().error("positions {}", positions.transpose());
        //             logger().error("grad {}", out.getGradient().transpose());
        //             logger().error("fgrad {}", fgrad.transpose());
        //             logger().error("fgrad1 {}", fgrad1.transpose());
        //             logger().error("fgrad2 {}", fgrad2.transpose());
        //             mut_edge.unlock();
        //         }
        //     }
        // }

        return out;
    }

    std::array<long, 8> SmoothEdgeEdge3Collision::vertex_ids(
        const Eigen::MatrixXi& _edges, const Eigen::MatrixXi& _faces) const
    {
        return vertices;
    }

    double SmoothEdgeEdge3Collision::operator()(const Vector<double, -1, 24>& positions, 
        const ParameterType &params) const
    {
        assert(positions.size() == 24);
        return evaluate_quadrature<double>(positions, params);
    }

    Vector<double, -1, 24> SmoothEdgeEdge3Collision::gradient(
        const Vector<double, -1, 24>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(24);
        using Diff=AutodiffScalarGrad<24>;
        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax<double, 24, 24> SmoothEdgeEdge3Collision::hessian(
        const Vector<double, -1, 24>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(24);
        using Diff=AutodiffScalarHessian<24>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}