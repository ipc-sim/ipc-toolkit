#include "edge_edge_face.hpp"
#include "smooth_point_face.hpp"
#include "smooth_edge_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <iterator>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/finitediff.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/line_line.hpp>
#include <mutex>
std::mutex mut;

namespace ipc {
    namespace {
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

        std::array<Vector3<scalar>, 4> normals;
        for (int i = 0; i < 4; i++)
            normals[i] = (points[face_to_vertex(i, 1)] - points[face_to_vertex(i, 0)]).cross(points[face_to_vertex(i, 2)] - points[face_to_vertex(i, 0)]).normalized();
        
        scalar out = smooth_edge_edge_potential_single_point<scalar>(
            points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
            points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)],
            normals[0], normals[1], normals[2], normals[3], params, dtype);
        
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
        //     finite_gradient(positions, f, fgrad, FD_RULE::CENTRAL, 1e-8);

        //     if (out.getGradient().norm() > 1e-8)
        //     {
        //         double err = (out.getGradient() - fgrad).norm() / out.getGradient().norm();
        //         if (err > 1e-2)
        //         {
        //             mut.lock();
        //             finite_gradient(positions, f, fgrad1, FD_RULE::LEFT, 1e-8);
        //             finite_gradient(positions, f, fgrad2, FD_RULE::RIGHT, 1e-8);

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
        //             mut.unlock();
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