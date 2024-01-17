#include "edge_edge_3d.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <iterator>
#include <ipc/utils/logger.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/edge_edge.hpp>

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
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhats_, mesh)
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
        
        Vector<double, 24> positions = dof(V, mesh.edges(), mesh.faces());
        has_edge_edge = compute_edge_edge_types(positions, param);
        has_point_edge = compute_point_edge_types(positions, param);
        is_active_ = has_edge_edge || has_point_edge;
    }

    bool SmoothEdgeEdge3Collision::compute_edge_edge_types(
        const Vector<double, 24>& positions, 
        const ParameterType &params)
    {
        std::array<Vector3<double>, 8> points_double = slice_positions<double, 8, 3>(positions);

        dtype = edge_edge_distance_type(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)],
            points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)]);
        
        if (dtype != EdgeEdgeDistanceType::EA_EB)
            return false;

        const Eigen::Vector3d direc = edge_edge_closest_point_direction<double>(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)], points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)], dtype);
        if (direc.norm() >= std::max(get_dhat(0), get_dhat(1)))
            return false;

        // normal term
        {
            bool all_skip = true;
            for (int e : {0, 1})
            {
                bool skip = true;
                for (int f : {2*e+0, 2*e+1})
                {
                    const Eigen::Vector3d normal = (points_double[face_to_vertex(f, 1)] - points_double[face_to_vertex(f, 0)]).cross(points_double[face_to_vertex(f, 2)] - points_double[face_to_vertex(f, 0)]).normalized();
                    const double val = (f >= 2 ? -1. : 1.) * direc.dot(normal) / direc.norm();
                    if (val < -params.alpha)
                        normal_types[f] = HEAVISIDE_TYPE::ZERO;
                    else
                    {
                        skip = false;
                        if (val > 0)
                            normal_types[f] = HEAVISIDE_TYPE::ONE;
                        else
                            normal_types[f] = HEAVISIDE_TYPE::VARIANT;
                    }
                }
                if (!skip)
                    all_skip = false;
            }
            if (all_skip)
                return false;
        }

        // tangent term
        {
            bool all_skip = true;
            for (int e : {0, 1})
            {
                bool skip = false;
                for (int f : {2*e+0, 2*e+1})
                {
                    const double val = (e > 0 ? 1. : -1.) * direc.dot(
                        point_line_closest_point_direction<double>(
                        points_double[face_to_vertex(f, 0)],
                        points_double[face_to_vertex(f, 1)],
                        points_double[face_to_vertex(f, 2)]).normalized()) / direc.norm();
                    if (val < -params.alpha)
                    {
                        tangent_types[f] = HEAVISIDE_TYPE::ZERO;
                        skip = true;
                    }
                    else if (val > 0)
                        tangent_types[f] = HEAVISIDE_TYPE::ONE;
                    else
                        tangent_types[f] = HEAVISIDE_TYPE::VARIANT;
                }
                if (!skip)
                    all_skip = false;
            }
            if (all_skip)
                return false;
        }

        int id = 0;
        for (int e : {0, 1})
            for (int v : {0, 1})
                edge_dtypes[id++] = point_edge_distance_type(points_double[face_to_vertex(2*e, v+1)], points_double[face_to_vertex(2*(1-e), 1)], points_double[face_to_vertex(2*(1-e), 2)]);
        
        // logger().debug("dtype {}, edge_types {} {} {} {}, tangent_types {} {} {} {}, normal_types {} {} {} {}",
        //     static_cast<int>(dtype), 
        //     static_cast<int>(edge_dtypes[0]),static_cast<int>(edge_dtypes[1]),static_cast<int>(edge_dtypes[2]),static_cast<int>(edge_dtypes[3]),
        //     static_cast<int>(tangent_types[0]),static_cast<int>(tangent_types[1]),static_cast<int>(tangent_types[2]),static_cast<int>(tangent_types[3]),
        //     static_cast<int>(normal_types[0]),static_cast<int>(normal_types[1]),static_cast<int>(normal_types[2]),static_cast<int>(normal_types[3]));
        return true;
    }

    double SmoothEdgeEdge3Collision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        std::array<Vector3<double>, 8> points = slice_positions<double, 8, 3>(positions);

        return edge_edge_distance(points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
            points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)], dtype);
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_quadrature(const Vector<double, 24>& positions, ParameterType params) const
    {
        scalar out(0.);
        if (has_edge_edge)
            out = out + evaluate_edge_edge_quadrature<scalar>(positions, params);
        if (has_point_edge)
            out = out + evaluate_point_edge_quadrature<scalar>(positions, params);
        return out;
    }

    bool SmoothEdgeEdge3Collision::compute_point_edge_types(
        const Vector<double, 24>& positions, 
        ParameterType params)
    {
        std::array<Vector3<double>, 8> points = slice_positions<double, 8, 3>(positions);
        bool skip = true;
        for (int le : {0, 2})
        {
            params.eps = pow(get_dhat(le/2), 2);
            for (int lv : {1, 2})
            {
                const auto &p = points[face_to_vertex(le, lv)];
                const auto &e0 = points[face_to_vertex(2-le, 1)];
                const auto &e1 = points[face_to_vertex(2-le, 2)];
                const auto &f0 = points[face_to_vertex(2-le+0, 0)];
                const auto &f1 = points[face_to_vertex(2-le+1, 0)];

                Vector3<double> direc = point_edge_closest_point_direction<double>(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
                const double dist_sqr = direc.squaredNorm();
                if (dist_sqr > params.eps)
                    continue;
                
                direc = direc / sqrt(dist_sqr);

                Vector3<double> t = point_line_closest_point_direction<double>(f0, e0, e1);
                if (-direc.dot(t) / t.norm() / params.alpha < -1)
                    continue;
                
                t = point_line_closest_point_direction<double>(f1, e1, e0);
                if (-direc.dot(t) / t.norm() / params.alpha < -1)
                    continue;

                const Vector3<double> normal1 = (e0 - f0).cross(e1 - f0);
                const Vector3<double> normal2 = (e1 - f1).cross(e0 - f1);
                if (direc.dot(normal1) / normal1.norm() / params.alpha < -1 && direc.dot(normal2) / normal2.norm() / params.alpha < -1)
                    continue;
                
                skip = false;
            }
        }
        return !skip;
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_point_edge_quadrature(const Vector<double, 24>& positions, ParameterType params) const
    {
        std::array<Vector3<scalar>, 8> points = slice_positions<scalar, 8, 3>(positions);
        scalar out(0.);
        for (int le : {0, 2})
        {
            params.eps = pow(get_dhat(le/2), 2);
            for (int lv : {1, 2})
            {
                out += smooth_point_edge_potential_single_point_3d<scalar>(
                    points[face_to_vertex(le, lv)],
                    points[face_to_vertex(2-le, 1)],points[face_to_vertex(2-le, 2)],
                    points[face_to_vertex(2-le+0, 0)],points[face_to_vertex(2-le+1, 0)], params);
            }
        }
        return out;
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_edge_edge_quadrature(const Vector<double, 24>& positions, ParameterType params) const
    {
        std::array<Vector3<scalar>, 8> points = slice_positions<scalar, 8, 3>(positions);
        const scalar out = smooth_edge_edge_potential_single_point<scalar>(
            points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
            points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)],
            points[face_to_vertex(0, 0)], points[face_to_vertex(1, 0)],
            points[face_to_vertex(2, 0)], points[face_to_vertex(3, 0)], 
            params, dhats, dtype, edge_dtypes, normal_types, tangent_types);
        
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

    double SmoothEdgeEdge3Collision::operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        assert(positions.size() == 24);
        return evaluate_quadrature<double>(positions, params);
    }

    Vector<double, -1, 3*max_vert_3d> SmoothEdgeEdge3Collision::gradient(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(24);
        using Diff=AutodiffScalarGrad<24>;
        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothEdgeEdge3Collision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(24);
        using Diff=AutodiffScalarHessian<24>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}