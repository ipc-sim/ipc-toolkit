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
    }

    SmoothEdgeEdge3Collision::SmoothEdgeEdge3Collision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhats_, mesh)
    {
        auto a = mesh.find_edge_adjacent_vertices(primitive0_);
        auto b = mesh.find_edge_adjacent_vertices(primitive1_);
        vertices = {{a[0], a[1], b[0], b[1], a[2], a[3], b[2], b[3]}};
        face_to_vertex << 4, 0, 1,
                          5, 1, 0,
                          6, 2, 3,
                          7, 3, 2;
        Vector<double, 24> positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothEdgeEdge3Collision::compute_types(
        const Vector<double, 24>& positions, 
        const ParameterType &params)
    {
        std::array<Vector3<double>, 8> points_double = slice_positions<double, 8, 3>(positions);

        dtype = edge_edge_distance_type(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)],
            points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)]);
        // return true;
        if (dtype != EdgeEdgeDistanceType::EA_EB)
            return false;

        const Eigen::Vector3d direc = edge_edge_closest_point_direction<double>(points_double[face_to_vertex(0, 1)], points_double[face_to_vertex(0, 2)], points_double[face_to_vertex(2, 1)], points_double[face_to_vertex(2, 2)], dtype);
        if (direc.squaredNorm() >= get_eps())
            return false;

        // if ((primitive0 == 1 && primitive1 == 8) || (primitive0 == 8 && primitive1 == 1))
        //     logger().error("dist {}, dtype {}", direc.squaredNorm(), static_cast<int>(dtype));

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
            bool all_skip = false;
            for (int e : {0, 1})
            {
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
                        all_skip = true;
                    }
                    else if (val > 0)
                        tangent_types[f] = HEAVISIDE_TYPE::ONE;
                    else
                        tangent_types[f] = HEAVISIDE_TYPE::VARIANT;
                }
            }
            if (all_skip)
                return false;
        }

        int id = 0;
        for (int e : {0, 1})
            for (int v : {0, 1})
                edge_dtypes[id++] = point_edge_distance_type(points_double[face_to_vertex(2*e, v+1)], points_double[face_to_vertex(2*(1-e), 1)], points_double[face_to_vertex(2*(1-e), 2)]);
        
        // logger().debug("before: edge {} {}, dtype {}, edge_types {} {} {} {}, tangent_types {} {} {} {}, normal_types {} {} {} {}",
        //     primitive0, primitive1,
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
        params.eps = get_eps();
        return evaluate_edge_edge_quadrature<scalar>(positions, params);
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_edge_edge_quadrature(const Vector<double, 24>& positions, ParameterType params) const
    {
        std::array<Vector3<scalar>, 8> points = slice_positions<scalar, 8, 3>(positions);
        params.eps = get_eps();
        const scalar out = smooth_edge_edge_potential_single_point<scalar>(
            points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
            points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)],
            points[face_to_vertex(0, 0)], points[face_to_vertex(1, 0)],
            points[face_to_vertex(2, 0)], points[face_to_vertex(3, 0)], 
            params, dtype, edge_dtypes, tangent_types, normal_types);

        return out;
    }

    double SmoothEdgeEdge3Collision::operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        assert(positions.size() == 24);

        // auto func = [&](const Eigen::VectorXd &x)
        // {
        //     return evaluate_quadrature<double>(positions, params);
        // };

        // Eigen::VectorXd g, gc, gl, gr;
        // my_finite_gradient(positions, func, gc, FD_RULE::CENTRAL);
        // my_finite_gradient(positions, func, gl, FD_RULE::LEFT);
        // my_finite_gradient(positions, func, gr, FD_RULE::RIGHT);
        // g = gradient(positions, params);
        
        // Eigen::VectorXd max_ = gr.array().max(gc.array().max(gl.array()));
        // Eigen::VectorXd min_ = gr.array().min(gc.array().min(gl.array()));
        // if ((max_ - min_).maxCoeff() > 1e-3 * max_.norm())
        // {
        //     logger().error("[edge-edge] {}: {} {}, {}, {}", (max_ - min_).maxCoeff(), g.transpose(), gc.transpose(), gl.transpose(), gr.transpose());
        // }
        
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