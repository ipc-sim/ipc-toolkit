#include "edge_vertex_3d.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <iterator>
#include <ipc/utils/logger.hpp>

namespace ipc {
    SmoothEdgeVertex3Collision::SmoothEdgeVertex3Collision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhats_, mesh)
    {
        vertices[0] = primitive1;
        {
            auto tmp = mesh.find_edge_adjacent_vertices(primitive0);
            for (int i = 0; i < tmp.size(); i++)
                vertices[1+i] = tmp[i];
        }
        int v_id = 5;
        {
            auto neighbors = mesh.find_vertex_adjacent_vertices(primitive1);
            n_neighbors = neighbors.size();

            if (vertices.size() < v_id + neighbors.size())
                throw std::runtime_error("Max vertex size too small!");
            for (const auto &lv : neighbors)
                vertices[v_id++] = lv;
        }
        
        Eigen::VectorXd positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothEdgeVertex3Collision::compute_types(
        const Eigen::VectorXd& positions, 
        ParameterType params)
    {
        auto points = slice_positions_large<double, 3>(positions);
        dtype = point_edge_distance_type(points.row(0), points.row(1), points.row(2));

        // return true;
        // mollifier
        if (dtype != PointEdgeDistanceType::P_E)
            return false;
        
        params.eps = get_eps();
        return smooth_point_edge_potential_single_point_3d_type(points.row(0), points.bottomRows(n_neighbors), 
            points.row(1), points.row(2), points.row(3), points.row(4), params);
    }

    double SmoothEdgeVertex3Collision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        auto points = slice_positions_large<double, 3>(positions);
        return point_edge_distance(points.row(0), points.row(1), points.row(2));
    }

    template <typename scalar> 
    scalar SmoothEdgeVertex3Collision::evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const
    {
        auto points = slice_positions_large<scalar, 3>(positions);
        params.eps = get_eps();
        return smooth_point_edge_potential_single_point_3d<scalar>(points.row(0), points.bottomRows(n_neighbors), points.row(1), points.row(2), points.row(3), points.row(4), params);
    }

    double SmoothEdgeVertex3Collision::operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        assert(positions.size() == ndofs());

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

    Vector<double, -1, 3*max_vert_3d> SmoothEdgeVertex3Collision::gradient(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarGrad<-1>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(ndofs());
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothEdgeVertex3Collision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarHessian<-1>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(ndofs(), ndofs());
    }
}