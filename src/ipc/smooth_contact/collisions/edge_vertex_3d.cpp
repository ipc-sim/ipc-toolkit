#include "edge_vertex_3d.hpp"
#include <ipc/smooth_contact/pairs/smooth_point_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
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
                throw std::runtime_error("[edge-vert] Max vertex size too small! " + std::to_string(v_id + neighbors.size()));
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
        bool return_val = true;
        if (dtype != PointEdgeDistanceType::P_E)
            return_val = false;
        
        params.dhat = get_dhat();
        return_val = return_val && smooth_point_edge_potential_single_point_3d_type(points.row(0), points.bottomRows(n_neighbors), 
            points.row(1), points.row(2), points.row(3), points.row(4), params);

        // double dist = sqrt(point_edge_distance(points.row(0), points.row(1), points.row(2), dtype));
        // if (dist < 1e-10)
        //     logger().warn("[edge-vert] dist {}, active {}", dist, return_val);

        // if (return_val || (abs(evaluate_quadrature<double>(positions, params)) > 1e-15 ))
        // {
        //     if (!return_val)
        //     {
        //         Vector3<double> direc = point_edge_closest_point_direction<double>(points.row(0), points.row(1), points.row(2), PointEdgeDistanceType::AUTO).normalized();
        //         logger().error("[edge-vert] Wrong type! error {}, dist {}, edge_term {}, vert_term {}", 
        //             abs(evaluate_quadrature<double>(positions, params)), dist, smooth_edge3_term_type(direc, points.row(1), points.row(2), points.row(3), points.row(4), params.alpha, params.beta), smooth_point3_term_type(points.row(0), direc, points.bottomRows(n_neighbors), params.alpha, params.beta));
        //     }
        //     return true;
        // }

        return return_val;
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
        params.dhat = get_dhat();
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
        // if (std::max((g - max_).maxCoeff(), (max_ - min_).maxCoeff()) > 1e-3 * std::max(max_.norm(), min_.norm()))
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
        using Diff=AutodiffScalarGrad<-1, max_vert_3d*3>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(ndofs());
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothEdgeVertex3Collision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarHessian<-1, max_vert_3d*3>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(ndofs(), ndofs());
    }
}