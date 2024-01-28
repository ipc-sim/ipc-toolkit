#include "face_vertex.hpp"
#include <ipc/smooth_contact/pairs/smooth_point_face.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <ipc/utils/logger.hpp>

namespace ipc {
    SmoothFaceVertexCollision::SmoothFaceVertexCollision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhats_, mesh)
    {
        for (int lv : {0,1,2})
            vertices[lv+1] = mesh.faces()(primitive0, lv);

        vertices[0] = primitive1;

        {
            int v_id = 4;
            auto neighbors = mesh.find_vertex_adjacent_vertices(vertices[0]);
            n_neighbors = neighbors.size();

            if (vertices.size() < v_id + n_neighbors)
                throw std::runtime_error("[face-vert] Max vertex size too small!" + std::to_string(v_id + n_neighbors));
            for (const auto &lv : neighbors)
                vertices[v_id++] = lv;
        }
        
        Eigen::VectorXd positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothFaceVertexCollision::compute_types(
        const Eigen::VectorXd& positions, 
        ParameterType params)
    {
        auto points = slice_positions_large<double, 3>(positions);
        dtype = point_triangle_distance_type(points.row(0), points.row(1), points.row(2), points.row(3));
        const double dist = sqrt(point_triangle_distance(points.row(0), points.row(1), points.row(2), points.row(3), dtype));

        bool return_val = true;
        if (dist >= get_dhat())
            return_val = false;

        // return true;
        if (dtype != PointTriangleDistanceType::P_T)
            return_val = false;

        const Vector3<double> normal = (points.row(2) - points.row(1)).cross(points.row(3) - points.row(1));
        const double Phi = 1 - (points.row(0) - points.row(1)).dot(normal) / dist / normal.norm();
        if (Phi >= params.alpha)
            return_val = false;

        params.dhat = get_dhat();
        Vector3<double> direc = point_triangle_closest_point_direction<double>(points.row(0), points.row(1), points.row(2), points.row(3), dtype) / dist;
        return_val = return_val && smooth_point3_term_type(points.row(0), direc / direc.norm(), points.bottomRows(n_neighbors), params.alpha, params.beta);

        // if (dist < 1e-10)
        // {
        //     logger().warn("[face-vert] dist {}, active {}, face term {}, point term {}, error {}", dist, return_val, !(Phi >= params.alpha), smooth_point3_term_type(points.row(0), direc / direc.norm(), points.bottomRows(n_neighbors), params.alpha, params.beta), abs(evaluate_quadrature<double>(positions, params)));
        //     // return true;
        // }

        // if (return_val || (abs(evaluate_quadrature<double>(positions, params)) > 1e-15 ))
        // {
        //     if (!return_val)
        //     {
        //         logger().error("[face-vert] dist {}, active {}, face term {}, point term {}, error {}", dist, return_val, !(Phi >= params.alpha), smooth_point3_term_type(points.row(0), direc / direc.norm(), points.bottomRows(n_neighbors), params.alpha, params.beta), abs(evaluate_quadrature<double>(positions, params)));
        //         return true;
        //     }
        // }

        return return_val;
    }

    double SmoothFaceVertexCollision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        auto points = slice_positions_large<double, 3>(positions);
        return point_triangle_distance(points.row(0), points.row(1), points.row(2), points.row(3));
    }

    template <typename scalar> 
    scalar SmoothFaceVertexCollision::evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const
    {
        auto points = slice_positions_large<scalar, 3>(positions);
        params.dhat = get_dhat();
        return smooth_point_face_potential_single_point<scalar>(points.row(0), points.bottomRows(n_neighbors),
        points.row(1), points.row(2), points.row(3), params, dtype);
        return scalar(0.);
    }

    double SmoothFaceVertexCollision::operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        assert(positions.size() == ndofs());

        return evaluate_quadrature<double>(positions, params);
    }

    Vector<double, -1, 3*max_vert_3d> SmoothFaceVertexCollision::gradient(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarGrad<-1>;

        // auto func = [&](const Eigen::VectorXd &x)
        // {
        //     return evaluate_quadrature<double>(positions, params);
        // };

        // Eigen::VectorXd g, gc, gl, gr;
        // my_finite_gradient(positions, func, gc, FD_RULE::CENTRAL, 1e-8);
        // my_finite_gradient(positions, func, gl, FD_RULE::LEFT, 1e-8);
        // my_finite_gradient(positions, func, gr, FD_RULE::RIGHT, 1e-8);
        // g = evaluate_quadrature<Diff>(positions, params).getGradient().head(ndofs());
        
        // Eigen::VectorXd max_ = gr.array().max(gc.array().max(gl.array()));
        // Eigen::VectorXd min_ = gr.array().min(gc.array().min(gl.array()));
        // if (std::max((g - max_).maxCoeff(), (max_ - min_).maxCoeff()) > 1e-2 * std::max(max_.norm(), min_.norm()))
        // {
        //     logger().error("[face-vert] {}: {} {}, {}, {}", (max_ - min_).maxCoeff(), g.transpose(), gc.transpose(), gl.transpose(), gr.transpose());
        // }

        return evaluate_quadrature<Diff>(positions, params).getGradient().head(ndofs());
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothFaceVertexCollision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarHessian<-1>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(ndofs(), ndofs());
    }
}