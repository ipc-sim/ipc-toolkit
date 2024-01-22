#include "vertex_vertex_3d.hpp"
#include "smooth_point_point.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <iterator>
#include <ipc/utils/logger.hpp>

namespace ipc {
    SmoothVertexVertex3Collision::SmoothVertexVertex3Collision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const std::array<double, 2> &dhats_,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhats_, mesh)
    {
        vertices[0] = primitive0;
        vertices[1] = primitive1;
        int id = 0;
        int v_id = 2;
        for (auto v : {primitive0, primitive1})
        {
            auto neighbors = mesh.find_vertex_adjacent_vertices(v);
            n_neighbors[id++] = neighbors.size();

            if (vertices.size() < v_id + neighbors.size())
                throw std::runtime_error("Max vertex size too small!");
            for (const auto &lv : neighbors)
                vertices[v_id++] = lv;
        }
        
        Eigen::VectorXd positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothVertexVertex3Collision::compute_types(
        const Eigen::VectorXd& positions, 
        ParameterType params)
    {
        auto points = slice_positions_large<double, 3>(positions);
        const Eigen::Ref<const RowVector3<double>>& va = points.row(0);
        const Eigen::Ref<const RowVector3<double>>& vb = points.row(1);
        const Eigen::Matrix<double, -1, 3> &ra = points.middleRows(2, n_neighbors[0]); 
        const Eigen::Matrix<double, -1, 3> &rb = points.bottomRows(n_neighbors[1]);

        RowVector3<double> direc = va - vb;
        const double dist = direc.norm();
        direc = direc / dist;
        // RowVector3<double> t, t_prev;

        if (dist*dist > get_eps())
            return false;
        
        params.eps = get_eps();
        bool return_val = smooth_point3_term_type(va, direc, ra, params.alpha, params.beta) &&
               smooth_point3_term_type(vb, -direc, rb, params.alpha, params.beta);

        if (dist < 1e-10)
            logger().warn("[vert-vert] dist {}, active {}, type 1 {}, type 2 {}", dist, return_val, smooth_point3_term_type(va, direc, ra, params.alpha, params.beta), smooth_point3_term_type(vb, -direc, rb, params.alpha, params.beta));

        if (return_val || (abs(evaluate_quadrature<double>(positions, params)) > 1e-12 ))
        {
            if (!return_val)
            {
                logger().error("[vert-vert] dist {}, active {}, type 1 {}, type 2 {}", dist, return_val, smooth_point3_term_type(va, direc, ra, params.alpha, params.beta), smooth_point3_term_type(vb, -direc, rb, params.alpha, params.beta));
                return true;
            }
        }

        return return_val;

        // return true;
        // assert(ra.rows() > 2);
        // assert(rb.rows() > 2);

        // bool normal_term = false;
        // bool tangent_term1 = true, tangent_term2 = true;
        // t_prev = ra.row(ra.rows()-1) - va;
        // for (int a = 0; a < ra.rows(); a++)
        // {
        //     t = ra.row(a) - va;
        //     tangent_term1 = tangent_term1 && direc.dot(t) / t.norm() / params.alpha > -1;
        //     normal_term = normal_term || -direc.dot(t_prev.cross(t).normalized()) / params.alpha > -1;
        //     std::swap(t, t_prev);
        // }

        // if (!normal_term)
        //     return false;

        // normal_term = false;
        // direc = -direc;
        // t_prev = rb.row(rb.rows()-1) - vb;
        // for (int b = 0; b < rb.rows(); b++)
        // {
        //     t = rb.row(b) - vb;
        //     tangent_term2 = tangent_term2 && direc.dot(t) / t.norm() / params.alpha > -1;
        //     normal_term = normal_term || -direc.dot(t_prev.cross(t).normalized()) / params.alpha > -1;
        //     std::swap(t, t_prev);
        // }

        // if (!normal_term || !tangent_term1 || !tangent_term2)
        //     return false;

        return true;
    }

    double SmoothVertexVertex3Collision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        auto points = slice_positions_large<double, 3>(positions);
        return (points.row(0) - points.row(1)).squaredNorm();
    }

    template <typename scalar> 
    scalar SmoothVertexVertex3Collision::evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const
    {
        auto points = slice_positions_large<scalar, 3>(positions);
        params.eps = get_eps();
        return smooth_point_point_potential_3d<scalar>(points.row(0), points.row(1), 
        points.middleRows(2, n_neighbors[0]), points.bottomRows(n_neighbors[1]), params);
    }

    double SmoothVertexVertex3Collision::operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
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
        //     logger().error("[vert-vert] {}: {} {}, {}, {}", (max_ - min_).maxCoeff(), g.transpose(), gc.transpose(), gl.transpose(), gr.transpose());
        // }

        return evaluate_quadrature<double>(positions, params);
    }

    Vector<double, -1, 3*max_vert_3d> SmoothVertexVertex3Collision::gradient(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarGrad<-1>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(ndofs());
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothVertexVertex3Collision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=AutodiffScalarHessian<-1>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(ndofs(), ndofs());
    }
}