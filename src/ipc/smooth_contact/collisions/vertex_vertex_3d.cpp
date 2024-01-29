#include "vertex_vertex_3d.hpp"
#include <ipc/smooth_contact/pairs/smooth_point_point.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <ipc/utils/logger.hpp>

namespace ipc {
    SmoothVertexVertex3Collision::SmoothVertexVertex3Collision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const double &dhat,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhat, mesh)
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
                throw std::runtime_error("[vert-vert] Max vertex size too small!" + std::to_string(v_id + neighbors.size()));
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
        auto points = slice_positions<double, -1, 3>(positions);
        const Eigen::Ref<const RowVector3<double>>& va = points.row(0);
        const Eigen::Ref<const RowVector3<double>>& vb = points.row(1);
        const Eigen::Matrix<double, -1, 3> &ra = points.middleRows(2, n_neighbors[0]); 
        const Eigen::Matrix<double, -1, 3> &rb = points.bottomRows(n_neighbors[1]);

        RowVector3<double> direc = va - vb;
        const double dist = direc.norm();
        direc = direc / dist;
        // RowVector3<double> t, t_prev;

        if (dist > get_dhat())
            return false;
        
        params.dhat = get_dhat();
        // otypes are computed here, if the return value is false, otypes are not initialized. It's not a bug because the potential is not computed in that case.
        bool return_val = smooth_point3_term_type(va, direc, ra, params.alpha, params.beta, otypes[0]) &&
               smooth_point3_term_type(vb, -direc, rb, params.alpha, params.beta, otypes[1]);

        // if (dist < 1e-10)
        //     logger().warn("[vert-vert] dist {}, active {}, type 1 {}, type 2 {}", dist, return_val, smooth_point3_term_type(va, direc, ra, params.alpha, params.beta), smooth_point3_term_type(vb, -direc, rb, params.alpha, params.beta));

        // if (return_val || (abs(evaluate_quadrature<double>(positions, params)) > 1e-15 ))
        // {
        //     if (!return_val)
        //     {
        //         logger().error("[vert-vert] dist {}, active {}, type 1 {}, type 2 {}", dist, return_val, smooth_point3_term_type(va, direc, ra, params.alpha, params.beta), smooth_point3_term_type(vb, -direc, rb, params.alpha, params.beta));
        //         return true;
        //     }
        // }

        return return_val;
    }

    double SmoothVertexVertex3Collision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        auto points = slice_positions<double, -1, 3>(positions);
        return (points.row(0) - points.row(1)).squaredNorm();
    }

    template <typename scalar> 
    scalar SmoothVertexVertex3Collision::evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const
    {
        auto points = slice_positions<scalar, -1, 3>(positions);
        params.dhat = get_dhat();
        return smooth_point_point_potential_3d<scalar>(points.row(0), points.row(1), 
        points.middleRows(2, n_neighbors[0]), points.bottomRows(n_neighbors[1]), params, otypes);
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
        using Diff=ADGrad<-1, max_vert_3d * 3>;
        return evaluate_quadrature<Diff>(positions, params).getGradient().head(ndofs());
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothVertexVertex3Collision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(ndofs());
        using Diff=ADHessian<-1, max_vert_3d * 3>;
        return evaluate_quadrature<Diff>(positions, params).getHessian().topLeftCorner(ndofs(), ndofs());
    }
}