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
            std::unordered_map<long, long> map;
            for (auto f : mesh.vertices_to_faces()[v])
            {
                for (int lv = 0; lv < 3; lv++)
                {
                    if (mesh.faces()(f, lv) == v)
                    {
                        map[mesh.faces()(f, (lv+1)%3)] = mesh.faces()(f, (lv+2)%3);
                        break;
                    }
                }
            }
            if (mesh.vertices_to_faces()[v].size() != map.size())
                throw std::runtime_error("Non-manifold vertex! Map size smaller than neighbor!");
            
            std::vector<long> neighbors;
            auto iter = map.find(map.begin()->first);
            while (neighbors.empty() || iter->first != neighbors.front())
            {
                neighbors.push_back(iter->first);
                iter = map.find(iter->second);
                if (iter == map.end())
                {
                    logger().error("neighbor faces {}, map {}", mesh.vertices_to_faces()[v].size(), map);
                    throw std::runtime_error("Non-manifold vertex! Cannot find next neighbor!");
                }
            }
            if (neighbors.size() != map.size())
                throw std::runtime_error("Non-manifold vertex!");
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
        const ParameterType &params)
    {
        auto points = slice_positions_large<double, 3>(positions);
        const Eigen::Ref<const RowVector3<double>>& va = points.row(0);
        const Eigen::Ref<const RowVector3<double>>& vb = points.row(1);
        const Eigen::Matrix<double, -1, 3> &ra = points.middleRows(2, n_neighbors[0]); 
        const Eigen::Matrix<double, -1, 3> &rb = points.bottomRows(n_neighbors[1]);

        RowVector3<double> direc = va - vb;
        const double dist = direc.norm();
        direc = direc / dist;
        RowVector3<double> t, t_prev;

        if (dist > std::max(dhats[0], dhats[1]))
            return false;

        assert(ra.rows() > 2);
        assert(rb.rows() > 2);

        bool normal_term = false;
        bool tangent_term1 = true, tangent_term2 = true;
        t_prev = ra.row(ra.rows()-1) - va;
        for (int a = 0; a < ra.rows(); a++)
        {
            t = ra.row(a) - va;
            tangent_term1 = tangent_term1 && direc.dot(t) / t.norm() / params.alpha > -1;
            normal_term = normal_term || -direc.dot(t_prev.cross(t).normalized()) / params.alpha > -1;
            std::swap(t, t_prev);
        }

        if (!normal_term)
            return false;

        normal_term = false;
        direc = -direc;
        t_prev = rb.row(rb.rows()-1) - vb;
        for (int b = 0; b < rb.rows(); b++)
        {
            t = rb.row(b) - vb;
            tangent_term2 = tangent_term2 && direc.dot(t) / t.norm() / params.alpha > -1;
            normal_term = normal_term || -direc.dot(t_prev.cross(t).normalized()) / params.alpha > -1;
            std::swap(t, t_prev);
        }

        if (!normal_term)
            return false;
        
        if (!tangent_term1 && !tangent_term2)
            return false;

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
        return smooth_point_point_potential_3d<scalar>(points.row(0), points.row(1), 
        points.middleRows(2, n_neighbors[0]), points.bottomRows(n_neighbors[1]), params, dhats);
    }

    double SmoothVertexVertex3Collision::operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        assert(positions.size() == ndofs());
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