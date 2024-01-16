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
            std::vector<long> neighbors;
            auto iter = map.find(map.begin()->first);
            while (neighbors.empty() || iter->first != neighbors.front())
            {
                neighbors.push_back(iter->first);
                iter = map.find(iter->second);
                assert (iter != map.end());
            }
            assert(neighbors.size() == map.size());
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
        // auto points = slice_positions_large<double, 3>(positions);
        return true;
    }

    double SmoothVertexVertex3Collision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        auto points = slice_positions_large<double, 3>(positions);
        return std::numeric_limits<double>::max();
    }

    template <typename scalar> 
    scalar SmoothVertexVertex3Collision::evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const
    {
        auto points = slice_positions_large<scalar, 3>(positions);
        return smooth_point_point_potential_3d<scalar>(points.row(0), points.row(1), 
        points.middleRows(2, n_neighbors[0]), points.bottomRows(n_neighbors[1]), params, dhats);
    }

    std::array<long, max_vert_3d> SmoothVertexVertex3Collision::vertex_ids(
        const Eigen::MatrixXi& _edges, const Eigen::MatrixXi& _faces) const
    {
        return vertices;
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