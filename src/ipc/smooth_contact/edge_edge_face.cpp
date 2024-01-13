#include "edge_edge_face.hpp"
#include "smooth_point_face.hpp"
#include "smooth_edge_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <ipc/utils/logger.hpp>

namespace ipc {

    SmoothEdgeEdge3Collision::SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh)
    : SmoothCollision(primitive0_, primitive1_, mesh)
    {
        faces = {{mesh.edges_to_faces()(primitive0, 0), mesh.edges_to_faces()(primitive0, 1),
                  mesh.edges_to_faces()(primitive1, 0), mesh.edges_to_faces()(primitive1, 1)}};
        vertices = {{mesh.edges()(primitive0, 0), mesh.edges()(primitive0, 1),
                    mesh.edges()(primitive1, 0), mesh.edges()(primitive1, 1), 
                    -1, -1, -1, -1}};

        for (int j : {0, 1, 2, 3})
        {
            const auto f = faces[j];
            for (int i : {0, 1, 2})
            {
                if (mesh.faces()(f, i) != vertices[j/2 + 0] && mesh.faces()(f, i) != vertices[j/2 + 1])
                    vertices[j + 4] = mesh.faces()(f, i);
            }
            assert(vertices[j + 4] >= 0);
        }
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_quadrature(const Vector<double, 24>& positions, const ParameterType &params) const
    {
        return scalar(0.);
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