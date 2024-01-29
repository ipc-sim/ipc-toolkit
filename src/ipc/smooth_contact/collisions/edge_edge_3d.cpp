#include "edge_edge_3d.hpp"
#include <ipc/smooth_contact/pairs/smooth_edge_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/distance/edge_edge.hpp>

namespace ipc {

    SmoothEdgeEdge3Collision::SmoothEdgeEdge3Collision(
    long primitive0_,
    long primitive1_,
    const CollisionMesh &mesh,
    const ParameterType &param,
    const double &dhat,
    const Eigen::MatrixXd &V): SmoothCollision<max_vert_3d>(primitive0_, primitive1_, dhat, mesh)
    {
        auto a = mesh.find_edge_adjacent_vertices(primitive0_);
        auto b = mesh.find_edge_adjacent_vertices(primitive1_);
        // vertices = {{a[0], a[1], b[0], b[1], a[2], a[3], b[2], b[3]}};
        vertices.fill(-1);
        vertices[0] = a[0];
        vertices[1] = a[1];
        vertices[4] = a[2];
        vertices[5] = a[3];
        vertices[2] = b[0];
        vertices[3] = b[1];
        vertices[6] = b[2];
        vertices[7] = b[3];
        face_to_vertex << 4, 0, 1,
                          5, 1, 0,
                          6, 2, 3,
                          7, 3, 2;
        Vector<double, 24> positions = dof(V, mesh.edges(), mesh.faces());
        is_active_ = compute_types(positions, param);
    }

    bool SmoothEdgeEdge3Collision::compute_types(
        const Vector<double, 24>& positions, 
        ParameterType params)
    {
        auto points = slice_positions<double, 8, 3>(positions);

        // already computed in SmoothCollisionsBuilder before class construction
        dtype = EdgeEdgeDistanceType::EA_EB;

        params.dhat = get_dhat();
        return smooth_edge_edge_potential_type(
            points.row(face_to_vertex(0, 1)), points.row(face_to_vertex(0, 2)),
            points.row(face_to_vertex(2, 1)), points.row(face_to_vertex(2, 2)),
            points.row(face_to_vertex(0, 0)), points.row(face_to_vertex(1, 0)),
            points.row(face_to_vertex(2, 0)), points.row(face_to_vertex(3, 0)), 
            params, dtype, otypes, mtypes);

        // double dist = sqrt(edge_edge_distance(points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
        //     points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)], dtype));

        // if (dist < 1e-10)
        // {
        //     logger().warn("[edge-edge] Dtype {}, Dist {}, active {}, residual {}", static_cast<int>(dtype), dist, return_val, abs(evaluate_quadrature<double>(positions, params, true)));
        //     // return true;
        // }

        // if (return_val || (abs(evaluate_quadrature<double>(positions, params)) > 1e-15 ))
        // {
        //     if (!return_val)
        //     {
        //         logger().error("[edge-edge] Dtype {}, Dist {}, active {}, residual {}", static_cast<int>(dtype), dist, return_val, abs(evaluate_quadrature<double>(positions, params, true)));
        //         return true;
        //     }
        // }

        // logger().debug("before: edge {} {}, dtype {}, edge_types {} {} {} {}, tangent_types {} {} {} {}, normal_types {} {} {} {}",
        //     primitive0, primitive1,
        //     static_cast<int>(dtype), 
        //     static_cast<int>(edge_dtypes[0]),static_cast<int>(edge_dtypes[1]),static_cast<int>(edge_dtypes[2]),static_cast<int>(edge_dtypes[3]),
        //     static_cast<int>(tangent_types[0]),static_cast<int>(tangent_types[1]),static_cast<int>(tangent_types[2]),static_cast<int>(tangent_types[3]),
        //     static_cast<int>(normal_types[0]),static_cast<int>(normal_types[1]),static_cast<int>(normal_types[2]),static_cast<int>(normal_types[3]));
    }

    double SmoothEdgeEdge3Collision::compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const
    {
        auto points = slice_positions<double, 8, 3>(positions);

        return edge_edge_distance(points.row(face_to_vertex(0, 1)), points.row(face_to_vertex(0, 2)),
            points.row(face_to_vertex(2, 1)), points.row(face_to_vertex(2, 2)), dtype);
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_quadrature(const Vector<double, 24>& positions, ParameterType params) const
    {
        auto points = slice_positions<scalar, 8, 3>(positions);
        params.dhat = get_dhat();
        return smooth_edge_edge_potential_single_point<scalar>(
            points.row(face_to_vertex(0, 1)), points.row(face_to_vertex(0, 2)),
            points.row(face_to_vertex(2, 1)), points.row(face_to_vertex(2, 2)),
            points.row(face_to_vertex(0, 0)), points.row(face_to_vertex(1, 0)),
            points.row(face_to_vertex(2, 0)), points.row(face_to_vertex(3, 0)), 
            params, dtype, otypes, mtypes);
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
        using Diff=ADGrad<24>;

        // auto func = [&](const Eigen::VectorXd &x)
        // {
        //     return evaluate_quadrature<double>(positions, params);
        // };

        // Eigen::VectorXd g, gc, gl, gr;
        // my_finite_gradient(positions, func, gc, FD_RULE::CENTRAL, 1e-8);
        // my_finite_gradient(positions, func, gl, FD_RULE::LEFT, 1e-8);
        // my_finite_gradient(positions, func, gr, FD_RULE::RIGHT, 1e-8);
        // g = evaluate_quadrature<Diff>(positions, params).getGradient();
        
        // Eigen::VectorXd max_ = gr.array().max(gc.array().max(gl.array()));
        // Eigen::VectorXd min_ = gr.array().min(gc.array().min(gl.array()));
        // if (std::max((g - max_).maxCoeff(), (max_ - min_).maxCoeff()) > 1e-2 * std::max(std::max(max_.norm(), min_.norm()), 1e-2))
        // {
        //     logger().error("[edge-edge] Not C1! fd err {:.5e}, grad err {:.5e}, norm {:.5e}", (max_ - min_).maxCoeff(), std::max((g - max_).maxCoeff(), (max_ - min_).maxCoeff()), std::max(max_.norm(), min_.norm()));
        // }

        return evaluate_quadrature<Diff>(positions, params).getGradient();
    }

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> SmoothEdgeEdge3Collision::hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const
    {
        DiffScalarBase::setVariableCount(24);
        using Diff=ADHessian<24>;
        return evaluate_quadrature<Diff>(positions, params).getHessian();
    }
}