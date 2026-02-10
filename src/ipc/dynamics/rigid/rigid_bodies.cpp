#include "rigid_bodies.hpp"

#include <ipc/utils/logger.hpp>

#include <Eigen/Core>

#include <memory>

namespace ipc::rigid {

std::shared_ptr<RigidBodies> RigidBodies::build_from_meshes(
    const std::vector<Eigen::MatrixXd>& rest_positions,
    const std::vector<Eigen::MatrixXi>& edges,
    const std::vector<Eigen::MatrixXi>& faces,
    const std::vector<double>& densities,
    std::vector<Pose>& initial_poses)
{
    assert(rest_positions.size() == edges.size());
    assert(rest_positions.size() == faces.size());
    assert(rest_positions.size() == densities.size());

    size_t num_vertices = 0, num_edges = 0, num_faces = 0;
    std::vector<index_t> body_vertex_starts(rest_positions.size() + 1);
    std::vector<index_t> body_edge_starts(edges.size() + 1);
    std::vector<index_t> body_face_starts(faces.size() + 1);
    body_vertex_starts[0] = body_edge_starts[0] = body_face_starts[0] = 0;

    for (size_t i = 0; i < rest_positions.size(); ++i) {
        body_vertex_starts[i + 1] = (num_vertices += rest_positions[i].rows());
        body_edge_starts[i + 1] = (num_edges += edges[i].rows());
        body_face_starts[i + 1] = (num_faces += faces[i].rows());
    }

    Eigen::MatrixXd concat_rest_positions(
        num_vertices, rest_positions[0].cols());
    Eigen::MatrixXi concat_edges(num_edges, 2);
    Eigen::MatrixXi concat_faces(num_faces, 3);

    for (size_t i = 0; i < rest_positions.size(); ++i) {
        assert(rest_positions[i].size() > 0);
        concat_rest_positions.middleRows(
            body_vertex_starts[i], rest_positions[i].rows()) =
            rest_positions[i];
        if (edges[i].size() > 0) {
            concat_edges.middleRows(body_edge_starts[i], edges[i].rows()) =
                edges[i].array() + body_vertex_starts[i];
        }
        if (faces[i].size() > 0) {
            concat_faces.middleRows(body_face_starts[i], faces[i].rows()) =
                faces[i].array() + body_vertex_starts[i];
        }
    }

    return std::make_shared<RigidBodies>(
        concat_rest_positions, concat_edges, concat_faces, body_vertex_starts,
        body_edge_starts, body_face_starts, densities, initial_poses);
}

RigidBodies::RigidBodies(
    Eigen::ConstRef<Eigen::MatrixXd> _rest_positions,
    Eigen::ConstRef<Eigen::MatrixXi> _edges,
    Eigen::ConstRef<Eigen::MatrixXi> _faces,
    const std::vector<index_t>& _body_vertex_starts,
    const std::vector<index_t>& _body_edge_starts,
    const std::vector<index_t>& _body_face_starts,
    const std::vector<double>& densities,
    std::vector<Pose>& initial_poses)
    : CollisionMesh(_rest_positions, _edges, _faces)
    , body_vertex_starts(_body_vertex_starts)
    , body_edge_starts(_body_edge_starts)
    , body_face_starts(_body_face_starts)
{
    assert(body_vertex_starts.size() == body_edge_starts.size());
    assert(body_edge_starts.size() == body_face_starts.size());
    assert(body_vertex_starts.back() == num_vertices());
    assert(body_edge_starts.back() == num_edges());
    assert(body_face_starts.back() == num_faces());
    assert(initial_poses.size() == body_vertex_starts.size() - 1);

    bodies.reserve(body_vertex_starts.size() - 1);
    for (size_t i = 0; i < body_vertex_starts.size() - 1; ++i) {
        bodies.emplace_back(
            m_rest_positions.middleRows(
                body_vertex_starts[i],
                body_vertex_starts[i + 1] - body_vertex_starts[i]),
            edges().middleRows(
                       body_edge_starts[i],
                       body_edge_starts[i + 1] - body_edge_starts[i])
                    .array()
                - body_vertex_starts[i],
            faces().middleRows(
                       body_face_starts[i],
                       body_face_starts[i + 1] - body_face_starts[i])
                    .array()
                - body_vertex_starts[i],
            densities[i], initial_poses[i]);
        logger().info(
            "Initial pose: position={}, rotation={}",
            initial_poses[i].position.transpose(),
            initial_poses[i].rotation.transpose());
    }

    assert(body_vertex_starts.size() > 0);
    this->can_collide = [this](size_t vi, size_t vj) {
        assert(body_vertex_starts.size() > 0);
        return this->vertex_to_body(vi) != this->vertex_to_body(vj);
    };
    assert(!this->can_collide(0, 0));
}

Eigen::VectorXd RigidBodies::to_rigid_dof(
    const std::vector<Pose>& poses, Eigen::ConstRef<Eigen::VectorXd> g) const
{
    assert(poses.size() == num_bodies());
    assert(g.size() == ndof());

    const int ndof_per_body = poses[0].ndof();
    Eigen::VectorXd result(poses.size() * ndof_per_body);

    // Apply the chain rule to map the gradient from the collision mesh to the
    // rigid degrees of freedom.
    for (size_t i = 0; i < bodies.size(); ++i) {
        assert(poses[i].ndof() == ndof_per_body);
        const index_t start = body_vertex_starts[i];
        const index_t end = body_vertex_starts[i + 1];
        const index_t nV = end - start;

        result.segment(i * ndof_per_body, ndof_per_body) =
            poses[i]
                .transform_vertices_jacobian(
                    rest_positions().middleRows(start, end - start))
                .transpose()
            * g.segment(dim() * start, dim() * nV);
    }

    return result;
}

Eigen::MatrixXd RigidBodies::to_rigid_dof(
    const std::vector<Pose>& poses,
    Eigen::ConstRef<Eigen::VectorXd> g,
    const Eigen::SparseMatrix<double>& H) const
{
    assert(poses.size() == num_bodies());
    assert(g.size() == ndof());
    assert(H.rows() == ndof() && H.cols() == ndof());
    assert(H.isApprox(H.transpose())); // Hessian should be symmetric

    const int ndof_per_body = poses[0].ndof();
    Eigen::MatrixXd result(
        poses.size() * ndof_per_body, poses.size() * ndof_per_body);

    std::vector<Eigen::MatrixXd> dV_dx(bodies.size());
    std::vector<Eigen::MatrixXd> d2V_dx2(bodies.size());
    for (size_t i = 0; i < bodies.size(); ++i) {
        assert(poses[i].ndof() == ndof_per_body);

        // Get the rest positions of the vertices for the i-th body
        const auto VR = body_rest_positions(i);

        // Compute the Jacobian and Hessian of the transformed vertices with
        // respect to the pose for the i-th body.
        dV_dx[i] = poses[i].transform_vertices_jacobian(VR);
        assert(
            dV_dx[i].rows() == dim() * body_num_vertices(i)
            && dV_dx[i].cols() == ndof_per_body);

        d2V_dx2[i] = poses[i].transform_vertices_hessian(VR);
        assert(d2V_dx2[i].rows() == dim() * body_num_vertices(i));
        assert(d2V_dx2[i].cols() == ndof_per_body * ndof_per_body);
    }

    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i; j < bodies.size(); ++j) {
            auto result_ij = result.block(
                i * ndof_per_body, j * ndof_per_body, ndof_per_body,
                ndof_per_body);

            result_ij = dV_dx[i].transpose()
                * H.block(
                    dim() * body_vertex_starts[i],
                    dim() * body_vertex_starts[j], //
                    dim() * body_num_vertices(i),  //
                    dim() * body_num_vertices(j))
                * dV_dx[j];

            if (i == j) {
                const index_t start = body_vertex_starts[i];
                const index_t end = body_vertex_starts[i + 1];
                auto gi = g.segment(dim() * start, dim() * (end - start));

                // Add the contribution from the second derivative of the
                // transformation (only for the diagonal blocks)
                for (int jj = 0; jj < ndof_per_body; ++jj) {
                    for (int ii = 0; ii < ndof_per_body; ++ii) {
                        auto d2Vi_dxii_dxjj =
                            d2V_dx2[i].col(ii * ndof_per_body + jj);
                        result_ij(ii, jj) += gi.dot(d2Vi_dxii_dxjj);
                    }
                }
            } else {
                // Fill in the symmetric block
                result.block(
                    j * ndof_per_body, i * ndof_per_body, ndof_per_body,
                    ndof_per_body) = result_ij.transpose();
            }
        }
    }

    // return project_to_psd(hess, project_hessian_to_psd);
    return result;
}

} // namespace ipc::rigid