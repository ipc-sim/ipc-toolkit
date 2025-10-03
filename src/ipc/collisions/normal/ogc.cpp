#include "ogc.hpp"

#include <ipc/tangent/closest_point.hpp>

namespace ipc::ogc {

bool check_vertex_feasible_region(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<VectorMax3d> x,
    const index_t vi)
{
    assert(mesh.are_adjacencies_initialized());
    assert(x.size() == vertices.cols());

    const VectorMax3d x_vi = vertices.row(vi);
    const VectorMax3d v_to_x = x - x_vi;

    return std::all_of(
        mesh.vertex_vertex_adjacencies()[vi].begin(),
        mesh.vertex_vertex_adjacencies()[vi].end(), [&](const index_t vj) {
            return v_to_x.dot(x_vi - vertices.row(vj)) >= 0;
        });
}

bool check_edge_feasible_region(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const index_t xi,
    const index_t ei)
{
    assert(mesh.are_adjacencies_initialized());
    assert(vertices.cols() == 3); // 3D only

    const Eigen::Vector3d x = vertices.row(xi);
    const Eigen::Vector3d x0 = vertices.row(mesh.edges()(ei, 0));
    const Eigen::Vector3d x1 = vertices.row(mesh.edges()(ei, 1));

    const double alpha = point_edge_closest_point(x, x0, x1);
    if (alpha <= 0 || alpha >= 1) {
        // Vertex xi is not in the feasible region of edge ei
        return false;
    }

    return std::all_of(
        mesh.edge_vertex_adjacencies()[ei].begin(),
        mesh.edge_vertex_adjacencies()[ei].end(), [&](const index_t x2i) {
            const Eigen::Vector3d x2 = vertices.row(x2i);
            // Perpendicular foot of for x2
            const Eigen::Vector3d p =
                (x1 - x0) * point_edge_closest_point(x2, x0, x1) + x0;
            return (x - p).dot(p - x2) > 0;
        });
}

bool is_edge_vertex_feasible(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const EdgeVertexCandidate& candidate,
    PointEdgeDistanceType dtype)
{
    const auto& [vi, ei] = candidate;
    const auto [_, e0i, e1i, __] =
        candidate.vertex_ids(mesh.edges(), mesh.faces());

    if (dtype == PointEdgeDistanceType::AUTO) {
        const auto [v, e0, e1, ___] =
            candidate.vertices(vertices, mesh.edges(), mesh.faces());
        dtype = point_edge_distance_type(v, e0, e1);
    }

    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        return ogc::check_vertex_feasible_region(mesh, vertices, vi, e0i)
            && ogc::check_vertex_feasible_region(mesh, vertices, e0i, vi);

    case PointEdgeDistanceType::P_E1:
        return ogc::check_vertex_feasible_region(mesh, vertices, vi, e1i)
            && ogc::check_vertex_feasible_region(mesh, vertices, e1i, vi);

    case PointEdgeDistanceType::P_E:
        // vi must be in the feasible region in this case
        return true;

    case PointEdgeDistanceType::AUTO:
    default:
        throw std::invalid_argument(
            fmt::format(
                "Invalid PointEdgeDistanceType: {}", static_cast<int>(dtype)));
    }
}

bool is_edge_edge_feasible(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const EdgeEdgeCandidate& candidate,
    EdgeEdgeDistanceType dtype)
{
    const auto [ea0i, ea1i, eb0i, eb1i] =
        candidate.vertex_ids(mesh.edges(), mesh.faces());

    if (dtype == EdgeEdgeDistanceType::AUTO) {
        const auto [ea0, ea1, eb0, eb1] =
            candidate.vertices(vertices, mesh.edges(), mesh.faces());
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return ogc::check_vertex_feasible_region(mesh, vertices, ea0i, eb0i)
            && ogc::check_vertex_feasible_region(mesh, vertices, eb0i, ea0i);

    case EdgeEdgeDistanceType::EA0_EB1:
        return ogc::check_vertex_feasible_region(mesh, vertices, ea0i, eb1i)
            && ogc::check_vertex_feasible_region(mesh, vertices, eb1i, ea0i);

    case EdgeEdgeDistanceType::EA1_EB0:
        return ogc::check_vertex_feasible_region(mesh, vertices, ea1i, eb0i)
            && ogc::check_vertex_feasible_region(mesh, vertices, eb0i, ea1i);

    case EdgeEdgeDistanceType::EA1_EB1:
        return ogc::check_vertex_feasible_region(mesh, vertices, ea1i, eb1i)
            && ogc::check_vertex_feasible_region(mesh, vertices, eb1i, ea1i);

    case EdgeEdgeDistanceType::EA_EB0: {
        const double alpha = point_edge_closest_point(
            vertices.row(eb0i), vertices.row(ea0i), vertices.row(ea1i));
        const Eigen::Vector3d xc =
            (vertices.row(ea1i) - vertices.row(ea0i)) * alpha
            + vertices.row(ea0i);
        return ogc::check_vertex_feasible_region(mesh, vertices, xc, eb1i);
    }

    case EdgeEdgeDistanceType::EA_EB1: {
        const double alpha = point_edge_closest_point(
            vertices.row(eb1i), vertices.row(ea0i), vertices.row(ea1i));
        const Eigen::Vector3d xc =
            (vertices.row(ea1i) - vertices.row(ea0i)) * alpha
            + vertices.row(ea0i);
        return ogc::check_vertex_feasible_region(mesh, vertices, xc, eb1i);
    }

    case EdgeEdgeDistanceType::EA0_EB: {
        const double alpha = point_edge_closest_point(
            vertices.row(eb0i), vertices.row(ea0i), vertices.row(ea1i));
        const Eigen::Vector3d xc =
            (vertices.row(eb1i) - vertices.row(eb0i)) * alpha
            + vertices.row(eb0i);
        return ogc::check_vertex_feasible_region(mesh, vertices, xc, ea0i);
    }

    case EdgeEdgeDistanceType::EA1_EB: {
        const double alpha = point_edge_closest_point(
            vertices.row(eb1i), vertices.row(ea0i), vertices.row(ea1i));
        const Eigen::Vector3d xc =
            (vertices.row(eb1i) - vertices.row(eb0i)) * alpha
            + vertices.row(eb0i);
        return ogc::check_vertex_feasible_region(mesh, vertices, xc, ea1i);
    }

    case EdgeEdgeDistanceType::EA_EB:
        return true;

    case EdgeEdgeDistanceType::AUTO:
    default:
        throw std::invalid_argument(
            fmt::format(
                "Invalid EdgeEdgeDistanceType: {}", static_cast<int>(dtype)));
    }
}

bool is_face_vertex_feasible(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const FaceVertexCandidate& candidate,
    PointTriangleDistanceType dtype)
{
    const auto& [fi, vi] = candidate;
    const auto [_, f0i, f1i, f2i] =
        candidate.vertex_ids(mesh.edges(), mesh.faces());

    if (dtype == PointTriangleDistanceType::AUTO) {
        const auto [v, f0, f1, f2] =
            candidate.vertices(vertices, mesh.edges(), mesh.faces());
        dtype = point_triangle_distance_type(v, f0, f1, f2);
    }

    switch (dtype) {
    case PointTriangleDistanceType::P_T0:
        return ogc::check_vertex_feasible_region(mesh, vertices, vi, f0i)
            && ogc::check_vertex_feasible_region(mesh, vertices, f0i, vi);

    case PointTriangleDistanceType::P_T1:
        return ogc::check_vertex_feasible_region(mesh, vertices, vi, f1i)
            && ogc::check_vertex_feasible_region(mesh, vertices, f1i, vi);

    case PointTriangleDistanceType::P_T2:
        return ogc::check_vertex_feasible_region(mesh, vertices, vi, f2i)
            && ogc::check_vertex_feasible_region(mesh, vertices, f2i, vi);

    case PointTriangleDistanceType::P_E0:
        return ogc::check_edge_feasible_region(
            mesh, vertices, vi, mesh.faces_to_edges()(fi, 0));

    case PointTriangleDistanceType::P_E1:
        return ogc::check_edge_feasible_region(
            mesh, vertices, vi, mesh.faces_to_edges()(fi, 1));

    case PointTriangleDistanceType::P_E2:
        return ogc::check_edge_feasible_region(
            mesh, vertices, vi, mesh.faces_to_edges()(fi, 2));

    case PointTriangleDistanceType::P_T:
        // vi must be in the feasible region in this case
        return true;

    case PointTriangleDistanceType::AUTO:
    default:
        throw std::invalid_argument(
            fmt::format(
                "Invalid PointTriangleDistanceType: {}",
                static_cast<int>(dtype)));
    }
}

} // namespace ipc::ogc