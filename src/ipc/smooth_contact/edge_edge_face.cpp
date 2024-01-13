#include "edge_edge_face.hpp"
#include "smooth_point_face.hpp"
#include "smooth_edge_edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>
#include <iterator>
#include <ipc/utils/logger.hpp>

namespace ipc {
    namespace {
        template <class T>
        std::array<Vector3<T>, 8> slice_positions(const Vector<double, 24> &positions)
        {
            std::array<Vector3<T>, 8> points;
            points.fill(Vector3<T>::Zero(3));
            
            for (int i = 0, id = 0; i < 8; i++)
                for (int d = 0; d < 3; d++, id++)
                    if constexpr (std::is_same<T, double>::value)
                        points[i](d) = positions(id);
                    else
                        points[i](d) = T(id, positions(id));

            return points;
        }

        template <typename Iter>
        size_t index_of(Iter first, Iter last, const typename std::iterator_traits<Iter>::value_type& x)
        {
            size_t i = 0;
            while (first != last && *first != x)
            ++first, ++i;
            return i;
        }

        bool not_in(const long &a, const std::array<long, 2> &b)
        {
            return a != b[0] && a != b[1];
        }
    }

    SmoothEdgeEdge3Collision::SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh)
    : SmoothCollision(primitive0_, primitive1_, mesh)
    {
        std::array<long, 4> faces = {{mesh.edges_to_faces()(primitive0, 0), mesh.edges_to_faces()(primitive0, 1),
                  mesh.edges_to_faces()(primitive1, 0), mesh.edges_to_faces()(primitive1, 1)}};
        std::array<std::array<long, 2>, 2> edges = {{ {{mesh.edges()(primitive0, 0), mesh.edges()(primitive0, 1)}}, {{mesh.edges()(primitive1, 0), mesh.edges()(primitive1, 1)}} }};
        vertices.fill(-1);

        face_to_vertex.setConstant(-1);

        for (int j : {0, 1, 2, 3})
        {
            const int le = j / 2;
            int i;
            for (i = 0; i < 3; i++)
            {
                const auto va = mesh.faces()(faces[j], i);

                if (not_in(va, edges[le]))
                {
                    face_to_vertex(j, 0) = 4 + j;
                    vertices[4 + j] = va;
                    break;
                }
            }
            const auto vb = mesh.faces()(faces[j], (i+1) % 3);
            const auto vc = mesh.faces()(faces[j], (i+2) % 3);

            if (j%2 == 0)
            {
                vertices[j + 0] = vb;
                vertices[j + 1] = vc;
            }

            face_to_vertex(j, 1) = index_of(vertices.begin(), vertices.begin() + 4, vb);
            face_to_vertex(j, 2) = index_of(vertices.begin(), vertices.begin() + 4, vc);
            assert(face_to_vertex(j, 1) < 4 && face_to_vertex(j, 2) < 4);
        }
    }

    template <typename scalar> 
    scalar SmoothEdgeEdge3Collision::evaluate_quadrature(const Vector<double, 24>& positions, const ParameterType &params) const
    {
        std::array<Vector3<scalar>, 8> points = slice_positions<scalar>(positions);
        std::array<Vector3<double>, 8> points_double = slice_positions<double>(positions);
        scalar out = scalar(0.);

        const EdgeEdgeDistanceType dtype = edge_edge_distance_type(points_double[0], points_double[1], points_double[2], points_double[3]);

        std::array<Vector3<scalar>, 4> normals;
        for (int i = 0; i < 4; i++)
            normals[i] = (points[face_to_vertex(i, 1)] - points[face_to_vertex(i, 0)]).cross(points[face_to_vertex(i, 2)] - points[face_to_vertex(i, 0)]).normalized();
        
        out += smooth_edge_edge_potential_single_point<scalar>(
            points[face_to_vertex(0, 1)], points[face_to_vertex(0, 2)],
            points[face_to_vertex(2, 1)], points[face_to_vertex(2, 2)],
            normals[0], normals[1], normals[2], normals[3], params, dtype);

        return out;
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