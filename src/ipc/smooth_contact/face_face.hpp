#pragma once

#include <ipc/candidates/face_face.hpp>
#include <ipc/collisions/collision.hpp>
#include <ipc/collision_mesh.hpp>

namespace ipc {

class SmoothFaceFaceCollision : public FaceFaceCandidate, public Collision<6> {
public:
    SmoothFaceFaceCollision(
        long _face0_id, 
        long _face1_id, 
        const CollisionMesh &mesh)
    : FaceFaceCandidate(_face0_id, _face1_id)
    { 
        vertices = vertex_ids(mesh.edges(), mesh.faces());
    }
    SmoothFaceFaceCollision(
        long _face0_id, 
        long _face1_id, 
        std::array<long, 6> _vertices)
    : FaceFaceCandidate(_face0_id, _face1_id), vertices(_vertices)
    { }
    virtual ~SmoothFaceFaceCollision() { }

    int num_vertices() const override
    {
        return 6;
    }

    std::array<long, 6> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const override
    {
        return {{faces(face0_id, 0), faces(face0_id, 1), faces(face0_id, 2),
                faces(face1_id, 0), faces(face1_id, 1), faces(face1_id, 2)}};
    }
    
    double compute_distance(const Vector<double, -1, 18>& positions) const override;

    Vector<double, -1, 18>
    compute_distance_gradient(const Vector<double, -1, 18>& positions) const override;

    MatrixMax<double, 18, 18>
    compute_distance_hessian(const Vector<double, -1, 18>& positions) const override;

    double operator()(const Vector<double, -1, 18>& positions, 
        const ParameterType &params) const override;

    Vector<double, -1, 18> gradient(
        const Vector<double, -1, 18>& positions, 
        const ParameterType &params) const override;

    MatrixMax<double, 18, 18> hessian(
        const Vector<double, -1, 18>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

    template <typename H>
    friend H AbslHashValue(H h, const SmoothFaceFaceCollision& ff)
    {
        long min_fi = std::min(ff.face0_id, ff.face1_id);
        long max_fi = std::max(ff.face0_id, ff.face1_id);
        return H::combine(std::move(h), min_fi, max_fi, ff.vertices);
    }

    void set_adaptive_dhat(const CollisionMesh &mesh, const double &dhat)
    {
        // dhat0 = std::min(dhat, mesh.min_distance_in_rest_config(edge0_id));
        // dhat1 = std::min(dhat, mesh.min_distance_in_rest_config(edge1_id));
    }

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector<double, -1, 18>& positions, const ParameterType &params) const;

    std::array<long, 6> vertices;

    double dhat0 = 0, dhat1 = 0;
};

}