#pragma once

#include <ipc/collisions/collision.hpp>
#include <ipc/collision_mesh.hpp>

namespace ipc {

template <int nvert>
class SmoothCollision : public Collision<nvert> {
public:
    SmoothCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh)
    : primitive0(primitive0_), primitive1(primitive1_)
    {
        vertices.fill(-1);
    }
    virtual ~SmoothCollision() { }

    std::array<long, nvert> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return vertices;
    }

    double compute_distance(const Vector<double, -1, 3*nvert>& positions) const override
    {
        return 0.;
    }

    Vector<double, -1, 3*nvert>
    compute_distance_gradient(const Vector<double, -1, 3*nvert>& positions) const override
    {
        return Vector<double, -1, 3*nvert>::Zero(3*nvert);
    }

    MatrixMax<double, 3*nvert, 3*nvert>
    compute_distance_hessian(const Vector<double, -1, 3*nvert>& positions) const override
    {
        return MatrixMax<double, 3*nvert, 3*nvert>::Zero(3*nvert, 3*nvert);
    }

    virtual void set_adaptive_dhat(const CollisionMesh &mesh, const double &dhat) {}

    long primitive0, primitive1;
    double dhat0 = 0, dhat1 = 0;
    std::array<long, nvert> vertices;
};

}