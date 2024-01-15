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

    bool is_active() const { return is_active_; }

    std::array<long, nvert> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return vertices;
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

    bool operator==(const SmoothCollision& other) const
    {
        return (primitive0 == other.primitive0 && primitive1 == other.primitive1) ||
            (primitive1 == other.primitive0 && primitive0 == other.primitive1);
    }
    bool operator!=(const SmoothCollision& other) const
    {
        return !(*this == other);
    }
    bool operator<(const SmoothCollision& other) const
    {
        long this_min = std::min(this->primitive0, this->primitive1);
        long other_min = std::min(other.primitive0, other.primitive1);
        if (this_min == other_min) {
            return std::max(this->primitive0, this->primitive1)
                < std::max(other.primitive0, other.primitive1);
        }
        return this_min < other_min;
    }
    const long& operator[](int idx) const
    {
        if (idx == 0)
            return primitive0;
        else if (idx == 1)
            return primitive1;
        else
            throw std::runtime_error("Invalid index in smooth_collision!");
    }

    template <typename H>
    friend H AbslHashValue(H h, const SmoothCollision& other)
    {
        long min_ei = std::min(other.primitive0, other.primitive1);
        long max_ei = std::max(other.primitive0, other.primitive1);
        return H::combine(std::move(h), min_ei, max_ei);
    }

    void set_adaptive_dhat(const double &dhat0_, const double &dhat1_)
    {
        // dhat0 = dhat0_;
        // dhat1 = dhat1_;
    }

protected:
    bool is_active_ = true;
    long primitive0, primitive1;
    double dhat0 = 0, dhat1 = 0;
    std::array<long, nvert> vertices;
};

}