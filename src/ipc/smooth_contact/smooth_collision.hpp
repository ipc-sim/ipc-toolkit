#pragma once

#include <ipc/collisions/collision.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/utils/unordered_tuple.hpp>

namespace ipc {

constexpr static int max_vert_2d = 6;
constexpr static int max_vert_3d = 24;

template <int nvert>
class SmoothCollision : public Collision<nvert> {
protected:
    SmoothCollision(
        long primitive0_,
        long primitive1_,
        const std::array<double, 2> &dhats_,
        const CollisionMesh &mesh)
    : primitive0(primitive0_), primitive1(primitive1_), dhats(dhats_)
    {
        vertices.fill(-1);
    }
    SmoothCollision(
        long primitive0_,
        long primitive1_,
        const double &dhat,
        const CollisionMesh &mesh)
    : SmoothCollision(primitive0_, primitive1_, {{dhat, dhat}}, mesh)
    { }
public:
    virtual ~SmoothCollision() { }

    bool is_active() const { return is_active_; }

    virtual int ndofs() const = 0;

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

    // template <typename H>
    // friend H AbslHashValue(H h, const SmoothCollision& other)
    // {
    //     long min_ei = std::min(other.primitive0, other.primitive1);
    //     long max_ei = std::max(other.primitive0, other.primitive1);
    //     return H::combine(std::move(h), min_ei, max_ei);
    // }
    unordered_tuple get_hash() const { return unordered_tuple(primitive0, primitive1); }

    double get_dhat(const int &id) const
    {
        return dhats[id];
    }

protected:
    bool is_active_ = true;
    long primitive0, primitive1;
    std::array<double, 2> dhats;
    std::array<long, nvert> vertices;
};

}