#pragma once

#include <ipc/collisions/collision.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/smooth_contact/common.hpp>

namespace ipc {

template <int nvert>
class SmoothCollision : public Collision<nvert> {
protected:
    SmoothCollision(
        long primitive0_,
        long primitive1_,
        const double &dhat,
        const CollisionMesh &mesh)
    : primitive0(primitive0_), primitive1(primitive1_), dhat_(dhat)
    {
        vertices.fill(-1);
    }
public:
    constexpr static int max_size = Collision<nvert>::max_size;
    virtual ~SmoothCollision() { }

    bool is_active() const { return is_active_; }

    virtual int ndofs() const = 0;

    std::array<long, nvert> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return vertices;
    }

    // ---- non distance type potential ----

    virtual double operator()(
        const Vector<double, -1, max_size>& positions, 
        const ParameterType &params) const { return 0.; }

    virtual Vector<double, -1, max_size> gradient(
        const Vector<double, -1, max_size>& positions, 
        const ParameterType &params) const { return Vector<double, -1, max_size>::Zero(positions.size()); }

    virtual MatrixMax<double, max_size, max_size> hessian(
        const Vector<double, -1, max_size>& positions, 
        const ParameterType &params) const { return MatrixMax<double, max_size, max_size>::Zero(positions.size(), positions.size()); }

    // ---- distance ----

    Vector<double, -1, max_size>
    compute_distance_gradient(const Vector<double, -1, max_size>& positions) const override
    {
        return Vector<double, -1, max_size>::Zero(max_size);
    }

    MatrixMax<double, max_size, max_size>
    compute_distance_hessian(const Vector<double, -1, max_size>& positions) const override
    {
        return MatrixMax<double, max_size, max_size>::Zero(max_size, max_size);
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
    std::pair<long, long> get_hash() const { return std::make_pair(primitive0, primitive1); }

    double get_dhat() const
    {
        return dhat_;
    }

    virtual std::string name() const { return ""; }

protected:
    bool is_active_ = true;
    long primitive0, primitive1;
    double dhat_;
    std::array<long, nvert> vertices;
};

}