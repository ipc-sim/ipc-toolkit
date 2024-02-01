#pragma once
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/smooth_contact/common.hpp>

namespace ipc {
    class Primitive
    {
    public:
        Primitive(const long &id, 
        const double &alpha,
        const double &beta)
        : _id(id), _alpha(alpha), _beta(beta) { }
        virtual ~Primitive() = default;

        bool operator==(const Primitive& other) const
        {
            return _id == other._id;
        }

        long id() const { return _id; }
        bool is_active() const { return is_active_; }
        virtual int n_vertices() const = 0;
        virtual int n_dofs() const = 0;
        const std::vector<long> &vertex_ids() const { return _vert_ids; }

    protected:
        std::vector<long> _vert_ids;
        long _id;
        double _alpha, _beta;

        bool is_active_ = true;
    };
}
