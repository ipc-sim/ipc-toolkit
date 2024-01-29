#pragma once
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/smooth_contact/common.hpp>

namespace ipc {
    class Primitive
    {
    public:
        Primitive(const long &id)
        : _id(id)
        {
        }
        virtual ~Primitive() = default;

        virtual bool is_active() = 0;
        virtual int n_vertices() const = 0;

        // assume the following functions are only called if active
        virtual double potential() const = 0;
        virtual Eigen::VectorXd grad() const = 0;
        virtual Eigen::MatrixXd hessian() const = 0;

    protected:
        long _id;
    };
}
