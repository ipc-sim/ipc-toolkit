#pragma once

#include <ipc/dynamics/rigid/pose.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace ipc::rigid {

class ImplicitEuler {
public:
    Eigen::VectorXd x_prev;
    Eigen::VectorXd v_prev;
    Eigen::VectorXd a_prev;
    double dt;
    size_t num_bodies;
    size_t pos_ndof;
    size_t rot_ndof;

    ImplicitEuler(
        Eigen::ConstRef<Eigen::VectorXd> _x,
        Eigen::ConstRef<Eigen::VectorXd> _v,
        Eigen::ConstRef<Eigen::VectorXd> _a,
        const double _dt,
        const size_t _num_bodies)
        : x_prev(_x)
        , v_prev(_v)
        , a_prev(_a)
        , dt(_dt)
        , num_bodies(_num_bodies)
    {
        assert(x_prev.size() == v_prev.size());
        assert(v_prev.size() == a_prev.size());

        pos_ndof = 3, rot_ndof = 9; // Default to 3D case
        if (x_prev.size() != num_bodies * (pos_ndof + rot_ndof)) {
            pos_ndof = 2, rot_ndof = 1; // 2D case
        }
        assert(x_prev.size() == num_bodies * (pos_ndof + rot_ndof));
    }

    void update(Eigen::ConstRef<Eigen::VectorXd> x)
    {
        assert(x.size() == x_prev.size());
        const Eigen::VectorXd v = (x - x_prev) / dt;
        a_prev = (v - v_prev) / dt;
        v_prev = v;
        x_prev = x;
    }

    std::vector<AffinePose> predicted_pose() const
    {
        const Eigen::VectorXd x_hat = x_prev + dt * v_prev;

        std::vector<AffinePose> predicted(num_bodies);
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, num_bodies),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    predicted[i] = pose(x_hat, i);
                }
            });

        return predicted;
    }

    AffinePose pose(size_t i) const { return pose(x_prev, i); }

    AffinePose pose(Eigen::ConstRef<Eigen::VectorXd> x, size_t i) const
    {
        AffinePose pose;

        pose.position = x.segment(i * (pos_ndof + rot_ndof), pos_ndof);
        if (rot_ndof == 1) {
            pose.rotation.resize(1, 1);
            pose.rotation(0, 0) = x(i * (pos_ndof + rot_ndof) + pos_ndof);
        } else {
            pose.rotation =
                x.segment(i * (pos_ndof + rot_ndof) + pos_ndof, rot_ndof)
                    .reshaped(3, 3);
        }

        return pose;
    }
};

} // namespace ipc::rigid