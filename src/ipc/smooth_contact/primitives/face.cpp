#include "face.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    // namespace {
    //     template <typename T>
    //     std::array<Vector3<T>, 3> double_to_autodiff(
    //         const Eigen::Ref<const Eigen::Vector3d>& v0,
    //         const Eigen::Ref<const Eigen::Vector3d>& v1,
    //         const Eigen::Ref<const Eigen::Vector3d>& v2)
    //     {
    //         Vector3<T> v0_, v1_, v2_;
    //         for (int d = 0; d < 3; d++)
    //         {
    //             v0_(d) = T(  d, v0(d));
    //             v1_(d) = T(3+d, v1(d));
    //             v2_(d) = T(6+d, v2(d));
    //         }
    //         return {{v0_, v1_, v2_}};
    //     }
    // }
    Face::Face(const long &fid,
        const Eigen::Ref<const Eigen::Vector3d>& d,
        const Eigen::Ref<const Eigen::Vector3d>& v0,
        const Eigen::Ref<const Eigen::Vector3d>& v1,
        const Eigen::Ref<const Eigen::Vector3d>& v2)
    : Primitive(fid), _d(d), _v0(v0), _v1(v1), _v2(v2)
    {
        normal = (v1 - v0).cross(v2 - v0);
    }
    bool Face::is_active()
    {
        return normal.dot(_d) > 0;
    }
    int Face::n_vertices() const
    {
        return 1 + n_face_neighbors_3d;
    }
    double Face::potential() const
    {
        return smooth_face_term<double>(_v0, _v1, _v2);
    }
    Eigen::VectorXd Face::grad() const
    {
        Eigen::VectorXd g;
        g.setZero(n_vertices() * 3);
        Eigen::VectorXd tmp(9);
        tmp << _v0, _v1, _v2;
        auto X = slice_positions<ADGrad<9>, 3, 3>(tmp);
        g.tail<9>() = smooth_face_term<ADGrad<9>>(X.row(0), X.row(1), X.row(2)).getGradient();
        return g;
    }
    Eigen::MatrixXd Face::hessian() const
    {
        Eigen::MatrixXd h;
        h.setZero(n_vertices() * 3, n_vertices() * 3);
        Eigen::VectorXd tmp(9);
        tmp << _v0, _v1, _v2;
        auto X = slice_positions<ADHessian<9>, 3, 3>(tmp);
        h.bottomRightCorner<9, 9>() = smooth_face_term<ADHessian<9>>(X.row(0), X.row(1), X.row(2)).getHessian();
        return h;
    }
}
