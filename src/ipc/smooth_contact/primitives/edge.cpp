#include "edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    namespace {
        template <typename T>
        std::array<Vector3<T>, 4> double_to_autodiff(
            const Eigen::Ref<const Eigen::Vector3d>& v0,
            const Eigen::Ref<const Eigen::Vector3d>& v1,
            const Eigen::Ref<const Eigen::Vector3d>& f0,
            const Eigen::Ref<const Eigen::Vector3d>& f1)
        {
            Vector3<T> v0_, v1_, f0_, f1_;
            for (int d = 0; d < 3; d++)
            {
                v0_(d) = T(  d, v0(d));
                v1_(d) = T(3+d, v1(d));
                f0_(d) = T(6+d, f0(d));
                f1_(d) = T(9+d, f1(d));
            }
            return {{v0_, v1_, f0_, f1_}};
        }
    }
    // d is a vector from any point on the edge to the point outside of the edge
    Edge3::Edge3(const long &eid,
        const Eigen::Ref<const Eigen::Vector3d>& d,
        const Eigen::Ref<const Eigen::Vector3d>& v0,
        const Eigen::Ref<const Eigen::Vector3d>& v1,
        const Eigen::Ref<const Eigen::Vector3d>& f0,
        const Eigen::Ref<const Eigen::Vector3d>& f1,
        const double alpha, const double beta)
    : Primitive(eid), _d(d), _v0(v0), _v1(v1), _f0(f0), _f1(f1), _alpha(alpha), _beta(beta)
    {
        otypes.set_size(2);
    }
    
    bool Edge3::is_active()
    {
        return smooth_edge3_term_type(_d, _v0, _v1, _f0, _f1, _alpha, _beta, otypes);
    }
    int Edge3::n_vertices() const
    {
        return 1 + n_edge_neighbors_3d;
    }
    
    double Edge3::potential() const
    {
        return smooth_edge3_term<double>(_d, _v0, _v1, _f0, _f1, _alpha, _beta, otypes);
    }
    Eigen::VectorXd Edge3::grad() const
    {
        Eigen::VectorXd tmp(15);
        tmp << _d, _v0, _v1, _f0, _f1;
        auto X = slice_positions<ADGrad<15>, 5, 3>(tmp);
        return smooth_edge3_term<ADGrad<15>>(X.row(0), X.row(1), X.row(2), X.row(3), X.row(4), _alpha, _beta, otypes).getGradient();
    }
    Eigen::MatrixXd Edge3::hessian() const
    {
        Eigen::MatrixXd h;
        Eigen::VectorXd tmp(15);
        tmp << _d, _v0, _v1, _f0, _f1;
        auto X = slice_positions<ADHessian<15>, 5, 3>(tmp);
        return smooth_edge3_term<ADHessian<15>>(X.row(0), X.row(1), X.row(2), X.row(3), X.row(4), _alpha, _beta, otypes).getHessian();
    }
}
