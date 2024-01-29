#include "point.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    Point3::Point3(const long &pid,
        const Eigen::Ref<const Vector3<double>>& d,
        const Eigen::Ref<const Vector3<double>>& v,
        const Eigen::Ref<const Eigen::Matrix<double, 3, -1>>& neighbors,
        const double alpha, const double beta)
    : Primitive(pid), _d(d), _v(v), _neighbors(neighbors), _alpha(alpha), _beta(beta)
    {
        otypes.set_size(n_neighbors());
    }
    
    bool Point3::is_active()
    {
        return smooth_point3_term_type(_v, -_d, _neighbors, _alpha, _beta, otypes);
    }
    int Point3::n_vertices() const
    {
        return n_neighbors() + 2;
    }
    
    double Point3::potential() const
    {
        return smooth_point3_term<double>(_v, -_d, _neighbors, _alpha, _beta, otypes);
    }
    Eigen::VectorXd Point3::grad() const
    {
        Eigen::VectorXd tmp(_neighbors.size() + 6);
        tmp.head(3) = _d;
        tmp.segment<3>(3) = _v;
        tmp.tail(_neighbors.size()) = Eigen::Map<const Eigen::VectorXd>(_neighbors.data(), _neighbors.size());
        auto X = slice_positions<ADGrad<-1>, -1, 3>(tmp);
        return smooth_point3_term<ADGrad<-1>>(X.row(1), -X.row(0), X.bottomRows(n_neighbors()), _alpha, _beta, otypes).getGradient();
    }
    Eigen::MatrixXd Point3::hessian() const
    {
        Eigen::VectorXd tmp(_neighbors.size() + 6);
        tmp.head(3) = _d;
        tmp.segment<3>(3) = _v;
        tmp.tail(_neighbors.size()) = Eigen::Map<const Eigen::VectorXd>(_neighbors.data(), _neighbors.size());
        auto X = slice_positions<ADHessian<-1>, -1, 3>(tmp);
        return smooth_point3_term<ADHessian<-1>>(X.row(1), -X.row(0), X.bottomRows(n_neighbors()), _alpha, _beta, otypes).getHessian();
    }
}
