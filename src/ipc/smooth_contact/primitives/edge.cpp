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

        // apply chain rule of "d -> d / |d|" on grad g
        // VectorMax3d normalize_vec_jacobian(VectorMax3d d, VectorMax3d g)
        // {
        //     g /= d.norm();
        //     d.normalize();
        //     return g - d * (d.transpose() * g);
        // }
    }
    // d is a vector from any point on the edge to the point outside of the edge
    Edge3::Edge3(const long &id, 
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const double &alpha,
        const double &beta)
    : Primitive(id, alpha, beta)
    {
        auto ids = mesh.find_edge_adjacent_vertices(id);
        _vert_ids = std::vector<long>(ids.begin(), ids.begin() + ids.size());
        otypes.set_size(2);

        ORIENTATION_TYPES tmp;
        is_active_ = smooth_edge3_term_type(d.normalized(), vertices.row(_vert_ids[0]), vertices.row(_vert_ids[1]), vertices.row(_vert_ids[2]), vertices.row(_vert_ids[3]), _alpha, _beta, tmp);
    }
    
    int Edge3::n_vertices() const
    {
        return n_edge_neighbors_3d;
    }
    
    double Edge3::potential(const Vector3d &d, const Vector12d &x) const
    {
        return smooth_edge3_term<double>(d.normalized(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), _alpha, _beta, otypes);
    }
    Vector15d Edge3::grad(const Vector3d &d, const Vector12d &x) const
    {
        Vector15d tmp;
        tmp << d, x;
        DiffScalarBase::setVariableCount(15);
        auto X = slice_positions<ADGrad<15>, 5, 3>(tmp);
        return smooth_edge3_term<ADGrad<15>>(X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3), X.row(4), _alpha, _beta, otypes).getGradient();
    }
    Matrix15d Edge3::hessian(const Vector3d &d, const Vector12d &x) const
    {
        Vector15d tmp;
        tmp << d, x;
        DiffScalarBase::setVariableCount(15);
        auto X = slice_positions<ADHessian<15>, 5, 3>(tmp);
        return smooth_edge3_term<ADHessian<15>>(X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3), X.row(4), _alpha, _beta, otypes).getHessian();
    }
}
