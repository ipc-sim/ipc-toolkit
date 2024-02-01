#include "point.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    Point3::Point3(const long &id, 
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const double &alpha,
        const double &beta)
    : Primitive(id, alpha, beta)
    {
        auto neighbor_ids = mesh.find_vertex_adjacent_vertices(id);
        n_neighbors = neighbor_ids.size();
        _vert_ids.reserve(1 + n_neighbors);
        _vert_ids.push_back(id);
        _vert_ids.insert( _vert_ids.end(), neighbor_ids.begin(), neighbor_ids.end() );

        if (_vert_ids.size() > n_vert_neighbors_3d)
            logger().error("Too many neighbors for point3 primitive! {} > {}", _vert_ids.size(), n_vert_neighbors_3d);
        
        otypes.set_size(n_neighbors);

        Eigen::Vector3d v = vertices.row(id);
        MatrixMax<double, n_vert_neighbors_3d, 3> neighbors(n_neighbors, dim);
        int k = 0;
        for (long i : neighbor_ids)
            neighbors.row(k++) = vertices.row(i);

        ORIENTATION_TYPES tmp;
        is_active_ = smooth_point3_term_type(v, -d.normalized(), neighbors, _alpha, _beta, tmp);
    }

    int Point3::n_vertices() const
    {
        return n_neighbors + 1;
    }
    
    double Point3::potential(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
    {
        auto X = slice_positions<double, -1, dim>(x);
        return smooth_point3_term<double>(X.row(0), -d.normalized(), X.bottomRows(n_neighbors), _alpha, _beta, otypes);
    }
    Vector<double, -1, Point3::max_size+Point3::dim> Point3::grad(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
    {
        Vector<double, -1, max_size+dim> tmp(d.size() + x.size());
        tmp << d, x;
        DiffScalarBase::setVariableCount(tmp.size());
        auto X = slice_positions<ADGrad<-1, max_size+dim>, -1, dim>(tmp);
        return smooth_point3_term<ADGrad<-1, max_size+dim>>(X.row(1), -X.row(0) / X.row(0).norm(), X.bottomRows(n_neighbors), _alpha, _beta, otypes).getGradient();
    }
    MatrixMax<double, Point3::max_size+Point3::dim, Point3::max_size+Point3::dim> Point3::hessian(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
    {
        Vector<double, -1, max_size+dim> tmp(d.size() + x.size());
        tmp << d, x;
        DiffScalarBase::setVariableCount(tmp.size());
        auto X = slice_positions<ADHessian<-1, max_size+dim>, -1, dim>(tmp);
        return smooth_point3_term<ADHessian<-1, max_size+dim>>(X.row(1), -X.row(0) / X.row(0).norm(), X.bottomRows(n_neighbors), _alpha, _beta, otypes).getHessian();
    }
}
