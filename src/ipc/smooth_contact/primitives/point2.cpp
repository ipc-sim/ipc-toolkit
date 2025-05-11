#include "point2.hpp"

namespace ipc {

namespace {
    bool smooth_point2_term_type(
        const Eigen::Ref<const Vector2d>& v,
        const Eigen::Ref<const Vector2d>& direc,
        const Eigen::Ref<const Vector2d>& e0,
        const Eigen::Ref<const Vector2d>& e1,
        const ParameterType& param,
        const bool orientable)
    {
        const Vector2d dn = -direc.normalized();
        const Vector2d t0 = (e0 - v).normalized(), t1 = (e1 - v).normalized();
    
        if (dn.dot(t0) <= -param.alpha_t || dn.dot(t1) <= -param.alpha_t)
            return false;
    
        if (orientable)
        {
            const double tmp = Math<double>::smooth_heaviside(-Math<double>::cross2(dn, t0), param.alpha_n, param.beta_n) + 
                                Math<double>::smooth_heaviside(Math<double>::cross2(dn, t1), param.alpha_n, param.beta_n);
            if (tmp <= 1. - param.alpha_n)
                return false;
        }
    
        return true;
    }

    template <class scalar>
    scalar smooth_point2_term(
        const Eigen::Ref<const Vector2<scalar>>& v,
        const Eigen::Ref<const Vector2<scalar>>& direc,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const ParameterType& param,
        const bool orientable)
    {
        const Vector2<scalar> dn = -direc.normalized();
        const Vector2<scalar> t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();
    
        scalar val = Math<scalar>::smooth_heaviside(dn.dot(t0), param.alpha_t, param.beta_t) *
                            Math<scalar>::smooth_heaviside(-dn.dot(t1), param.alpha_t, param.beta_t);
    
        if (orientable)
        {
            const scalar tmp = Math<scalar>::smooth_heaviside(-Math<scalar>::cross2(dn, t0), param.alpha_n, param.beta_n) + 
                                    Math<scalar>::smooth_heaviside(-Math<scalar>::cross2(dn, t1), param.alpha_n, param.beta_n);
            val = val * Math<scalar>::smooth_heaviside(tmp - 1., param.alpha_n, 0);
        }
    
        return val * ((e0 - v).norm() + (e1 - v).norm()) / 2.;
    }

    template <class scalar>
    scalar smooth_point2_term_one_side(
        const Eigen::Ref<const Vector2<scalar>>& v,
        const Eigen::Ref<const Vector2<scalar>>& direc,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const ParameterType& param)
    {
        const Vector2<scalar> dn = -direc.normalized();
        const Vector2<scalar> t0 = e0 - v;
    
        const scalar tangent_term = Math<scalar>::smooth_heaviside(dn.dot(t0) / t0.norm(), param.alpha_t, param.beta_t);
        
        return tangent_term * t0.norm();
    }
}

Point2::Point2(const long &id, 
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const ParameterType& param)
: Primitive(id, param)
{
    orientable = !mesh.is_codim_vertex(id);
    auto neighbor_verts = mesh.find_vertex_adjacent_vertices(id);
    has_neighbor_1 = neighbor_verts[0] >= 0;
    has_neighbor_2 = neighbor_verts[1] >= 0;

    if (has_neighbor_1 && has_neighbor_2)
    {
        _vert_ids = {{id, neighbor_verts[0], neighbor_verts[1]}};
        is_active_ = smooth_point2_term_type(vertices.row(id), d, vertices.row(_vert_ids[1]), vertices.row(_vert_ids[2]), _param, orientable);
    }
    else if (has_neighbor_1 || has_neighbor_2)
    {
        _vert_ids = {{id, has_neighbor_1 ? neighbor_verts[0] : neighbor_verts[1]}};

        const Vector2d dn = -d.normalized();
        const Vector2d t0 = (vertices.row(_vert_ids[1]) - vertices.row(_vert_ids[0])).normalized();
    
        is_active_ = dn.dot(t0) > -param.alpha_t;
    }
    else
    {
        _vert_ids = {{id}};
        is_active_ = true;
    }
}

int Point2::n_vertices() const
{
    return _vert_ids.size();
}

double Point2::potential(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
{
    if (has_neighbor_1 && has_neighbor_2)
        return smooth_point2_term<double>(x.segment<dim>(0), d, x.segment<dim>(dim), x.segment<dim>(2 * dim), _param, orientable);
    else if (has_neighbor_1 || has_neighbor_2)
        return smooth_point2_term_one_side<double>(x.segment<dim>(0), d, x.segment<dim>(dim), _param);
    else
        return 1.;
}
Vector<double, -1, Point2::max_size+Point2::dim> Point2::grad(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
{
    if (has_neighbor_1 && has_neighbor_2)
    {
        DiffScalarBase::setVariableCount(4*dim);
        using T = ADGrad<4*dim>;
        Vector<double, 4*dim> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 4, dim> X = slice_positions<T, 4, dim>(tmp);
        return smooth_point2_term<T>(X.row(1), X.row(0), X.row(2), X.row(3), _param, orientable).getGradient();
    }
    else if (has_neighbor_1 || has_neighbor_2)
    {
        DiffScalarBase::setVariableCount(3*dim);
        using T = ADGrad<3*dim>;
        Vector<double, 3*dim> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 3, dim> X = slice_positions<T, 3, dim>(tmp);
        return smooth_point2_term_one_side<T>(X.row(1), X.row(0), X.row(2), _param).getGradient();
    }
    else
        return Vector<double, -1, Point2::max_size+Point2::dim>::Zero(x.size() + d.size());
}
MatrixMax<double, Point2::max_size+Point2::dim, Point2::max_size+Point2::dim> Point2::hessian(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
{
    if (has_neighbor_1 && has_neighbor_2)
    {
        DiffScalarBase::setVariableCount(4*dim);
        using T = ADHessian<4*dim>;
        Vector<double, 4*dim> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 4, dim> X = slice_positions<T, 4, dim>(tmp);
        return smooth_point2_term<T>(X.row(1), X.row(0), X.row(2), X.row(3), _param, orientable).getHessian();
    }
    else if (has_neighbor_1 || has_neighbor_2)
    {
        DiffScalarBase::setVariableCount(3*dim);
        using T = ADHessian<3*dim>;
        Vector<double, 3*dim> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 3, dim> X = slice_positions<T, 3, dim>(tmp);
        return smooth_point2_term_one_side<T>(X.row(1), X.row(0), X.row(2), _param).getHessian();
    }
    else
        return MatrixMax<double, -1, Point2::max_size+Point2::dim>::Zero(x.size() + d.size(), x.size() + d.size());
}
}