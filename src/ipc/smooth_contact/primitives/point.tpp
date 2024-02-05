#pragma once
#include "point.hpp"

namespace ipc {
template <class scalar>
scalar smooth_point2_term(
    const Eigen::Ref<const Vector2<scalar>>& v,
    const Eigen::Ref<const Vector2<scalar>>& direc,
    const Eigen::Ref<const Vector2<scalar>>& e0,
    const Eigen::Ref<const Vector2<scalar>>& e1,
    const double &alpha,
    const double &beta)
{
    const Vector2<scalar> t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();

    const scalar tangent_term = Math<scalar>::smooth_heaviside(direc.dot(t0), alpha, beta) *
                        Math<scalar>::smooth_heaviside(-direc.dot(t1), alpha, beta);

    const scalar tmp = Math<scalar>::smooth_heaviside(-Math<scalar>::cross2(direc, t0), alpha, beta) + 
                         Math<scalar>::smooth_heaviside(-Math<scalar>::cross2(direc, t1), alpha, beta);
    const scalar normal_term = Math<scalar>::smooth_heaviside(tmp - 1., alpha, 0);

    return tangent_term * normal_term * ((e0 - v).norm() + (e1 - v).norm()) / 2.;
}

template <typename scalar, int n_neighbors = -1>
scalar smooth_point3_term(
    const Eigen::Ref<const RowVector3<scalar>>& v,
    const Eigen::Ref<const RowVector3<scalar>>& direc,
    const Eigen::Ref<const Eigen::Matrix<scalar, n_neighbors, 3>>& neighbors,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    RowVector3<scalar> t, t_prev;
    assert(neighbors.rows() > 2);
    assert(otypes.size() == neighbors.rows());

    const RowVector3<scalar> dn = direc.normalized();
    scalar tangent_term(1.);
    scalar weight(0.);
    scalar normal_term(0.);
    t_prev = neighbors.row(neighbors.rows() - 1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
            tangent_term = tangent_term * Math<scalar>::smooth_heaviside(-dn.dot(t) / t.norm(), alpha, beta);

        if (otypes.normal_type(a) == HEAVISIDE_TYPE::VARIANT)
            normal_term = normal_term + Math<scalar>::smooth_heaviside(dn.dot(t_prev.cross(t).normalized()), alpha, beta);
        
        weight = weight + t.squaredNorm();
        std::swap(t, t_prev);
    }

    weight /= 3.;
    
    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE)
        normal_term = scalar(1.);
    else
        normal_term = Math<scalar>::smooth_heaviside(normal_term - 1, alpha, 0);

    return weight * normal_term * tangent_term;
}
}
