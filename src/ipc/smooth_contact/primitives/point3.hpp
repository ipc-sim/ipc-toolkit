#pragma once

#include "primitive.hpp"
#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
    class Point3 : public Primitive
    {
    public:
        constexpr static int n_core_points = 1;
        constexpr static int dim = 3;
        constexpr static int max_size = n_vert_neighbors_3d * dim;
        // d is a vector from this point to the other primitive
        Point3(const long &id, 
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const double &alpha,
        const double &beta);

        Point3(const long &id, const CollisionMesh& mesh, const Eigen::MatrixXd& vertices);

        int n_vertices() const override;
        int n_dofs() const override { return n_vertices() * dim; }
        
        // assume the following functions are only called if active
        double potential(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const;
        // derivatives including wrt. d (the closest direction) in front
        Vector<double, -1, max_size+dim> grad(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const;
        MatrixMax<double, max_size+dim, max_size+dim> hessian(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const;
    private:
        int n_neighbors;
        ORIENTATION_TYPES otypes;
    };

/// @brief 
/// @tparam scalar 
/// @param direc normalized
/// @param v 
/// @param direc points from v to the other point
/// @param neighbors follow counter-clockwise order
/// @param params 
/// @return 
template <typename scalar, int n_neighbors = -1>
scalar smooth_point3_term(
    const Eigen::Ref<const RowVector3<scalar>>& v,
    const Eigen::Ref<const RowVector3<scalar>>& direc,
    const Eigen::Ref<const Eigen::Matrix<scalar, n_neighbors, 3>>& neighbors,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

bool smooth_point3_term_type(
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Matrix<double, -1, 3> &neighbors,
    const double &alpha,
    const double &beta,
    ORIENTATION_TYPES &otypes);

std::tuple<double, Eigen::VectorXd> smooth_point3_term_gradient(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& neighbors,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

std::tuple<double, Eigen::VectorXd, Eigen::MatrixXd> smooth_point3_term_hessian(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& neighbors,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

double smooth_point3_term_tangent(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

std::tuple<double, Eigen::VectorXd> smooth_point3_term_tangent_gradient(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

std::tuple<double, Eigen::VectorXd, Eigen::MatrixXd> smooth_point3_term_tangent_hessian(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

double smooth_point3_term_normal(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

std::tuple<double, Eigen::VectorXd> smooth_point3_term_normal_gradient(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);

std::tuple<double, Eigen::VectorXd, Eigen::MatrixXd> smooth_point3_term_normal_hessian(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes);
}
