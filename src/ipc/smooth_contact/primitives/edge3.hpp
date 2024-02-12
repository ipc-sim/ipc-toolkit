#pragma once

#include "primitive.hpp"
#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Edge3 : public Primitive {
public:
    constexpr static int n_core_points = 2;
    constexpr static int dim = 3;
    // d is a vector from closest point on the edge to the point outside of the
    // edge
    Edge3(
        const long& id,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const double& alpha,
        const double& beta);

    int n_vertices() const override;
    int n_dofs() const override { return n_vertices() * 3; }

    double potential(
        const Eigen::Ref<const Eigen::Vector3d>& d,
        const Eigen::Ref<const Vector12d>& x) const;
    Vector15d grad(
        const Eigen::Ref<const Eigen::Vector3d>& d,
        const Eigen::Ref<const Vector12d>& x) const;
    Matrix15d hessian(
        const Eigen::Ref<const Eigen::Vector3d>& d,
        const Eigen::Ref<const Vector12d>& x) const;

private:
    ORIENTATION_TYPES otypes;
};

/// @brief
/// @tparam scalar
/// @param dn from edge to point outside, normalized
/// @param e0
/// @param e1
/// @param f0 face [f0, e0, e1]
/// @param f1 face [f1, e1, e0]
/// @param alpha
/// @return
template <typename scalar>
inline scalar smooth_edge3_term_template(
    const Eigen::Ref<const Vector3<scalar>>& dn,
    const Eigen::Ref<const Vector3<scalar>>& e0,
    const Eigen::Ref<const Vector3<scalar>>& e1,
    const Eigen::Ref<const Vector3<scalar>>& f0,
    const Eigen::Ref<const Vector3<scalar>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

bool smooth_edge3_term_type(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    ORIENTATION_TYPES& otypes);

/// @brief
/// @param dn should be unit size, from edge to point outside
/// @param e0
/// @param e1
/// @param f0
/// @param f1
/// @param alpha
/// @param beta
/// @param otypes
/// @return
double smooth_edge3_normal_term(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

std::tuple<double, Vector<double, 15>> smooth_edge3_normal_term_gradient(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

std::tuple<double, Vector<double, 15>, Eigen::Matrix<double, 15, 15>>
smooth_edge3_normal_term_hessian(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

double smooth_edge3_tangent_term(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

std::tuple<double, Vector<double, 15>> smooth_edge3_tangent_term_gradient(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

std::tuple<double, Vector<double, 15>, Eigen::Matrix<double, 15, 15>>
smooth_edge3_tangent_term_hessian(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

double smooth_edge3_term(
    const Eigen::Ref<const Vector3<double>>& direc,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

std::tuple<double, Vector<double, 15>> smooth_edge3_term_gradient(
    const Eigen::Ref<const Vector3<double>>& direc,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);

std::tuple<double, Vector<double, 15>, Eigen::Matrix<double, 15, 15>>
smooth_edge3_term_hessian(
    const Eigen::Ref<const Vector3<double>>& direc,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes);
} // namespace ipc