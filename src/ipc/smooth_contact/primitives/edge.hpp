#pragma once

#include "primitive.hpp"
#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
    class Edge3 : public Primitive
    {
    public:
        constexpr static int n_core_points = 2;
        constexpr static int dim = 3;
        // d is a vector from closest point on the edge to the point outside of the edge
        Edge3(const long &id, 
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const double &alpha,
        const double &beta);
        
        int n_vertices() const override;
        int n_dofs() const override { return n_vertices() * 3; }

        double potential(const Vector3d &d, const Vector12d &x) const;
        Vector15d grad(const Vector3d &d, const Vector12d &x) const;
        Matrix15d hessian(const Vector3d &d, const Vector12d &x) const;
    private:
        ORIENTATION_TYPES otypes;
    };

    template <typename scalar>
    inline scalar smooth_edge2_term(
        const Eigen::Ref<const Vector2<scalar>>& direc,
        const Eigen::Ref<const Vector2<scalar>>& tangent)
    {
        return tangent.norm();
    }

    /// @brief 
    /// @tparam scalar 
    /// @param direc from edge to point outside, normalized
    /// @param e0 
    /// @param e1 
    /// @param f0 face [f0, e0, e1]
    /// @param f1 face [f1, e1, e0]
    /// @param alpha 
    /// @return 
    template <typename scalar>
    inline scalar smooth_edge3_term(
        const Eigen::Ref<const Vector3<scalar>>& direc,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        scalar tangent_term = scalar(1.); 
        // if (otypes.tangent_type(0) != HEAVISIDE_TYPE::ONE)
        {
            const Vector3<scalar> t0 = point_line_closest_point_direction<scalar>(f0, e0, e1);
            tangent_term = tangent_term * Math<scalar>::smooth_heaviside(-direc.dot(t0) / t0.norm(), alpha, beta);
        }
        // if (otypes.tangent_type(1) != HEAVISIDE_TYPE::ONE)
        {
            const Vector3<scalar> t1 = point_line_closest_point_direction<scalar>(f1, e0, e1);
            tangent_term = tangent_term * Math<scalar>::smooth_heaviside(-direc.dot(t1) / t1.norm(), alpha, beta);
        }

        scalar normal_term = scalar(0.);
        // if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE || otypes.normal_type(1) == HEAVISIDE_TYPE::ONE)
            // normal_term = scalar(1.);
        // else
        {
            const Vector3<scalar> n0 = (e0 - f0).cross(e1 - f0);
            const Vector3<scalar> n1 = -(e0 - f1).cross(e1 - f1);
            normal_term = Math<scalar>::smooth_heaviside( (Math<scalar>::smooth_heaviside(direc.dot(n0) / n0.norm(), alpha, beta) +
            Math<scalar>::smooth_heaviside(direc.dot(n1) / n1.norm(), alpha, beta) - 1), alpha, 0);
        }

        return (e1 - e0).squaredNorm() * tangent_term * normal_term;
    }

    inline bool smooth_edge3_term_type(
        const Eigen::Ref<const Vector3<double>>& direc,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        ORIENTATION_TYPES &otypes)
    {
        otypes.set_size(2);

        const Vector3<double> t0 = point_line_closest_point_direction<double>(f0, e0, e1);
        const Vector3<double> t1 = point_line_closest_point_direction<double>(f1, e0, e1);
        otypes.tangent_type(0) = otypes.compute_type(-direc.dot(t0) / t0.norm(), alpha, beta);
        otypes.tangent_type(1) = otypes.compute_type(-direc.dot(t1) / t1.norm(), alpha, beta);
        if (otypes.tangent_type(0) == HEAVISIDE_TYPE::ZERO || otypes.tangent_type(1) == HEAVISIDE_TYPE::ZERO)
            return false;

        const Vector3<double> n0 = (e0 - f0).cross(e1 - f0);
        const Vector3<double> n1 = -(e0 - f1).cross(e1 - f1);
        const double tmp0 = direc.dot(n0) / n0.norm();
        const double tmp1 = direc.dot(n1) / n1.norm();
        otypes.normal_type(0) = otypes.compute_type(tmp0, alpha, beta);
        otypes.normal_type(1) = otypes.compute_type(tmp1, alpha, beta);
        const double sum = Math<double>::smooth_heaviside(tmp0, alpha, beta) + 
                            Math<double>::smooth_heaviside(tmp1, alpha, beta);
        if (sum <= 1 - alpha)
            return false;
        else if (sum >= 1)
        {
            otypes.normal_type(0) = HEAVISIDE_TYPE::ONE;
            otypes.normal_type(1) = HEAVISIDE_TYPE::ONE;
        }

        return true;
    }
}