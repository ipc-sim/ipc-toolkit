#include "edge.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
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
        is_active_ = smooth_edge3_term_type(d.normalized(), vertices.row(_vert_ids[0]), vertices.row(_vert_ids[1]), vertices.row(_vert_ids[2]), vertices.row(_vert_ids[3]), _alpha, _beta, otypes);
    }
    
    int Edge3::n_vertices() const
    {
        return n_edge_neighbors_3d;
    }
    
    double Edge3::potential(const Vector3d &d, const Vector12d &x) const
    {
        return smooth_edge3_term_template<double>(d.normalized(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), _alpha, _beta, otypes);
    }
    Vector15d Edge3::grad(const Vector3d &d, const Vector12d &x) const
    {
        // Vector15d tmp;
        // tmp << d, x;
        // DiffScalarBase::setVariableCount(15);
        // auto X = slice_positions<ADGrad<15>, 5, 3>(tmp);
        // return smooth_edge3_term_template<ADGrad<15>>(X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3), X.row(4), _alpha, _beta, otypes).getGradient();

        return std::get<1>(smooth_edge3_term_gradient(d, x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), _alpha, _beta, otypes));
    }
    Matrix15d Edge3::hessian(const Vector3d &d, const Vector12d &x) const
    {
        Vector15d tmp;
        tmp << d, x;
        DiffScalarBase::setVariableCount(15);
        auto X = slice_positions<ADHessian<15>, 5, 3>(tmp);
        return smooth_edge3_term_template<ADHessian<15>>(X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3), X.row(4), _alpha, _beta, otypes).getHessian();
    }

    bool smooth_edge3_term_type(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        ORIENTATION_TYPES &otypes)
    {
        otypes.set_size(2);

        const Vector3<double> t0 = PointEdgeDistance<double, 3>::point_line_closest_point_direction(f0, e0, e1);
        const Vector3<double> t1 = PointEdgeDistance<double, 3>::point_line_closest_point_direction(f1, e0, e1);
        otypes.tangent_type(0) = otypes.compute_type(-dn.dot(t0) / t0.norm(), alpha, beta);
        otypes.tangent_type(1) = otypes.compute_type(-dn.dot(t1) / t1.norm(), alpha, beta);
        if (otypes.tangent_type(0) == HEAVISIDE_TYPE::ZERO || otypes.tangent_type(1) == HEAVISIDE_TYPE::ZERO)
            return false;

        const Vector3<double> n0 = (e0 - f0).cross(e1 - f0);
        const Vector3<double> n1 = -(e0 - f1).cross(e1 - f1);
        const double tmp0 = dn.dot(n0) / n0.norm();
        const double tmp1 = dn.dot(n1) / n1.norm();
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

    double smooth_edge3_normal_term(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE || otypes.normal_type(1) == HEAVISIDE_TYPE::ONE)
            return 1.;
        
        return Math<double>::smooth_heaviside(
            negative_orientation_penalty(e0 - f0, e1 - f0, dn, alpha, beta) +
            negative_orientation_penalty(e0 - f1, e1 - f1,-dn, alpha, beta) - 1, alpha, 0);
    }

    std::tuple<double, Vector<double, 15>>
    smooth_edge3_normal_term_gradient(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        double val = 1.;
        Vector<double, 15> gradient = Vector<double, 15>::Zero();

        if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE || otypes.normal_type(1) == HEAVISIDE_TYPE::ONE)
            return std::make_tuple(val, gradient);
        
        {
            const auto [y, dy] = negative_orientation_penalty_grad(e0 - f0, e1 - f0, dn, alpha, beta);

            val += y;

            // dn
            gradient.head<3>() += dy.tail<3>();
            // e0, e1
            gradient.segment<6>(3) += dy.head<6>();
            // f0
            gradient.segment<3>(9) -= dy.segment<3>(0) + dy.segment<3>(3);
        }

        {
            const auto [y, dy] = negative_orientation_penalty_grad(e0 - f1, e1 - f1, -dn, alpha, beta);

            val += y;

            // dn
            gradient.head<3>() -= dy.tail<3>();
            // e0, e1
            gradient.segment<6>(3) += dy.head<6>();
            // f1
            gradient.segment<3>(12) -= dy.segment<3>(0) + dy.segment<3>(3);
        }

        gradient *= Math<double>::smooth_heaviside_grad(val - 1, alpha, 0);
        val = Math<double>::smooth_heaviside(val - 1, alpha, 0);
        return std::make_tuple(val, gradient);
    }

    std::tuple<double, Vector<double, 15>, Eigen::Matrix<double, 15, 15>>
    smooth_edge3_normal_term_hessian(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        assert(false);
        double value = 1.;
        Vector<double, 15> gradient = Vector<double, 15>::Zero();
        Eigen::Matrix<double, 15, 15> hessian = Eigen::Matrix<double, 15, 15>::Zero();
        return std::make_tuple(value, gradient, hessian);
    }

    double smooth_edge3_tangent_term(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        double tangent_term = 1.;
        if (otypes.tangent_type(0) != HEAVISIDE_TYPE::ONE)
        {
            const Vector3<double> t0 = PointEdgeDistance<double, 3>::point_line_closest_point_direction(f0, e0, e1);
            tangent_term *= opposite_direction_penalty(t0, -dn, alpha, beta);
        }
        if (otypes.tangent_type(1) != HEAVISIDE_TYPE::ONE)
        {
            const Vector3<double> t1 = PointEdgeDistance<double, 3>::point_line_closest_point_direction(f1, e0, e1);
            tangent_term *= opposite_direction_penalty(t1, -dn, alpha, beta);
        }

        return tangent_term;
    }

    std::tuple<double, Vector<double, 15>>
    smooth_edge3_tangent_term_gradient(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        double value = 1.;
        Vector<double, 15> gradient = Vector<double, 15>::Zero();
        if (otypes.tangent_type(0) != HEAVISIDE_TYPE::ONE)
        {
            const auto [t0, g] = PointEdgeDistanceDerivatives<3>::point_line_closest_point_direction_grad(f0, e0, e1);
            const auto [tmp_val, tmp_grad] = opposite_direction_penalty_grad(t0, -dn, alpha, beta);

            value *= tmp_val;
            
            // dn
            gradient.head<3>() = -tmp_grad.tail<3>();
            // e0, e1
            gradient.segment<6>(3) = tmp_grad.head<3>().transpose() * g.block<3, 6>(0, 3);
            // f0
            gradient.segment<3>(9) = tmp_grad.head<3>().transpose() * g.block<3, 3>(0, 0);
        }
        if (otypes.tangent_type(1) != HEAVISIDE_TYPE::ONE)
        {
            const auto [t1, g] = PointEdgeDistanceDerivatives<3>::point_line_closest_point_direction_grad(f1, e0, e1);
            const auto [tmp_val, tmp_grad] = opposite_direction_penalty_grad(t1, -dn, alpha, beta);

            value *= tmp_val;
            
            // dn
            gradient.head<3>() = -tmp_grad.tail<3>();
            // e0, e1
            gradient.segment<6>(3) = tmp_grad.head<3>().transpose() * g.block<3, 6>(0, 3);
            // f0
            gradient.segment<3>(12) = tmp_grad.head<3>().transpose() * g.block<3, 3>(0, 0);
        }

        return std::make_tuple(value, gradient);
    }

    std::tuple<double, Vector<double, 15>, Eigen::Matrix<double, 15, 15>>
    smooth_edge3_tangent_term_hessian(
        const Eigen::Ref<const Vector3<double>>& dn,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        assert(false);
        double value = 1.;
        Vector<double, 15> gradient = Vector<double, 15>::Zero();
        Eigen::Matrix<double, 15, 15> hessian = Eigen::Matrix<double, 15, 15>::Zero();
        return std::make_tuple(value, gradient, hessian);
    }

    double smooth_edge3_term(
        const Eigen::Ref<const Vector3<double>>& direc,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        const Vector3<double> dn = direc.normalized();
        double tangent_term = smooth_edge3_tangent_term(dn, e0, e1, f0, f1, alpha, beta, otypes);
        double normal_term = smooth_edge3_normal_term(dn, e0, e1, f0, f1, alpha, beta, otypes);

        return (e1 - e0).squaredNorm() * tangent_term * normal_term;
    }

    std::tuple<double, Vector<double, 15>>
    smooth_edge3_term_gradient(
        const Eigen::Ref<const Vector3<double>>& direc,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        assert(otypes.size() == 2);

        const auto [dn, dn_grad] = normalize_vector_grad(direc);

        auto [t_term, t_grad] = smooth_edge3_tangent_term_gradient(dn, e0, e1, f0, f1, alpha, beta, otypes);
        auto [n_term, n_grad] = smooth_edge3_normal_term_gradient(dn, e0, e1, f0, f1, alpha, beta, otypes);

        t_grad.head<3>() = dn_grad * t_grad.head<3>();
        n_grad.head<3>() = dn_grad * n_grad.head<3>();

        const double weight = (e1 - e0).squaredNorm();
        const double val = weight * t_term * n_term;
        
        Vector<double, 15> gradient = (weight * t_term) * n_grad + (weight * n_term) * t_grad;
        gradient.segment<3>(3) += (2 * t_term * n_term) * (e0 - e1);
        gradient.segment<3>(6) += (2 * t_term * n_term) * (e1 - e0);

        return std::make_tuple(val, gradient);
    }

    std::tuple<double, Vector<double, 15>, Eigen::Matrix<double, 15, 15>>
    smooth_edge3_term_hessian(
        const Eigen::Ref<const Vector3<double>>& direc,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta,
        const ORIENTATION_TYPES &otypes)
    {
        assert(false);
        double value = 1.;
        Vector<double, 15> gradient = Vector<double, 15>::Zero();
        Eigen::Matrix<double, 15, 15> hessian = Eigen::Matrix<double, 15, 15>::Zero();
        return std::make_tuple(value, gradient, hessian);
    }
}
