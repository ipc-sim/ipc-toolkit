#include "ground_contact.hpp"

#include "ipc/utils/eigen_ext.hpp"

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::rigid {

GroundContact::GroundContact(const double height)
    : m_ground_height(height)
    , m_dhat(1e-3)
{
    m_barrier = std::make_shared<ClampedLogBarrier>();
}

// ---- Cumulative functions ---------------------------------------------------

double GroundContact::operator()(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int ndof = x.size() / bodies.num_bodies();

    double energy = 0.0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        energy += operator()(bodies, i, Pose(x.segment(i * ndof, ndof)));
    }
    return energy;
}

Eigen::VectorXd GroundContact::gradient(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        grad.segment(i * ndof, ndof) =
            gradient(bodies, i, Pose(x.segment(i * ndof, ndof)));
    }
    return grad;
}

Eigen::MatrixXd GroundContact::hessian(
    const RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(x.size(), x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        hess.block(i * ndof, i * ndof, ndof, ndof) = hessian(
            bodies, i, Pose(x.segment(i * ndof, ndof)), project_hessian_to_psd);
    }
    return hess;
}

// ---- Per-body functions -----------------------------------------------------

double GroundContact::operator()(
    const RigidBodies& bodies, const size_t body_index, const Pose& pose) const
{
    double energy = 0.0;
    Eigen::MatrixXd vertices = bodies.body_vertices(body_index, pose);
    for (int i = 0; i < vertices.rows(); ++i) {
        const double height_above_ground = vertices(i, 1) - m_ground_height;
        energy += (*m_barrier)(height_above_ground, m_dhat);
    }
    return energy;
}

VectorMax6d GroundContact::gradient(
    const RigidBodies& bodies, const size_t body_index, const Pose& pose) const
{
#ifndef NDEBUG
    const int ndof = pose.position.size() + pose.rotation.size();
#endif

    const Eigen::MatrixXd vertices = bodies.body_vertices(body_index, pose);

    Eigen::VectorXd dE_dV = Eigen::VectorXd::Zero(vertices.size());
    for (int i = 0; i < vertices.rows(); ++i) {
        const double height_above_ground = vertices(i, 1) - m_ground_height;
        const double db_dh =
            m_barrier->first_derivative(height_above_ground, m_dhat);

        // Gradient wrt y position only
        dE_dV(vertices.rows() + i) = db_dh;
    }

    const Eigen::MatrixXd dV_dx = pose.transform_vertices_jacobian(
        bodies.body_rest_positions(body_index));
    assert(dV_dx.rows() == vertices.size() && dV_dx.cols() == ndof);

    return dV_dx.transpose() * dE_dV; // Chain rule
}

MatrixMax6d GroundContact::hessian(
    const RigidBodies& bodies,
    const size_t body_index,
    const Pose& pose,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const int ndof = pose.position.size() + pose.rotation.size();

    const Eigen::MatrixXd vertices = bodies.body_vertices(body_index, pose);

    // Compute dE/dV
    Eigen::VectorXd dE_dV = Eigen::VectorXd::Zero(vertices.size());
    for (int i = 0; i < vertices.rows(); ++i) {
        const double height_above_ground = vertices(i, 1) - m_ground_height;
        const double db_dh =
            m_barrier->first_derivative(height_above_ground, m_dhat);

        // Gradient wrt y position only
        dE_dV(vertices.rows() + i) = db_dh;
    }

    // Compute d2E/dV2
    Eigen::SparseMatrix<double> d2E_dV2(vertices.size(), vertices.size());
    {
        std::vector<Eigen::Triplet<double>> triplets;

        for (int i = 0; i < vertices.rows(); ++i) {
            const double height_above_ground = vertices(i, 1) - m_ground_height;
            const double db_dh =
                m_barrier->second_derivative(height_above_ground, m_dhat);

            if (db_dh != 0.0) {
                // Hessian wrt y position only
                triplets.emplace_back(
                    vertices.rows() + i, vertices.rows() + i, db_dh);
            }
        }

        d2E_dV2.setFromTriplets(triplets.begin(), triplets.end());
    }

    const Eigen::MatrixXd dV_dx = pose.transform_vertices_jacobian(
        bodies.body_rest_positions(body_index));
    assert(dV_dx.rows() == vertices.size() && dV_dx.cols() == ndof);
    const Eigen::MatrixXd d2V_dx2 =
        pose.transform_vertices_hessian(bodies.body_rest_positions(body_index));
    assert(d2V_dx2.rows() == vertices.size());
    assert(d2V_dx2.cols() == ndof * ndof);

    MatrixMax6d hess = dV_dx.transpose() * d2E_dV2 * dV_dx;
    assert(hess.rows() == ndof && hess.cols() == ndof);
    for (int j = 0; j < ndof; ++j) {
        for (int i = 0; i < ndof; ++i) {
            auto d2V_dxj_dxk = d2V_dx2.col(i * ndof + j);
            hess(i, j) += dE_dV.dot(d2V_dxj_dxk);
        }
    }

    return project_to_psd(hess, project_hessian_to_psd);
}

} // namespace ipc::rigid