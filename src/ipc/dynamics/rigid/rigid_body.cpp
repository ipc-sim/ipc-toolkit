#include "rigid_body.hpp"

#include <ipc/dynamics/rigid/mass.hpp>
#include <ipc/utils/sinc.hpp>

namespace ipc::rigid {

namespace {
    void center_vertices(
        Eigen::Ref<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        Pose<>& pose)
    {
        // compute the center of mass several times to get more accurate
        for (int i = 0; i < 10; i++) {
            double mass;
            VectorMax3d com;
            MatrixMax3d inertia;
            compute_mass_properties(
                vertices,
                (vertices.cols() == 2 || faces.size() == 0) ? edges : faces,
                mass, com, inertia);
            vertices.rowwise() -= com.transpose();
            pose.position += com;
            if (com.squaredNorm() < 1e-8) {
                break;
            }
        }
    }

    inline Eigen::DiagonalMatrix<double, 3> compute_J(const VectorMax3d& I)
    {
        if (I.size() == 1) {
            return I.asDiagonal();
        } else {
            assert(I.size() == 3);
            return Eigen::DiagonalMatrix<double, 3>(
                0.5 * (-I.x() + I.y() + I.z()), //
                0.5 * (I.x() - I.y() + I.z()),  //
                0.5 * (I.x() + I.y() - I.z()));
        }
    }
} // namespace

RigidBody::RigidBody(
    Eigen::Ref<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    Pose<>& initial_pose)
{
    const double density = 1.0; // Default density

    assert(vertices.size() > 0);
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    const int dim = vertices.cols();
    assert(dim == 2 || dim == 3);

    // 1. Center the vertices, so the mass properties are computed correctly
    center_vertices(vertices, edges, faces, initial_pose);

    // 2. Compute the mass properties
    VectorMax3d center_of_mass;
    MatrixMax3d inertia_tensor;
    compute_mass_properties(
        vertices, (dim == 2 || faces.size() == 0) ? edges : faces, m_mass,
        center_of_mass, inertia_tensor);

    // 3. The mass above is actually volume, so we need to scale it by the
    // density to get the mass.
    m_mass *= density;

    // 4. Convert the inertia tensor to the principal axes moments of inertia
    MatrixMax3d R0;
    if (dim == 3) {
        // This computation is taken from ProjectChrono: https://bit.ly/2RpbTl1
        // The eigen values of the inertia tensor are the principal moments
        // of inertia, which are the diagonal elements of the diagonalized
        // inertia tensor. The eigenvectors are the principal axes of the
        // inertia tensor, which are the columns of the rotation matrix R₀.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver;

        // Remove small values from the inertia tensor to avoid numerical
        // issues in the eigen decomposition
        const double threshold = 1e-16 * inertia_tensor.maxCoeff();
        inertia_tensor = (inertia_tensor.array().abs() < threshold)
                             .select(0.0, inertia_tensor);

        solver.compute(inertia_tensor);
        assert(solver.info() == Eigen::Success);

        // Multiply by density to get the units of moment of inertia
        m_moment_of_inertia = density * solver.eigenvalues();
        if ((m_moment_of_inertia.array() < 0).any()) {
            logger().warn(
                "Negative moments of inertia ({}), inverting.",
                m_moment_of_inertia);
            // This typically only happens with negative ε inertias
            m_moment_of_inertia = m_moment_of_inertia.array().abs();
        }

        // The rotation from the principal inertial frame to the input world
        // frame.
        R0 = solver.eigenvectors();

        // Ensure that we have an orientation preserving transform
        if (R0.determinant() < 0.0) {
            R0.col(0) *= -1.0;
        }
        assert(R0.isUnitary(1e-9));

        // TODO: Enable this code
        // int num_rot_dof_fixed =
        //     is_dof_fixed.tail(PoseD::dim_to_rot_ndof(dim())).count();
        // if (num_rot_dof_fixed == 2) {
        //     // Convert moment of inertia to world coordinates
        //     // https://physics.stackexchange.com/a/268812
        //     moment_of_inertia = -I.diagonal().array() + I.diagonal().sum();
        //     R0.setIdentity();
        // } else if (num_rot_dof_fixed == 1) {
        //     spdlog::warn(
        //         "Rigid body dynamics with two rotational DoF has "
        //         "not been tested thoroughly.");
        // }

        // TODO: this code below need to be updated

        // Remove the initial rotation from the rest vertices
        vertices = vertices * R0;

        // Store the initial rotation in the pose (R = RᵢR₀)
        Eigen::AngleAxisd r = Eigen::AngleAxisd(
            Eigen::Matrix3d(initial_pose.rotation_matrix() * R0));
        initial_pose.rotation = r.angle() * r.axis();

        // TODO:
        // ω = R₀ᵀω₀ (ω₀ expressed in body coordinates)
        // this->velocity.rotation = R0.transpose() * this->velocity.rotation;
        // Eigen::Matrix3d Q_t0 = this->pose.construct_rotation_matrix();
        // this->Qdot = Q_t0 * Hat(this->velocity.rotation);

        // τ = R₀ᵀτ₀ (τ₀ expressed in body coordinates)
        // NOTE: this transformation will be done later
        // this->force.rotation = R0.transpose() * this->force.rotation;
    } else {
        // For 2D, the inertia tensor is a scalar, and the rotation vector
        // is a single value.
        m_moment_of_inertia = density * inertia_tensor.diagonal();
        // The input orientation is already in the inertial frame
        R0 = Eigen::Matrix<double, 1, 1>::Identity();
    }

    m_J = compute_J(m_moment_of_inertia);

    // Zero out the velocity and forces of fixed dof
    // this->velocity.zero_dof(is_dof_fixed, R0);
    // this->force.zero_dof(is_dof_fixed, R0);

    // Compute and construct some useful constants
    // mass_matrix.resize(ndof());
    // mass_matrix.diagonal().head(pos_ndof()).setConstant(mass);
    // mass_matrix.diagonal().tail(rot_ndof()) = moment_of_inertia;

    // r_max = this->vertices.rowwise().norm().maxCoeff();

    // average_edge_length = 0;
    // for (long i = 0; i < edges.rows(); i++) {
    //     average_edge_length +=
    //         (this->vertices.row(edges(i, 0)) - this->vertices.row(edges(i,
    //         1)))
    //             .norm();
    // }
    // if (edges.rows() > 0) {
    //     average_edge_length /= edges.rows();
    // }
    // assert(std::isfinite(average_edge_length));

    // init_bvh();
}

// Eigen::MatrixXd
// RigidBody::transform_vertices(const Eigen::MatrixXd& rest_positions) const
// {
//     // Compute: R(θ) x̄ + p
//     // transpose because x is row-ordered
//     const int dim = rest_positions.cols();

//     assert(
//         (dim == 2 && rotation_vector.size() == 1)
//         || (dim == 3 && rotation_vector.size() == 3));

//     // Convert the rotation vector to a rotation matrix
//     MatrixMax3d R(dim, dim);
//     if (dim == 2) {
//         const double theta = rotation_vector(0);
//         R << sin(theta), -cos(theta), cos(theta), sin(theta);
//     } else {
//         assert(dim == 3);
//         const double sinc_angle = sinc_norm_x<double>(rotation_vector);
//         const double sinc_half_angle =
//             sinc_norm_x<double>((rotation_vector / 2).eval());
//         const Eigen::Matrix3d K = cross_product_matrix(rotation_vector);
//         const Eigen::Matrix3d K2 = K * K;
//         R = sinc_angle * K + 0.5 * sinc_half_angle * sinc_half_angle * K2;
//         R.diagonal().array() += 1.0;
//     }

//     return (rest_positions * R.transpose()).rowwise() + position.transpose();
// }

} // namespace ipc::rigid