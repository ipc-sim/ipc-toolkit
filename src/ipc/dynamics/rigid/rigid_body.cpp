#include "rigid_body.hpp"

#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/dynamics/rigid/mass.hpp>
#include <ipc/math/sinc.hpp>

namespace ipc::rigid {

namespace {
    void center_vertices(
        Eigen::Ref<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        Pose& pose)
    {
        // compute the center of mass several times to get more accurate
        for (int i = 0; i < 10; i++) {
            double mass;
            VectorMax3d com;
            MatrixMax3d inertia;
            compute_mass_properties(
                vertices,
                (vertices.cols() == 2 || faces.size() == 0) ? edges : faces,
                1.0, // density (1.0 because we only want the center of mass)
                mass, com, inertia);
            vertices.rowwise() -= com.transpose();
            pose.position += com;
            if (com.squaredNorm() < 1e-8) {
                break;
            }
        }
    }

    /// Convert the classical principal moments of inertia to the diagonal
    /// second-moment matrix J = ∫ρ x̄x̄ᵀ (3D only; in 2D the second moment is
    /// computed directly).
    inline Eigen::DiagonalMatrix<double, Eigen::Dynamic, 3>
    compute_J(Eigen::ConstRef<VectorMax3d> I) // NOLINT
    {
        assert(I.size() == 3);
        Eigen::DiagonalMatrix<double, Eigen::Dynamic, 3> J(3); // NOLINT
        J.diagonal() << 0.5 * (-I.x() + I.y() + I.z()),
            0.5 * (I.x() - I.y() + I.z()), 0.5 * (I.x() + I.y() - I.z());
        return J;
    }
} // namespace

RigidBody::RigidBody(
    Eigen::Ref<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double density,
    Pose& initial_pose,
    const VectorMax6b& is_dof_fixed)
{
    assert(vertices.size() > 0);
    assert(edges.size() == 0 || edges.cols() == 2);
    assert(faces.size() == 0 || faces.cols() == 3);

    const int dim = vertices.cols();
    assert(dim == 2 || dim == 3);

    const int ndof = dim == 2 ? 3 : 6;
    assert(is_dof_fixed.size() == 0 || is_dof_fixed.size() == ndof);
    m_is_dof_fixed = is_dof_fixed.size() == ndof
        ? is_dof_fixed
        : VectorMax6b(VectorMax6b::Zero(ndof));

    // 1. Center the vertices, so the mass properties are computed correctly
    // TODO: This should not be necessary. Determine why the mass properties
    // are not computed correctly without centering the vertices.
    Pose centering_pose = Pose::Identity(dim);
    center_vertices(vertices, edges, faces, centering_pose);

    // 2. Compute the mass properties
    VectorMax3d center_of_mass;
    MatrixMax3d inertia_tensor;
    compute_mass_properties(
        vertices, (dim == 2 || faces.size() == 0) ? edges : faces, density,
        m_mass, center_of_mass, inertia_tensor);
    m_volume = m_mass / density;

    // 3. Convert the inertia tensor to the principal axes moments of inertia
    const int num_rot_dof_fixed =
        dim == 3 ? int(m_is_dof_fixed.tail(3).count()) : 0;
    if (dim == 3 && num_rot_dof_fixed == 2) {
        // Rotation is confined to a single world axis: keep the world-frame
        // moments of inertia and skip the principal-axes rotation R₀ so the
        // free rotation-vector component stays aligned with its world axis.
        m_moment_of_inertia = inertia_tensor.diagonal();
        m_R0 = Eigen::Matrix3d::Identity();
    } else if (dim == 3) {
        if (num_rot_dof_fixed == 1) {
            logger().warn(
                "Rigid body dynamics with two free rotational DOF has not "
                "been tested thoroughly.");
        }
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

        // The principal moments of inertia are the eigenvalues of the inertia
        // tensor.
        m_moment_of_inertia = solver.eigenvalues();
        if ((m_moment_of_inertia.array() < 0).any()) {
            logger().warn(
                "Negative moments of inertia ({}), inverting.",
                m_moment_of_inertia);
            // This typically only happens with negative ε inertias
            m_moment_of_inertia = m_moment_of_inertia.array().abs();
        }

        // The rotation from the principal inertial frame to the input world
        // frame.
        m_R0 = solver.eigenvectors();

        // Ensure that we have an orientation preserving transform
        if (m_R0.determinant() < 0.0) {
            m_R0.col(0) *= -1.0;
        }
        assert(m_R0.isUnitary(1e-9));

        // Remove the initial rotation from the rest vertices
        vertices = vertices * m_R0;

        // Store the initial rotation in the pose (R = RᵢR₀)
        Eigen::AngleAxisd r = Eigen::AngleAxisd(Eigen::Matrix3d(m_R0));
        centering_pose.rotation = r.angle() * r.axis();

        // Compute the diagonal second moment in the principal frame from the
        // classical inertia tensor.
        m_J = compute_J(m_moment_of_inertia);
    } else {
        // 2D: the inertia is the full 2×2 second moment ∫ρ x̄x̄ᵀ (about the
        // COM). Mirror the 3D principal-frame handling: diagonalize, fold the
        // principal rotation into the pose, and rotate the rest vertices.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver;
        solver.compute(Eigen::Matrix2d(inertia_tensor));
        assert(solver.info() == Eigen::Success);

        Eigen::Matrix2d R0 = solver.eigenvectors(); // NOLINT
        if (R0.determinant() < 0.0) {
            R0.col(0) *= -1.0;
        }
        assert(R0.isUnitary(1e-9));
        m_R0 = R0;

        // Remove the initial rotation from the rest vertices
        vertices = vertices * m_R0;

        // Store the initial rotation angle in the pose
        centering_pose.rotation.resize(1);
        centering_pose.rotation(0) = std::atan2(R0(1, 0), R0(0, 0));

        // Diagonal second moment in the principal frame
        m_J = Eigen::DiagonalMatrix<double, Eigen::Dynamic, 3>(
            solver.eigenvalues().cwiseAbs());

        // Scalar moment about z: I_z = ∫ρ(x² + y²) = tr(∫ρ x̄x̄ᵀ)
        m_moment_of_inertia.resize(1);
        m_moment_of_inertia(0) = m_J.diagonal().sum();
    }

    // 4. Apply the centering pose to the initial pose
    initial_pose = initial_pose * centering_pose;

    m_external_force = Pose::Identity(dim);

    m_bounding_radius = vertices.rowwise().norm().maxCoeff();

    m_bvh = std::make_shared<LBVH>();
    m_bvh->build(vertices, edges, faces);
}

} // namespace ipc::rigid