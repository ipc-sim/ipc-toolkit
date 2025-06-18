#include "mass.hpp"

#include <ipc/utils/logger.hpp>

#include <igl/massmatrix.h>

namespace ipc::rigid {

namespace {

    void compute_mass_properties_2D(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const double density,
        double& total_mass,
        VectorMax3d& center,
        MatrixMax3d& inertia)
    {
        Eigen::SparseMatrix<double> mass_matrix;
        construct_mass_matrix(vertices, edges, mass_matrix);
        total_mass = mass_matrix.sum();

        if (total_mass == 0) {
            center.setZero(vertices.cols());
        } else {
            center = (mass_matrix * vertices).colwise().sum() / total_mass;
        }

        // ∑ mᵢ rᵢ ⋅ rᵢ
        inertia.resize(1, 1);
        inertia(0) = (mass_matrix * vertices.rowwise().squaredNorm()).sum();

        // Total mass above is the paremeter length of the edges, so we need to
        // multiply by the density to get the total mass in the correct units.
        total_mass *= density;
        // Same for the inertia.
        inertia *= density;
    }

    // Based on ChTriangleMeshConnected.cpp::ComputeMassProperties from Chrono:
    // Copyright (c) 2016, Project Chrono Development Team
    // All rights reserved.
    // https://github.com/projectchrono/chrono/blob/main/LICENSE
    //
    // This requires the mesh to be closed, watertight, with proper triangle
    // orientation.
    bool compute_mass_properties_3D(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& faces,
        const double density,
        double& total_mass,
        VectorMax3d& center,
        MatrixMax3d& inertia)
    {
        assert(vertices.cols() == 3);
        assert(faces.rows() > 0 && faces.cols() == 3);

        // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
        std::array<double, 10> integral = { 0.0 };

        for (int i = 0; i < faces.rows(); i++) {
            // Get vertices of triangle i.
            const Eigen::Vector3d& v0 = vertices.row(faces(i, 0));
            const Eigen::Vector3d& v1 = vertices.row(faces(i, 1));
            const Eigen::Vector3d& v2 = vertices.row(faces(i, 2));

            // Get cross product of edges and normal vector.
            const Eigen::Vector3d& V1mV0 = v1 - v0;
            const Eigen::Vector3d& V2mV0 = v2 - v0;
            const Eigen::Vector3d& N = V1mV0.cross(V2mV0);

            // Compute integral terms.
            double tmp0, tmp1, tmp2;
            double f1x, f2x, f3x, g0x, g1x, g2x;
            tmp0 = v0.x() + v1.x();
            f1x = tmp0 + v2.x();
            tmp1 = v0.x() * v0.x();
            tmp2 = tmp1 + v1.x() * tmp0;
            f2x = tmp2 + v2.x() * f1x;
            f3x = v0.x() * tmp1 + v1.x() * tmp2 + v2.x() * f2x;
            g0x = f2x + v0.x() * (f1x + v0.x());
            g1x = f2x + v1.x() * (f1x + v1.x());
            g2x = f2x + v2.x() * (f1x + v2.x());

            double f1y, f2y, f3y, g0y, g1y, g2y;
            tmp0 = v0.y() + v1.y();
            f1y = tmp0 + v2.y();
            tmp1 = v0.y() * v0.y();
            tmp2 = tmp1 + v1.y() * tmp0;
            f2y = tmp2 + v2.y() * f1y;
            f3y = v0.y() * tmp1 + v1.y() * tmp2 + v2.y() * f2y;
            g0y = f2y + v0.y() * (f1y + v0.y());
            g1y = f2y + v1.y() * (f1y + v1.y());
            g2y = f2y + v2.y() * (f1y + v2.y());

            double f1z, f2z, f3z, g0z, g1z, g2z;
            tmp0 = v0.z() + v1.z();
            f1z = tmp0 + v2.z();
            tmp1 = v0.z() * v0.z();
            tmp2 = tmp1 + v1.z() * tmp0;
            f2z = tmp2 + v2.z() * f1z;
            f3z = v0.z() * tmp1 + v1.z() * tmp2 + v2.z() * f2z;
            g0z = f2z + v0.z() * (f1z + v0.z());
            g1z = f2z + v1.z() * (f1z + v1.z());
            g2z = f2z + v2.z() * (f1z + v2.z());

            // Update integrals.
            integral[0] += N.x() * f1x;
            integral[1] += N.x() * f2x;
            integral[2] += N.y() * f2y;
            integral[3] += N.z() * f2z;
            integral[4] += N.x() * f3x;
            integral[5] += N.y() * f3y;
            integral[6] += N.z() * f3z;
            integral[7] += N.x() * (v0.y() * g0x + v1.y() * g1x + v2.y() * g2x);
            integral[8] += N.y() * (v0.z() * g0y + v1.z() * g1y + v2.z() * g2y);
            integral[9] += N.z() * (v0.x() * g0z + v1.x() * g1z + v2.x() * g2z);
        }

        constexpr double inv_6 = 1.0 / 6.0;
        constexpr double inv_24 = 1.0 / 24.0;
        constexpr double inv_60 = 1.0 / 60.0;
        constexpr double inv_120 = 1.0 / 120.0;
        integral[0] *= inv_6;
        integral[1] *= inv_24;
        integral[2] *= inv_24;
        integral[3] *= inv_24;
        integral[4] *= inv_60;
        integral[5] *= inv_60;
        integral[6] *= inv_60;
        integral[7] *= inv_120;
        integral[8] *= inv_120;
        integral[9] *= inv_120;

        // total_mass
        total_mass = integral[0];
        if (total_mass <= 0 || !std::isfinite(total_mass)) {
            // 3D mass computation only works for closed meshes!
            return false;
        }
        assert(total_mass > 0);

        // center of mass
        center =
            Eigen::Vector3d(integral[1], integral[2], integral[3]) / total_mass;

        // inertia relative to world origin
        inertia.resize(3, 3);
        inertia(0, 0) = integral[5] + integral[6];
        inertia(0, 1) = -integral[7];
        inertia(0, 2) = -integral[9];
        inertia(1, 0) = inertia(0, 1);
        inertia(1, 1) = integral[4] + integral[6];
        inertia(1, 2) = -integral[8];
        inertia(2, 0) = inertia(0, 2);
        inertia(2, 1) = inertia(1, 2);
        inertia(2, 2) = integral[4] + integral[5];

        // inertia relative to center of mass
        inertia(0, 0) -=
            total_mass * (center.y() * center.y() + center.z() * center.z());
        inertia(0, 1) += total_mass * center.x() * center.y();
        inertia(0, 2) += total_mass * center.z() * center.x();
        inertia(1, 0) = inertia(0, 1);
        inertia(1, 1) -=
            total_mass * (center.z() * center.z() + center.x() * center.x());
        inertia(1, 2) += total_mass * center.y() * center.z();
        inertia(2, 0) = inertia(0, 2);
        inertia(2, 1) = inertia(1, 2);
        inertia(2, 2) -=
            total_mass * (center.x() * center.x() + center.y() * center.y());

        // Total mass above is the volume of the mesh, so we need to multiply by
        // the density to get the total mass in the correct units.
        total_mass *= density;
        // Same for the inertia.
        inertia *= density;

        return true;
    }

    void compute_mass_properties_point_cloud(
        const Eigen::MatrixXd& vertices,
        const double density,
        double& total_mass,
        VectorMax3d& center,
        MatrixMax3d& inertia)
    {
        assert(vertices.size() > 0);
        total_mass = vertices.rows(); // Point cloud has unit mass per point
        center = vertices.colwise().mean();
        inertia = Eigen::Matrix3d::Zero();
        // https://i.ytimg.com/vi/9jOrDufoO50/hqdefault.jpg
        for (long i = 0; i < vertices.rows(); i++) {
            const Eigen::Vector3d v = vertices.row(i);
            const Eigen::Vector3d v_sqr = v.array().pow(2);
            inertia(0, 0) += v_sqr.y() + v_sqr.z();
            inertia(1, 1) += v_sqr.x() + v_sqr.z();
            inertia(2, 2) += v_sqr.x() + v_sqr.y();
            inertia(0, 1) += -v.x() * v.y();
            inertia(0, 2) += -v.x() * v.z();
            inertia(1, 2) += -v.y() * v.z();
        }
        inertia(1, 0) = inertia(0, 1);
        inertia(2, 0) = inertia(0, 2);
        inertia(2, 1) = inertia(1, 2);

        // Total mass above is the number of points, so we need to multiply by
        // the density to get the total mass in the correct units.
        total_mass *= density;
        // Same for the inertia.
        inertia *= density;
    }

} // namespace

void compute_mass_properties(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& facets,
    const double density,
    double& total_mass,
    VectorMax3d& center,
    MatrixMax3d& inertia)
{
    assert(facets.cols() <= 3);
    if (vertices.cols() == 2) {
        compute_mass_properties_2D(
            vertices, facets, density, total_mass, center, inertia);
    }
    if (facets.size() == 0 || facets.cols() != 3) {
        compute_mass_properties_point_cloud(
            vertices, density, total_mass, center, inertia);
    } else if (
        facets.rows() < 4 // Not enough faces to form a closed mesh
        || !compute_mass_properties_3D(
            vertices, facets, density, total_mass, center, inertia)) {
        // If the mesh is not closed, we fall back to treating it as a point
        // cloud.
        compute_mass_properties_point_cloud(
            vertices, density, total_mass, center, inertia);
    }
}

// Construct the sparse mass matrix for the given mesh (V, F).
void construct_mass_matrix(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& facets,
    Eigen::SparseMatrix<double>& mass_matrix)
{
    if (vertices.cols() == 2 || facets.cols() == 2) {
        assert(facets.cols() == 2);
        Eigen::VectorXd vertex_masses = Eigen::VectorXd::Zero(vertices.rows());
        for (long i = 0; i < facets.rows(); i++) {
            const double edge_length =
                (vertices.row(facets(i, 1)) - vertices.row(facets(i, 0)))
                    .norm();
            // Add Voronoi areas to the vertex weight
            vertex_masses(facets(i, 0)) += edge_length / 2;
            vertex_masses(facets(i, 1)) += edge_length / 2;
        }
        mass_matrix.resize(vertices.rows(), vertices.rows());
        mass_matrix.diagonal() = vertex_masses;
    } else if (facets.cols() == 3) {
        assert(vertices.cols() == 3); // Only use triangles in 3D
        igl::massmatrix(
            vertices, facets, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI,
            mass_matrix);
    } else {
        // Probably a point cloud
        mass_matrix.resize(vertices.rows(), vertices.rows());
        mass_matrix.setIdentity();
    }
}

} // namespace ipc::rigid
