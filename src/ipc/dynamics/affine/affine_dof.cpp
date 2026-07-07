#include "affine_dof.hpp"

namespace ipc::affine {

rigid::AffinePose
dof_to_pose(Eigen::ConstRef<Eigen::VectorXd> x, const size_t i, const int dim)
{
    const int ndof = dim + dim * dim;
    assert((i + 1) * ndof <= size_t(x.size()));

    rigid::AffinePose pose;
    pose.position = x.segment(i * ndof, dim);
    pose.rotation = x.segment(i * ndof + dim, dim * dim).reshaped(dim, dim);
    return pose;
}

std::vector<rigid::AffinePose>
dof_to_poses(Eigen::ConstRef<Eigen::VectorXd> x, const int dim)
{
    const int ndof = dim + dim * dim;
    assert(x.size() % ndof == 0);

    std::vector<rigid::AffinePose> poses(x.size() / ndof);
    for (size_t i = 0; i < poses.size(); i++) {
        poses[i] = dof_to_pose(x, i, dim);
    }
    return poses;
}

Eigen::MatrixXd
vertices(const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int dim = bodies.dim();

    Eigen::MatrixXd V(bodies.num_vertices(), dim);
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        V.middleRows(bodies.body_vertex_start(i), bodies.body_num_vertices(i)) =
            dof_to_pose(x, i, dim).transform_vertices(
                bodies.body_rest_positions(i));
    }
    return V;
}

Eigen::SparseMatrix<double> affine_jacobian(const rigid::RigidBodies& bodies)
{
    const int dim = bodies.dim();
    const int ndof = dim + dim * dim;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(bodies.num_vertices() * dim * (dim + 1));

    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        const Eigen::SparseMatrix<double> J_i =
            rigid::AffinePose::J(bodies.body_rest_positions(i));

        const index_t row_start = dim * bodies.body_vertex_start(i);
        const index_t col_start = ndof * i;
        for (int k = 0; k < J_i.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(J_i, k); it;
                 ++it) {
                triplets.emplace_back(
                    row_start + it.row(), col_start + it.col(), it.value());
            }
        }
    }

    Eigen::SparseMatrix<double> J(
        dim * bodies.num_vertices(), ndof * bodies.num_bodies());
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

} // namespace ipc::affine
