#include <common.hpp>

#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>

using namespace ipc;

void define_adaptive_stiffness(py::module_& m)
{
    m.def(
        "initial_barrier_stiffness",
        [](const double bbox_diagonal, const Barrier& barrier,
           const double dhat, const double average_mass,
           Eigen::ConstRef<Eigen::VectorXd> grad_energy,
           Eigen::ConstRef<Eigen::VectorXd> grad_barrier,
           const double min_barrier_stiffness_scale = 1e11,
           const double dmin = 0) {
            double max_barrier_stiffness;
            double r = initial_barrier_stiffness(
                bbox_diagonal, barrier, dhat, average_mass, grad_energy,
                grad_barrier, max_barrier_stiffness,
                min_barrier_stiffness_scale, dmin);
            return std::make_tuple(r, max_barrier_stiffness);
        },
        R"ipc_Qu8mg5v7(
        Compute an inital barrier stiffness using the barrier potential gradient.

        Parameters:
            bbox_diagonal: Length of the diagonal of the bounding box of the scene.
            barrier: Barrier function.
            dhat: Activation distance of the barrier.
            average_mass: Average mass of all bodies.
            grad_energy: Gradient of the elasticity energy function.
            grad_barrier: Gradient of the barrier potential.
            min_barrier_stiffness_scale: Scale used to premultiply the minimum barrier stiffness.
            dmin: Minimum distance between elements.

        Returns:
            Tuple of:
            The initial barrier stiffness.
            Maximum stiffness of the barrier.
        )ipc_Qu8mg5v7",
        "bbox_diagonal"_a, "barrier"_a, "dhat"_a, "average_mass"_a,
        "grad_energy"_a, "grad_barrier"_a,
        "min_barrier_stiffness_scale"_a = 1e11, "dmin"_a = 0);

    m.def(
        "update_barrier_stiffness", &update_barrier_stiffness,
        R"ipc_Qu8mg5v7(
        Update the barrier stiffness if the distance is decreasing and less than dhat_epsilon_scale * diag.

        Parameters:
            prev_min_distance: Previous minimum distance between elements.
            min_distance: Current minimum distance between elements.
            max_barrier_stiffness: Maximum stiffness of the barrier.
            barrier_stiffness: Current barrier stiffness.
            bbox_diagonal: Length of the diagonal of the bounding box of the scene.
            dhat_epsilon_scale: Update if distance is less than this fraction of the diagonal.
            dmin: Minimum distance between elements.

        Returns:
            The updated barrier stiffness.
        )ipc_Qu8mg5v7",
        "prev_min_distance"_a, "min_distance"_a, "max_barrier_stiffness"_a,
        "barrier_stiffness"_a, "bbox_diagonal"_a, "dhat_epsilon_scale"_a = 1e-9,
        "dmin"_a = 0);

    m.def(
        "semi_implicit_stiffness",
        static_cast<double (*)(
            const CollisionStencil&, Eigen::ConstRef<VectorMax12d>,
            Eigen::ConstRef<VectorMax4d>, Eigen::ConstRef<MatrixMax12d>,
            const double)>(&semi_implicit_stiffness),
        R"ipc_Qu8mg5v7(
        Compute the semi-implicit stiffness for a single collision.

        See [Ando 2024] for details.

        Parameters:
            stencil: Collision stencil.
            vertex_ids: Vertex indices of the collision.
            vertices: Vertex positions.
            mass: Vertex masses.
            local_hess: Local hessian of the elasticity energy function.
            dmin: Minimum distance between elements.

        Returns:
            The semi-implicit stiffness.
        )ipc_Qu8mg5v7",
        "stencil"_a, "vertices"_a, "mass"_a, "local_hess"_a, "dmin"_a);

    m.def(
        "semi_implicit_stiffness",
        static_cast<Eigen::VectorXd (*)(
            const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
            const NormalCollisions&, Eigen::ConstRef<Eigen::VectorXd>,
            const Eigen::SparseMatrix<double>&, const double)>(
            &semi_implicit_stiffness),
        R"ipc_Qu8mg5v7(
        Compute the semi-implicit stiffness's for all collisions.

        See [Ando 2024] for details.

        Parameters:
            mesh: Collision mesh.
            vertices: Vertex positions.
            collisions: Normal collisions.
            vertex_masses: Lumped vertex masses.
            hess: Hessian of the elasticity energy function.
            dmin: Minimum distance between elements.

        Returns:
            The semi-implicit stiffness's.
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "collisions"_a, "vertex_masses"_a, "hess"_a,
        "dmin"_a);

    m.def(
        "semi_implicit_stiffness",
        static_cast<Eigen::VectorXd (*)(
            const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
            const Candidates&, Eigen::ConstRef<Eigen::VectorXd>,
            const Eigen::SparseMatrix<double>&, const double)>(
            &semi_implicit_stiffness),
        R"ipc_Qu8mg5v7(
        Compute the semi-implicit stiffness's for all collisions.

        See [Ando 2024] for details.

        Parameters:
            mesh: Collision mesh.
            vertices: Vertex positions.
            collisions: Collisions candidates.
            vertex_masses: Lumped vertex masses.
            hess: Hessian of the elasticity energy function.
            dmin: Minimum distance between elements.

        Returns:
            The semi-implicit stiffness's.
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "collisions"_a, "vertex_masses"_a, "hess"_a,
        "dmin"_a);
}
