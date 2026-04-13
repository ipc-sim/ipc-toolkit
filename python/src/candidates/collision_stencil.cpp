#include <common.hpp>

#include <ipc/candidates/collision_stencil.hpp>

using namespace ipc;

void define_collision_stencil(py::module_& m)
{
    py::class_<CollisionStencil>(m, "CollisionStencil")
        .def(
            "num_vertices", &CollisionStencil::num_vertices,
            "Get the number of vertices in the collision stencil.")
        .def(
            "dim", &CollisionStencil::dim,
            R"ipc_Qu8mg5v7(
            Get the dimension of the collision stencil.

            Parameters:
                ndof: Number of degrees of freedom in the stencil.

            Returns:
                The dimension of the collision stencil.
            )ipc_Qu8mg5v7",
            "ndof"_a)
        .def(
            "vertex_ids", &CollisionStencil::vertex_ids,
            R"ipc_Qu8mg5v7(
            Get the vertex IDs of the collision stencil.

            Parameters:
                edges: Collision mesh edges
                faces: Collision mesh faces

            Returns:
                The vertex IDs of the collision stencil. Size is always 4, but elements i > num_vertices() are -1.
            )ipc_Qu8mg5v7",
            "edges"_a, "faces"_a)
        .def(
            "vertices", &CollisionStencil::vertices,
            R"ipc_Qu8mg5v7(
            Get the vertex attributes of the collision stencil.

            T Type of the attributes

            Parameters:
                vertices: Vertex attributes
                edges: Collision mesh edges
                faces: Collision mesh faces

            Returns:
                The vertex positions of the collision stencil. Size is always 4, but elements i > num_vertices() are NaN.
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a)
        .def(
            "dof", &CollisionStencil::dof,
            R"ipc_Qu8mg5v7(
            Select this stencil's DOF from the full matrix of DOF.

            T Type of the DOF

            Parameters:
                X: Full matrix of DOF (rowwise).
                edges: Collision mesh edges
                faces: Collision mesh faces

            Returns:
                This stencil's DOF.
            )ipc_Qu8mg5v7",
            "X"_a, "edges"_a, "faces"_a)
        .def(
            "compute_distance",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>>(
                &CollisionStencil::compute_distance, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance of the stencil.

            Parameters:
                vertices: Collision mesh vertices.
                edges: Collision mesh edges.
                faces: Collision mesh faces.

            Returns:
                Distance of the stencil.
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a)
        .def(
            "compute_distance_gradient",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>>(
                &CollisionStencil::compute_distance_gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.

            Parameters:
                vertices: Collision mesh vertices.
                edges: Collision mesh edges.
                faces: Collision mesh faces.

            Returns:
                Distance gradient of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a)
        .def(
            "compute_distance_hessian",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>>(
                &CollisionStencil::compute_distance_hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.

            Parameters:
                vertices: Collision mesh vertices.
                edges: Collision mesh edges.
                faces: Collision mesh faces.

            Returns:
                Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a)
        .def(
            "compute_coefficients",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>>(
                &CollisionStencil::compute_coefficients, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.

            Parameters:
                vertices: Collision mesh vertices.
                edges: Collision mesh edges.
                faces: Collision mesh faces.

            Returns:
                Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a)
        .def(
            "compute_distance",
            py::overload_cast<Eigen::ConstRef<VectorMax12d>>(
                &CollisionStencil::compute_distance, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance of the stencil.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance of the stencil.
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def(
            "compute_distance_gradient",
            py::overload_cast<Eigen::ConstRef<VectorMax12d>>(
                &CollisionStencil::compute_distance_gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance gradient of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def(
            "compute_distance_hessian",
            py::overload_cast<Eigen::ConstRef<VectorMax12d>>(
                &CollisionStencil::compute_distance_hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def(
            "compute_coefficients",
            py::overload_cast<Eigen::ConstRef<VectorMax12d>>(
                &CollisionStencil::compute_coefficients, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance coefficients of the stencil w.r.t. the stencil's vertex positions.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Distance of the stencil.
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def(
            "ccd",
            [](const CollisionStencil& self,
               Eigen::ConstRef<VectorMax12d> vertices_t0,
               Eigen::ConstRef<VectorMax12d> vertices_t1,
               const double min_distance, const double tmax,
               const NarrowPhaseCCD& narrow_phase_ccd) {
                double toi;
                bool r = self.ccd(
                    vertices_t0, vertices_t1, toi, min_distance, tmax,
                    narrow_phase_ccd);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
                Perform narrow-phase CCD on the candidate.

                Parameters:
                    vertices_t0: Stencil vertices at the start of the time step.
                    vertices_t1: Stencil vertices at the end of the time step.
                    min_distance: Minimum separation distance between primitives.
                    tmax: Maximum time (normalized) to look for collisions. Should be in [0, 1].
                    narrow_phase_ccd: The narrow phase CCD algorithm to use.

                Returns:
                    Tuple of:
                    If the candidate had a collision over the time interval.
                    Computed time of impact (normalized).
                )ipc_Qu8mg5v7",
            "vertices_t0"_a, "vertices_t1"_a, "min_distance"_a = 0.0,
            "tmax"_a = 1.0, "narrow_phase_ccd"_a = DEFAULT_NARROW_PHASE_CCD)
        .def(
            "print_ccd_query",
            [](const CollisionStencil& self,
               Eigen::ConstRef<VectorMax12d> vertices_t0,
               Eigen::ConstRef<VectorMax12d> vertices_t1) -> void {
                self.write_ccd_query(std::cout, vertices_t0, vertices_t1);
            },
            R"ipc_Qu8mg5v7(
                Print the CCD query to cout.

                Parameters:
                                    vertices_t0: Stencil vertices at the start of the time step.
                    vertices_t1: Stencil vertices at the end of the time step.
                )ipc_Qu8mg5v7",
            "vertices_t0"_a, "vertices_t1"_a)
        .def(
            "compute_distance_vector",
            py::overload_cast<Eigen::ConstRef<VectorMax12d>>(
                &CollisionStencil::compute_distance_vector, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance vector of the stencil: t = sum(c_i * x_i).

            The distance vector is the vector between the closest points on the
            collision primitives. Its squared norm equals the squared distance.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                The distance vector (dim-dimensional, i.e., 2D or 3D).
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def(
            "compute_distance_vector_with_coefficients",
            [](const CollisionStencil& self,
               Eigen::ConstRef<VectorMax12d> positions) {
                VectorMax4d coeffs;
                VectorMax3d dv =
                    self.compute_distance_vector(positions, coeffs);
                return std::make_tuple(dv, coeffs);
            },
            R"ipc_Qu8mg5v7(
            Compute the distance vector and the coefficients together.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                Tuple of:
                The distance vector (dim-dimensional).
                The computed coefficients c_i.
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def(
            "compute_distance_vector",
            py::overload_cast<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>>(
                &CollisionStencil::compute_distance_vector, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the distance vector of the stencil.

            Parameters:
                vertices: Collision mesh vertices.
                edges: Collision mesh edges.
                faces: Collision mesh faces.

            Returns:
                The distance vector (dim-dimensional).
            )ipc_Qu8mg5v7",
            "vertices"_a, "edges"_a, "faces"_a)
        .def(
            "compute_distance_vector_jacobian",
            &CollisionStencil::compute_distance_vector_jacobian,
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the distance vector w.r.t. positions.

            J = [c_0 I, c_1 I, ..., c_n I]^T where I is the dim x dim identity.

            Note:
                positions can be computed as stencil.dof(vertices, edges, faces)

            Parameters:
                positions: Stencil's vertex positions.

            Returns:
                The Jacobian dt/dx as a matrix of shape (ndof, dim).
            )ipc_Qu8mg5v7",
            "positions"_a)
        .def_static(
            "diag_distance_vector_outer",
            &CollisionStencil::diag_distance_vector_outer,
            R"ipc_Qu8mg5v7(
            Compute diag((dt/dx)(dt/dx)^T) efficiently (Eq. 11).

            Result is [c_0^2, c_0^2, c_0^2, c_1^2, ...] (each c_i^2 repeated
            dim times).

            Parameters:
                coeffs: The coefficients c_i (from compute_coefficients).
                dim: The spatial dimension (2 or 3).

            Returns:
                The diagonal of (dt/dx)(dt/dx)^T as a vector of size ndof.
            )ipc_Qu8mg5v7",
            "coeffs"_a, "dim"_a)
        .def_static(
            "diag_distance_vector_t_outer",
            &CollisionStencil::diag_distance_vector_t_outer,
            R"ipc_Qu8mg5v7(
            Compute diag((dt/dx * t)(dt/dx * t)^T) efficiently (Eq. 12).

            Result is element-wise square of [c_0*t^T, c_1*t^T, ..., c_n*t^T].

            Parameters:
                coeffs: The coefficients c_i.
                distance_vector: The distance vector t.

            Returns:
                The diagonal of (dt/dx*t)(dt/dx*t)^T as a vector of size ndof.
            )ipc_Qu8mg5v7",
            "coeffs"_a, "distance_vector"_a)
        .def_static(
            "contract_distance_vector_jacobian",
            &CollisionStencil::contract_distance_vector_jacobian,
            R"ipc_Qu8mg5v7(
            Compute p^T (dt/dx) efficiently as sum(c_i * p_i) (Eqs. 13-14).

            Given p = [p_0, p_1, ..., p_n]^T where p_i are dim-dimensional,
            this computes p^T (dt/dx) = sum(c_i * p_i) which is a
            dim-dimensional vector.

            Parameters:
                coeffs: The coefficients c_i.
                p: A vector of size ndof (the direction for the quadratic form).
                dim: The spatial dimension (2 or 3).

            Returns:
                p^T (dt/dx) as a dim-dimensional vector.
            )ipc_Qu8mg5v7",
            "coeffs"_a, "p"_a, "dim"_a);
}
