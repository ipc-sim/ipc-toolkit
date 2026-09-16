#pragma once

#include <common.hpp>

#include <ipc/candidates/collision_stencil.hpp>

#include <iostream>
#include <tuple>

namespace ipc {

/// @brief Bind the stencil operations onto a candidate class.
///
/// The candidate types are no longer CollisionStencils -- the broad phase
/// stores millions of them and a vptr would double their size -- so Python no
/// longer inherits these from that base and each candidate binds them itself.
/// Lambdas rather than py::overload_cast, because the candidates whose
/// distance type can be known a priori take an extra defaulted argument and so
/// do not match the CollisionStencil signatures.
///
/// The full docstrings live on CollisionStencil; these are the same functions.
///
/// @tparam C The candidate type.
/// @param cls The pybind11 class to add the methods to.
template <class C, class PyClass> void define_stencil_methods(PyClass& cls)
{
    using CR12 = Eigen::ConstRef<VectorMax12d>;
    using CRXd = Eigen::ConstRef<Eigen::MatrixXd>;
    using CRXi = Eigen::ConstRef<Eigen::MatrixXi>;

    cls.def(
           "num_vertices", [](const C& self) { return self.num_vertices(); },
           "Get the number of vertices in the collision stencil.")
        .def(
            "dim", [](const C& self, const int ndof) { return self.dim(ndof); },
            "Get the dimension of the collision stencil.", "ndof"_a)
        .def(
            "vertex_ids",
            [](const C& self, CRXi edges, CRXi faces) {
                return self.vertex_ids(edges, faces);
            },
            "Get the vertex IDs of the collision stencil.", "edges"_a,
            "faces"_a)
        .def(
            "vertices",
            [](const C& self, CRXd V, CRXi edges, CRXi faces) {
                return self.vertices(V, edges, faces);
            },
            "Get the vertex attributes of the collision stencil.", "vertices"_a,
            "edges"_a, "faces"_a)
        .def(
            "dof",
            [](const C& self, CRXd X, CRXi edges, CRXi faces) {
                return self.dof(X, edges, faces);
            },
            "Select this stencil's DOF from the full matrix of DOF.", "X"_a,
            "edges"_a, "faces"_a)
        // -- distance ------------------------------------------------------
        .def(
            "compute_distance",
            [](const C& self, CR12 positions) {
                return self.compute_distance(positions);
            },
            "Compute the distance of the stencil.", "positions"_a)
        .def(
            "compute_distance",
            [](const C& self, CRXd V, CRXi edges, CRXi faces) {
                return self.compute_distance(V, edges, faces);
            },
            "Compute the distance of the stencil.", "vertices"_a, "edges"_a,
            "faces"_a)
        .def(
            "compute_distance_gradient",
            [](const C& self, CR12 positions) {
                return self.compute_distance_gradient(positions);
            },
            "Compute the distance gradient of the stencil.", "positions"_a)
        .def(
            "compute_distance_gradient",
            [](const C& self, CRXd V, CRXi edges, CRXi faces) {
                return self.compute_distance_gradient(V, edges, faces);
            },
            "Compute the distance gradient of the stencil.", "vertices"_a,
            "edges"_a, "faces"_a)
        .def(
            "compute_distance_hessian",
            [](const C& self, CR12 positions) {
                return self.compute_distance_hessian(positions);
            },
            "Compute the distance Hessian of the stencil.", "positions"_a)
        .def(
            "compute_distance_hessian",
            [](const C& self, CRXd V, CRXi edges, CRXi faces) {
                return self.compute_distance_hessian(V, edges, faces);
            },
            "Compute the distance Hessian of the stencil.", "vertices"_a,
            "edges"_a, "faces"_a)
        .def(
            "compute_coefficients",
            [](const C& self, CR12 positions) {
                return self.compute_coefficients(positions);
            },
            "Compute the coefficients of the stencil.", "positions"_a)
        .def(
            "compute_coefficients",
            [](const C& self, CRXd V, CRXi edges, CRXi faces) {
                return self.compute_coefficients(V, edges, faces);
            },
            "Compute the coefficients of the stencil.", "vertices"_a, "edges"_a,
            "faces"_a)
        // -- distance vector -----------------------------------------------
        .def(
            "compute_distance_vector",
            [](const C& self, CR12 positions) {
                return self.compute_distance_vector(positions);
            },
            "Compute the distance vector of the stencil: t = sum(c_i * x_i).",
            "positions"_a)
        .def(
            "compute_distance_vector",
            [](const C& self, CRXd V, CRXi edges, CRXi faces) {
                return self.compute_distance_vector(V, edges, faces);
            },
            "Compute the distance vector of the stencil.", "vertices"_a,
            "edges"_a, "faces"_a)
        .def(
            "compute_distance_vector_with_coefficients",
            [](const C& self, CR12 positions) {
                VectorMax4d coeffs;
                VectorMax3d dv =
                    self.compute_distance_vector(positions, coeffs);
                return std::make_tuple(dv, coeffs);
            },
            "Compute the distance vector and the coefficients together.",
            "positions"_a)
        .def(
            "compute_distance_vector_jacobian",
            [](const C& self, CR12 positions) {
                return self.compute_distance_vector_jacobian(positions);
            },
            "Compute the Jacobian of the distance vector w.r.t. positions.",
            "positions"_a)
        // -- CCD -----------------------------------------------------------
        .def(
            "ccd",
            [](const C& self, CR12 vertices_t0, CR12 vertices_t1,
               const double min_distance, const double tmax,
               const NarrowPhaseCCD& narrow_phase_ccd) {
                double toi;
                bool r = self.ccd(
                    vertices_t0, vertices_t1, toi, min_distance, tmax,
                    narrow_phase_ccd);
                return std::make_tuple(r, toi);
            },
            "Perform narrow-phase CCD on the candidate.", "vertices_t0"_a,
            "vertices_t1"_a, "min_distance"_a = 0.0, "tmax"_a = 1.0,
            "narrow_phase_ccd"_a = DEFAULT_NARROW_PHASE_CCD)
        .def(
            "print_ccd_query",
            [](const C& self, CR12 vertices_t0, CR12 vertices_t1) -> void {
                self.write_ccd_query(std::cout, vertices_t0, vertices_t1);
            },
            "Print the CCD query to cout.", "vertices_t0"_a, "vertices_t1"_a)
        // -- type-independent helpers --------------------------------------
        .def_static(
            "diag_distance_vector_outer", &C::diag_distance_vector_outer,
            "Compute diag((dt/dx)(dt/dx)^T) efficiently (Eq. 11).", "coeffs"_a,
            "dim"_a)
        .def_static(
            "diag_distance_vector_t_outer", &C::diag_distance_vector_t_outer,
            "Compute diag((dt/dx * t)(dt/dx * t)^T) efficiently (Eq. 12).",
            "coeffs"_a, "distance_vector"_a)
        .def_static(
            "contract_distance_vector_jacobian",
            &C::contract_distance_vector_jacobian,
            "Compute p^T (dt/dx) efficiently as sum(c_i * p_i) (Eqs. 13-14).",
            "coeffs"_a, "p"_a, "dim"_a);
}

} // namespace ipc
