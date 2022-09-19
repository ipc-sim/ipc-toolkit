#include "../common.hpp"

#include <ipc/friction/friction_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_constraint(py::module_& m)
{
    py::class_<FrictionConstraint> friction_constraint(m, "FrictionConstraint");

    py::enum_<FrictionConstraint::DiffWRT>(friction_constraint, "DiffWRT")
        .value("X", FrictionConstraint::DiffWRT::X)
        .value("Ut", FrictionConstraint::DiffWRT::Ut)
        .value("U", FrictionConstraint::DiffWRT::U)
        .export_values();

    friction_constraint
        .def("num_vertices", &FrictionConstraint::num_vertices, "")
        .def(
            "vertex_indices", &FrictionConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential_gradient",
            &FrictionConstraint::compute_potential_gradient, "", py::arg("U"),
            py::arg("E"), py::arg("F"), py::arg("epsv_times_h"))
        .def(
            "compute_potential_hessian",
            &FrictionConstraint::compute_potential_hessian, "", py::arg("U"),
            py::arg("E"), py::arg("F"), py::arg("epsv_times_h"),
            py::arg("project_hessian_to_psd"))
        .def(
            "compute_force",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double, const double, const double, const bool>(
                &FrictionConstraint::compute_force, py::const_),
            "", py::arg("X"), py::arg("U"), py::arg("E"), py::arg("F"),
            py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("dmin") = 0,
            py::arg("no_mu") = false)
        .def(
            "compute_force",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, const double, const double,
                const double, const double, const bool>(
                &FrictionConstraint::compute_force, py::const_),
            "", py::arg("X"), py::arg("Ut"), py::arg("U"), py::arg("E"),
            py::arg("F"), py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("dmin") = 0,
            py::arg("no_mu") = false)
        .def(
            "compute_force_jacobian",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double, const double, FrictionConstraint::DiffWRT,
                const double>(
                &FrictionConstraint::compute_force_jacobian, py::const_),
            "", py::arg("X"), py::arg("U"), py::arg("E"), py::arg("F"),
            py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("wrt"), py::arg("dmin") = 0)
        .def(
            "compute_force_jacobian",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, const double, const double,
                const double, FrictionConstraint::DiffWRT, const double>(
                &FrictionConstraint::compute_force_jacobian, py::const_),
            "", py::arg("X"), py::arg("Ut"), py::arg("U"), py::arg("E"),
            py::arg("F"), py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("wrt"), py::arg("dmin") = 0)
        .def_readwrite(
            "closest_point", &FrictionConstraint::closest_point,
            "Barycentric coordinates of the closest point(s)")
        .def_readwrite(
            "tangent_basis", &FrictionConstraint::tangent_basis,
            "Tangent basis of the contact (max size 3×2)")
        .def_readwrite(
            "normal_force_magnitude",
            &FrictionConstraint::normal_force_magnitude,
            "Contact force magnitude")
        .def_readwrite("mu", &FrictionConstraint::mu, "Coefficient of friction")
        .def_readwrite(
            "weight", &FrictionConstraint::weight,
            "Weight in the final sum of potentials");

    py::class_<
        VertexVertexFrictionConstraint, VertexVertexCandidate,
        FrictionConstraint>(m, "VertexVertexFrictionConstraint")
        .def(
            py::init<long, long>(), "", py::arg("vertex0_index"),
            py::arg("vertex1_index"))
        .def(py::init<const VertexVertexCandidate&>(), "", py::arg("candidate"))
        .def(
            py::init<const VertexVertexConstraint&>(), "",
            py::arg("constraint"))
        .def(
            py::init<
                const VertexVertexConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            "", py::arg("constraint"), py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("dhat"), py::arg("barrier_stiffness"))
        .def("num_vertices", &VertexVertexFrictionConstraint::num_vertices, "")
        .def(
            "vertex_indices", &VertexVertexFrictionConstraint::vertex_indices,
            "", py::arg("E"), py::arg("F"))
        .def(
            "compute_potential",
            &VertexVertexFrictionConstraint::compute_potential<double>, "",
            py::arg("U"), py::arg("E"), py::arg("F"), py::arg("epsv_times_h"));

    py::class_<
        EdgeVertexFrictionConstraint, EdgeVertexCandidate, FrictionConstraint>(
        m, "EdgeVertexFrictionConstraint")
        .def(
            py::init<long, long>(), "", py::arg("edge_index"),
            py::arg("vertex_index"))
        .def(py::init<const EdgeVertexCandidate&>(), "", py::arg("constraint"))
        .def(py::init<const EdgeVertexConstraint&>(), "", py::arg("constraint"))
        .def(
            py::init<
                const EdgeVertexConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            "", py::arg("constraint"), py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("dhat"), py::arg("barrier_stiffness"))
        .def("num_vertices", &EdgeVertexFrictionConstraint::num_vertices, "")
        .def(
            "vertex_indices", &EdgeVertexFrictionConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential",
            &EdgeVertexFrictionConstraint::compute_potential<double>, "",
            py::arg("U"), py::arg("E"), py::arg("F"), py::arg("epsv_times_h"));

    py::class_<
        EdgeEdgeFrictionConstraint, EdgeEdgeCandidate, FrictionConstraint>(
        m, "EdgeEdgeFrictionConstraint")
        .def(
            py::init<long, long>(), "", py::arg("edge0_index"),
            py::arg("edge1_index"))
        .def(py::init<const EdgeEdgeCandidate&>(), "", py::arg("constraint"))
        .def(py::init<const EdgeEdgeConstraint&>(), "", py::arg("constraint"))
        .def(
            py::init<
                const EdgeEdgeConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            "", py::arg("constraint"), py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("dhat"), py::arg("barrier_stiffness"))
        .def("num_vertices", &EdgeEdgeFrictionConstraint::num_vertices, "")
        .def(
            "vertex_indices", &EdgeEdgeFrictionConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential",
            &EdgeEdgeFrictionConstraint::compute_potential<double>, "",
            py::arg("U"), py::arg("E"), py::arg("F"), py::arg("epsv_times_h"));

    py::class_<
        FaceVertexFrictionConstraint, FaceVertexCandidate, FrictionConstraint>(
        m, "FaceVertexFrictionConstraint")
        .def(
            py::init<long, long>(), "", py::arg("face_index"),
            py::arg("vertex_index"))
        .def(py::init<const FaceVertexCandidate&>(), "", py::arg("constraint"))
        .def(py::init<const FaceVertexConstraint&>(), "", py::arg("constraint"))
        .def(
            py::init<
                const FaceVertexConstraint&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double>(),
            "", py::arg("constraint"), py::arg("V"), py::arg("E"), py::arg("F"),
            py::arg("dhat"), py::arg("barrier_stiffness"))
        .def("num_vertices", &FaceVertexFrictionConstraint::num_vertices, "")
        .def(
            "vertex_indices", &FaceVertexFrictionConstraint::vertex_indices, "",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential",
            &FaceVertexFrictionConstraint::compute_potential<double>, "",
            py::arg("U"), py::arg("E"), py::arg("F"), py::arg("epsv_times_h"));

    py::class_<FrictionConstraints>(m, "FrictionConstraints")
        .def("__len__", &FrictionConstraints::size, "")
        .def("empty", &FrictionConstraints::empty, "")
        .def("clear", &FrictionConstraints::clear, "")
        .def(
            "__getitem__",
            [](FrictionConstraints& self, size_t idx) -> FrictionConstraint* {
                return &self[idx];
            })
        .def_readwrite(
            "vv_constraints", &FrictionConstraints::vv_constraints, "")
        .def_readwrite(
            "ev_constraints", &FrictionConstraints::ev_constraints, "")
        .def_readwrite(
            "ee_constraints", &FrictionConstraints::ee_constraints, "")
        .def_readwrite(
            "fv_constraints", &FrictionConstraints::fv_constraints, "");
}
