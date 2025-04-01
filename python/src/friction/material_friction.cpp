#include <common.hpp>
#include <ipc/friction/material_friction.hpp>

namespace py = pybind11;
using namespace ipc;

void define_material_friction(py::module_& m)
{
    py::class_<MaterialFriction>(m, "MaterialFriction")
        .def(py::init<>(), "Default constructor")
        .def(py::init<const Eigen::MatrixXd&>(), 
            "Constructor with a friction table",
            py::arg("friction_table"))
        .def("set_friction_table", &MaterialFriction::set_friction_table,
            "Set the friction coefficient table",
            py::arg("friction_table"))
        .def("get_friction_table", &MaterialFriction::get_friction_table,
            "Get the friction coefficient table")
        .def("get_coefficient", &MaterialFriction::get_coefficient,
            "Get the friction coefficient for a collision",
            py::arg("collision"), py::arg("default_mu"));
}
