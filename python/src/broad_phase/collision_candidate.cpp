#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/utils/logger.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_candidate_classes(py::module_& m)
{
    py::class_<VertexVertexCandidate>(m, "VertexVertexCandidate")
        .def(py::init<long, long>())
        .def(
            "__str__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.vertex0_index, ev.vertex1_index);
            })
        .def(
            "__repr__",
            [](const VertexVertexCandidate& ev) {
                return fmt::format(
                    "VertexVertexCandidate({:d}, {:d})", ev.vertex0_index,
                    ev.vertex1_index);
            })
        .def("__eq__", &VertexVertexCandidate::operator==)
        .def("__ne__", &VertexVertexCandidate::operator!=)
        .def("__lt__", &VertexVertexCandidate::operator<)
        .def_readwrite("vertex0_index", &VertexVertexCandidate::vertex0_index)
        .def_readwrite("vertex1_index", &VertexVertexCandidate::vertex1_index);

    py::class_<EdgeVertexCandidate>(m, "EdgeVertexCandidate")
        .def(py::init<long, long>())
        .def(
            "__str__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge_index, ev.vertex_index);
            })
        .def(
            "__repr__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "EdgeVertexCandidate({:d}, {:d})", ev.edge_index,
                    ev.vertex_index);
            })
        .def("__eq__", &EdgeVertexCandidate::operator==)
        .def("__ne__", &EdgeVertexCandidate::operator!=)
        .def("__lt__", &EdgeVertexCandidate::operator<)
        .def_readwrite("edge_index", &EdgeVertexCandidate::edge_index)
        .def_readwrite("vertex_index", &EdgeVertexCandidate::vertex_index);

    py::class_<EdgeEdgeCandidate>(m, "EdgeEdgeCandidate")
        .def(py::init<long, long>())
        .def(
            "__str__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge0_index, ev.edge1_index);
            })
        .def(
            "__repr__",
            [](const EdgeEdgeCandidate& ev) {
                return fmt::format(
                    "EdgeEdgeCandidate({:d}, {:d})", ev.edge0_index,
                    ev.edge1_index);
            })
        .def("__eq__", &EdgeEdgeCandidate::operator==)
        .def("__ne__", &EdgeEdgeCandidate::operator!=)
        .def("__lt__", &EdgeEdgeCandidate::operator<)
        .def_readwrite("edge0_index", &EdgeEdgeCandidate::edge0_index)
        .def_readwrite("edge1_index", &EdgeEdgeCandidate::edge1_index);

    py::class_<EdgeFaceCandidate>(m, "EdgeFaceCandidate")
        .def(py::init<long, long>())
        .def(
            "__str__",
            [](const EdgeFaceCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge_index, ev.face_index);
            })
        .def(
            "__repr__",
            [](const EdgeFaceCandidate& ev) {
                return fmt::format(
                    "EdgeFaceCandidate({:d}, {:d})", ev.edge_index,
                    ev.face_index);
            })
        .def("__eq__", &EdgeFaceCandidate::operator==)
        .def("__ne__", &EdgeFaceCandidate::operator!=)
        .def("__lt__", &EdgeFaceCandidate::operator<)
        .def_readwrite("edge_index", &EdgeFaceCandidate::edge_index)
        .def_readwrite("face_index", &EdgeFaceCandidate::face_index);

    py::class_<FaceVertexCandidate>(m, "FaceVertexCandidate")
        .def(py::init<long, long>())
        .def(
            "__str__",
            [](const FaceVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.face_index, ev.vertex_index);
            })
        .def(
            "__repr__",
            [](const FaceVertexCandidate& ev) {
                return fmt::format(
                    "FaceVertexCandidate({:d}, {:d})", ev.face_index,
                    ev.vertex_index);
            })
        .def("__eq__", &FaceVertexCandidate::operator==)
        .def("__ne__", &FaceVertexCandidate::operator!=)
        .def("__lt__", &FaceVertexCandidate::operator<)
        .def_readwrite("face_index", &FaceVertexCandidate::face_index)
        .def_readwrite("vertex_index", &FaceVertexCandidate::vertex_index);

    py::class_<Candidates>(m, "Candidates")
        .def(py::init())
        .def("size", &Candidates::size)
        .def("empty", &Candidates::empty)
        .def("clear", &Candidates::clear)
        .def_readwrite("ev_candidates", &Candidates::ev_candidates)
        .def_readwrite("ee_candidates", &Candidates::ee_candidates)
        .def_readwrite("fv_candidates", &Candidates::fv_candidates);
}
