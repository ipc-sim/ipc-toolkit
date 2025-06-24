#include <common.hpp>

#include <ipc/candidates/edge_face.hpp>
#include <ipc/utils/logger.hpp>

using namespace ipc;

void define_edge_face_candidate(py::module_& m)
{
    py::class_<EdgeFaceCandidate>(m, "EdgeFaceCandidate")
        .def(py::init<index_t, index_t>(), "edge_id"_a, "face_id"_a)
        .def(
            py::init([](std::tuple<index_t, index_t> edge_and_face_id) {
                return std::make_unique<EdgeFaceCandidate>(
                    std::get<0>(edge_and_face_id),
                    std::get<1>(edge_and_face_id));
            }),
            "edge_and_face_id"_a)
        .def(
            "__str__",
            [](const EdgeFaceCandidate& ev) {
                return fmt::format("[{:d}, {:d}]", ev.edge_id, ev.face_id);
            })
        .def(
            "__repr__",
            [](const EdgeFaceCandidate& ev) {
                return fmt::format(
                    "EdgeFaceCandidate({:d}, {:d})", ev.edge_id, ev.face_id);
            })
        .def("__eq__", &EdgeFaceCandidate::operator==, "other"_a)
        .def("__ne__", &EdgeFaceCandidate::operator!=, "other"_a)
        .def(
            "__lt__", &EdgeFaceCandidate::operator<,
            "Compare EdgeFaceCandidate for sorting.", "other"_a)
        .def_readwrite("edge_id", &EdgeFaceCandidate::edge_id, "ID of the edge")
        .def_readwrite(
            "face_id", &EdgeFaceCandidate::face_id, "ID of the face");

    py::implicitly_convertible<
        std::tuple<index_t, index_t>, EdgeFaceCandidate>();
}
