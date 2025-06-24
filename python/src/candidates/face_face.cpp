#include <common.hpp>

#include <ipc/candidates/face_face.hpp>
#include <ipc/utils/logger.hpp>

using namespace ipc;

void define_face_face_candidate(py::module_& m)
{
    py::class_<FaceFaceCandidate>(m, "FaceFaceCandidate")
        .def(py::init<index_t, index_t>(), "face0_id"_a, "face1_id"_a)
        .def(
            py::init([](std::tuple<index_t, index_t> face_ids) {
                return std::make_unique<FaceFaceCandidate>(
                    std::get<0>(face_ids), std::get<1>(face_ids));
            }),
            "face_ids"_a)
        .def(
            "__str__",
            [](const FaceFaceCandidate& ff) {
                return fmt::format("[{:d}, {:d}]", ff.face0_id, ff.face1_id);
            })
        .def(
            "__repr__",
            [](const FaceFaceCandidate& ff) {
                return fmt::format(
                    "FaceFaceCandidate({:d}, {:d})", ff.face0_id, ff.face1_id);
            })
        .def("__eq__", &FaceFaceCandidate::operator==, "other"_a)
        .def("__ne__", &FaceFaceCandidate::operator!=, "other"_a)
        .def(
            "__lt__", &FaceFaceCandidate::operator<,
            "Compare FaceFaceCandidate for sorting.", "other"_a)
        .def_readwrite(
            "face0_id", &FaceFaceCandidate::face0_id, "ID of the first face.")
        .def_readwrite(
            "face1_id", &FaceFaceCandidate::face1_id, "ID of the second face.");

    py::implicitly_convertible<
        std::tuple<index_t, index_t>, FaceFaceCandidate>();
}
