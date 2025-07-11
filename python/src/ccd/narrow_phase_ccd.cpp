#include <common.hpp>

#include <ipc/ccd/narrow_phase_ccd.hpp>

using namespace ipc;

class PyNarrowPhaseCCD : public NarrowPhaseCCD {
public:
    using NarrowPhaseCCD::NarrowPhaseCCD; // Inherit constructors
    bool point_point_ccd(
        Eigen::ConstRef<VectorMax3d> p0_t0,
        Eigen::ConstRef<VectorMax3d> p1_t0,
        Eigen::ConstRef<VectorMax3d> p0_t1,
        Eigen::ConstRef<VectorMax3d> p1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "point_point_ccd");
        if (override) {
            const py::tuple obj =
                override(p0_t0, p1_t0, p0_t1, p1_t1, min_distance, tmax);
            toi = obj[1].cast<double>();
            return obj[0].cast<bool>();
        }
        throw std::runtime_error("pure virtual function called");
    }
    bool point_edge_ccd(
        Eigen::ConstRef<VectorMax3d> p_t0,
        Eigen::ConstRef<VectorMax3d> e0_t0,
        Eigen::ConstRef<VectorMax3d> e1_t0,
        Eigen::ConstRef<VectorMax3d> p_t1,
        Eigen::ConstRef<VectorMax3d> e0_t1,
        Eigen::ConstRef<VectorMax3d> e1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "point_edge_ccd");
        if (override) {
            const py::tuple obj = override(
                p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, min_distance, tmax);
            toi = obj[1].cast<double>();
            return obj[0].cast<bool>();
        }
        throw std::runtime_error("pure virtual function called");
    }
    bool point_triangle_ccd(
        Eigen::ConstRef<Eigen::Vector3d> p_t0,
        Eigen::ConstRef<Eigen::Vector3d> t0_t0,
        Eigen::ConstRef<Eigen::Vector3d> t1_t0,
        Eigen::ConstRef<Eigen::Vector3d> t2_t0,
        Eigen::ConstRef<Eigen::Vector3d> p_t1,
        Eigen::ConstRef<Eigen::Vector3d> t0_t1,
        Eigen::ConstRef<Eigen::Vector3d> t1_t1,
        Eigen::ConstRef<Eigen::Vector3d> t2_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "point_triangle_ccd");
        if (override) {
            const py::tuple obj = override(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1,
                min_distance, tmax);
            toi = obj[1].cast<double>();
            return obj[0].cast<bool>();
        }
        throw std::runtime_error("pure virtual function called");
    }
    bool edge_edge_ccd(
        Eigen::ConstRef<Eigen::Vector3d> ea0_t0,
        Eigen::ConstRef<Eigen::Vector3d> ea1_t0,
        Eigen::ConstRef<Eigen::Vector3d> eb0_t0,
        Eigen::ConstRef<Eigen::Vector3d> eb1_t0,
        Eigen::ConstRef<Eigen::Vector3d> ea0_t1,
        Eigen::ConstRef<Eigen::Vector3d> ea1_t1,
        Eigen::ConstRef<Eigen::Vector3d> eb0_t1,
        Eigen::ConstRef<Eigen::Vector3d> eb1_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "edge_edge_ccd");
        if (override) {
            const py::tuple obj = override(
                ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
                min_distance, tmax);
            toi = obj[1].cast<double>();
            return obj[0].cast<bool>();
        }
        throw std::runtime_error("pure virtual function called");
    }
};

void define_narrow_phase_ccd(py::module_& m)
{
    py::class_<NarrowPhaseCCD, PyNarrowPhaseCCD>(m, "NarrowPhaseCCD")
        .def(py::init<>())
        .def(
            "point_point_ccd",
            [](const NarrowPhaseCCD& self, Eigen::ConstRef<VectorMax3d> p0_t0,
               Eigen::ConstRef<VectorMax3d> p1_t0,
               Eigen::ConstRef<VectorMax3d> p0_t1,
               Eigen::ConstRef<VectorMax3d> p1_t1,
               const double min_distance = 0.0, const double tmax = 1.0) {
                double toi;
                bool r = self.point_point_ccd(
                    p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            "p0_t0"_a, "p1_t0"_a, "p0_t1"_a, "p1_t1"_a, "min_distance"_a = 0.0,
            "tmax"_a = 1.0)
        .def(
            "point_edge_ccd",
            [](const NarrowPhaseCCD& self, Eigen::ConstRef<VectorMax3d> p_t0,
               Eigen::ConstRef<VectorMax3d> e0_t0,
               Eigen::ConstRef<VectorMax3d> e1_t0,
               Eigen::ConstRef<VectorMax3d> p_t1,
               Eigen::ConstRef<VectorMax3d> e0_t1,
               Eigen::ConstRef<VectorMax3d> e1_t1,
               const double min_distance = 0.0, const double tmax = 1.0) {
                double toi;
                bool r = self.point_edge_ccd(
                    p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi, min_distance,
                    tmax);
                return std::make_tuple(r, toi);
            },
            "p_t0"_a, "e0_t0"_a, "e1_t0"_a, "p_t1"_a, "e0_t1"_a, "e1_t1"_a,
            "min_distance"_a = 0.0, "tmax"_a = 1.0)
        .def(
            "point_triangle_ccd",
            [](const NarrowPhaseCCD& self,
               Eigen::ConstRef<Eigen::Vector3d> p_t0,
               Eigen::ConstRef<Eigen::Vector3d> t0_t0,
               Eigen::ConstRef<Eigen::Vector3d> t1_t0,
               Eigen::ConstRef<Eigen::Vector3d> t2_t0,
               Eigen::ConstRef<Eigen::Vector3d> p_t1,
               Eigen::ConstRef<Eigen::Vector3d> t0_t1,
               Eigen::ConstRef<Eigen::Vector3d> t1_t1,
               Eigen::ConstRef<Eigen::Vector3d> t2_t1,
               const double min_distance = 0.0, const double tmax = 1.0) {
                double toi;
                bool r = self.point_triangle_ccd(
                    p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi,
                    min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            "p_t0"_a, "t0_t0"_a, "t1_t0"_a, "t2_t0"_a, "p_t1"_a, "t0_t1"_a,
            "t1_t1"_a, "t2_t1"_a, "min_distance"_a = 0.0, "tmax"_a = 1.0)
        .def(
            "edge_edge_ccd",
            [](const NarrowPhaseCCD& self,
               Eigen::ConstRef<Eigen::Vector3d> ea0_t0,
               Eigen::ConstRef<Eigen::Vector3d> ea1_t0,
               Eigen::ConstRef<Eigen::Vector3d> eb0_t0,
               Eigen::ConstRef<Eigen::Vector3d> eb1_t0,
               Eigen::ConstRef<Eigen::Vector3d> ea0_t1,
               Eigen::ConstRef<Eigen::Vector3d> ea1_t1,
               Eigen::ConstRef<Eigen::Vector3d> eb0_t1,
               Eigen::ConstRef<Eigen::Vector3d> eb1_t1,
               const double min_distance = 0.0, const double tmax = 1.0) {
                double toi;
                bool r = self.edge_edge_ccd(
                    ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
                    eb1_t1, toi, min_distance, tmax);
                return std::make_tuple(r, toi);
            },
            "ea0_t0"_a, "ea1_t0"_a, "eb0_t0"_a, "eb1_t0"_a, "ea0_t1"_a,
            "ea1_t1"_a, "eb0_t1"_a, "eb1_t1"_a, "min_distance"_a = 0.0,
            "tmax"_a = 1.0);
}
