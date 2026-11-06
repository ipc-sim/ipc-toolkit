#include <common.hpp>

#include <ipc/ccd/narrow_phase_ccd.hpp>

#include <array>

using namespace ipc;

namespace {
/// @brief Whether a Python override accepts the `mode` keyword argument.
///
/// Overrides written before NarrowPhaseCCDMode existed take no `mode`
/// parameter, so we only pass it to overrides that declare it (or accept
/// **kwargs). Older overrides keep working; they just never learn that a
/// query comes from a piecewise driver.
bool accepts_mode(const py::function& override)
{
    try {
        const py::module_ inspect = py::module_::import("inspect");
        const py::object params =
            inspect.attr("signature")(override).attr("parameters");
        if (params.contains("mode")) {
            return true;
        }
        const py::object var_keyword =
            inspect.attr("Parameter").attr("VAR_KEYWORD");
        for (const py::handle param : params.attr("values")()) {
            const py::object kind = param.attr("kind");
            if (kind.equal(var_keyword)) {
                return true;
            }
        }
        return false;
    } catch (const py::error_already_set&) {
        return true; // No inspectable signature, so we pass the argument.
    }
}
} // namespace

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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "point_point_ccd");
        if (override) {
            const py::tuple obj = call_override(
                0, override, mode, p0_t0, p1_t0, p0_t1, p1_t1, min_distance,
                tmax);
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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "point_edge_ccd");
        if (override) {
            const py::tuple obj = call_override(
                1, override, mode, p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1,
                min_distance, tmax);
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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "point_triangle_ccd");
        if (override) {
            const py::tuple obj = call_override(
                2, override, mode, p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1,
                t1_t1, t2_t1, min_distance, tmax);
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
        const double tmax = 1.0,
        const NarrowPhaseCCDMode mode =
            NarrowPhaseCCDMode::STANDALONE) const override
    {
        py::gil_scoped_acquire gil; // Acquire GIL before calling Python code
        py::function override = py::get_override(this, "edge_edge_ccd");
        if (override) {
            const py::tuple obj = call_override(
                3, override, mode, ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1,
                ea1_t1, eb0_t1, eb1_t1, min_distance, tmax);
            toi = obj[1].cast<double>();
            return obj[0].cast<bool>();
        }
        throw std::runtime_error("pure virtual function called");
    }

private:
    /// @brief Call a Python override, passing `mode` only if it accepts it.
    /// @note Called with the GIL held, which also guards m_accepts_mode.
    template <typename... Args>
    py::object call_override(
        const int method,
        const py::function& override,
        const NarrowPhaseCCDMode mode,
        Args&&... args) const
    {
        int8_t& accepts = m_accepts_mode[method];
        if (accepts < 0) {
            accepts = accepts_mode(override) ? 1 : 0;
        }
        if (accepts) {
            return override(std::forward<Args>(args)..., "mode"_a = mode);
        }
        return override(std::forward<Args>(args)...);
    }

    /// @brief Per-method cache of accepts_mode() (-1 = unknown).
    mutable std::array<int8_t, 4> m_accepts_mode { { -1, -1, -1, -1 } };
};

void define_narrow_phase_ccd(py::module_& m)
{
    py::enum_<NarrowPhaseCCDMode>(
        m, "NarrowPhaseCCDMode", "How a narrow-phase CCD query is used.")
        .value(
            "STANDALONE", NarrowPhaseCCDMode::STANDALONE,
            "A standalone query: the method applies its own conservative "
            "strategy.")
        .value(
            "COMPONENT", NarrowPhaseCCDMode::COMPONENT,
            "The per-piece linear CCD inside a piecewise driver (e.g., "
            "NonlinearCCD): return the raw time of impact.")
        .value(
            "COMPONENT_NO_ZERO_TOI", NarrowPhaseCCDMode::COMPONENT_NO_ZERO_TOI,
            "Like COMPONENT, but a time of impact of exactly zero is "
            "unacceptable.");

    py::class_<
        NarrowPhaseCCD, PyNarrowPhaseCCD, std::shared_ptr<NarrowPhaseCCD>>(
        m, "NarrowPhaseCCD")
        .def(py::init<>())
        .def(
            "point_point_ccd",
            [](const NarrowPhaseCCD& self, Eigen::ConstRef<VectorMax3d> p0_t0,
               Eigen::ConstRef<VectorMax3d> p1_t0,
               Eigen::ConstRef<VectorMax3d> p0_t1,
               Eigen::ConstRef<VectorMax3d> p1_t1,
               const double min_distance = 0.0, const double tmax = 1.0,
               const NarrowPhaseCCDMode mode = NarrowPhaseCCDMode::STANDALONE) {
                double toi;
                bool r = self.point_point_ccd(
                    p0_t0, p1_t0, p0_t1, p1_t1, toi, min_distance, tmax, mode);
                return std::make_tuple(r, toi);
            },
            "p0_t0"_a, "p1_t0"_a, "p0_t1"_a, "p1_t1"_a, "min_distance"_a = 0.0,
            "tmax"_a = 1.0, "mode"_a = NarrowPhaseCCDMode::STANDALONE)
        .def(
            "point_edge_ccd",
            [](const NarrowPhaseCCD& self, Eigen::ConstRef<VectorMax3d> p_t0,
               Eigen::ConstRef<VectorMax3d> e0_t0,
               Eigen::ConstRef<VectorMax3d> e1_t0,
               Eigen::ConstRef<VectorMax3d> p_t1,
               Eigen::ConstRef<VectorMax3d> e0_t1,
               Eigen::ConstRef<VectorMax3d> e1_t1,
               const double min_distance = 0.0, const double tmax = 1.0,
               const NarrowPhaseCCDMode mode = NarrowPhaseCCDMode::STANDALONE) {
                double toi;
                bool r = self.point_edge_ccd(
                    p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi, min_distance,
                    tmax, mode);
                return std::make_tuple(r, toi);
            },
            "p_t0"_a, "e0_t0"_a, "e1_t0"_a, "p_t1"_a, "e0_t1"_a, "e1_t1"_a,
            "min_distance"_a = 0.0, "tmax"_a = 1.0,
            "mode"_a = NarrowPhaseCCDMode::STANDALONE)
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
               const double min_distance = 0.0, const double tmax = 1.0,
               const NarrowPhaseCCDMode mode = NarrowPhaseCCDMode::STANDALONE) {
                double toi;
                bool r = self.point_triangle_ccd(
                    p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi,
                    min_distance, tmax, mode);
                return std::make_tuple(r, toi);
            },
            "p_t0"_a, "t0_t0"_a, "t1_t0"_a, "t2_t0"_a, "p_t1"_a, "t0_t1"_a,
            "t1_t1"_a, "t2_t1"_a, "min_distance"_a = 0.0, "tmax"_a = 1.0,
            "mode"_a = NarrowPhaseCCDMode::STANDALONE)
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
               const double min_distance = 0.0, const double tmax = 1.0,
               const NarrowPhaseCCDMode mode = NarrowPhaseCCDMode::STANDALONE) {
                double toi;
                bool r = self.edge_edge_ccd(
                    ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
                    eb1_t1, toi, min_distance, tmax, mode);
                return std::make_tuple(r, toi);
            },
            "ea0_t0"_a, "ea1_t0"_a, "eb0_t0"_a, "eb1_t0"_a, "ea0_t1"_a,
            "ea1_t1"_a, "eb0_t1"_a, "eb1_t1"_a, "min_distance"_a = 0.0,
            "tmax"_a = 1.0, "mode"_a = NarrowPhaseCCDMode::STANDALONE);
}
