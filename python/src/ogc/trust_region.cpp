#include <common.hpp>

#include <ipc/ogc/trust_region.hpp>

using namespace ipc;

void define_trust_region(py::module_& m)
{
    py::class_<ogc::TrustRegion>(
        m, "TrustRegion",
        R"ipc_Qu8mg5v7(
         A trust region for filtering optimization steps.
         See "Offset Geometric Contact" by Chen et al. [2025]
         )ipc_Qu8mg5v7")
        .def(
            py::init<double>(),
            R"ipc_Qu8mg5v7(
            Construct a TrustRegion object.

            Parameters:
                dhat: The offset distance for contact.
            )ipc_Qu8mg5v7",
            "dhat"_a)
        .def(
            "warm_start_time_step", &ogc::TrustRegion::warm_start_time_step,
            R"ipc_Qu8mg5v7(
            Warm start the time step by moving towards the predicted positions.

            This also initializes the trust region.

            Parameters:
                mesh: The collision mesh.
                x: Current vertex positions. (Will be modified)
                pred_x: Predicted vertex positions.
                collisions: Collisions to be initialized.
                dhat: The offset distance for contact.
                min_distance: Minimum distance between elements.
                broad_phase: Broad phase collision detection.
            )ipc_Qu8mg5v7",
            "mesh"_a, "x"_a, "pred_x"_a, "collisions"_a, "dhat"_a,
            "min_distance"_a = 0.0,
            "broad_phase"_a = make_default_broad_phase())
        .def(
            "update", &ogc::TrustRegion::update,
            R"ipc_Qu8mg5v7(
            Update the trust region based on the current positions.

            Parameters:
                mesh: The collision mesh.
                x: Current vertex positions.
                collisions: Collisions to be updated.
                min_distance: Minimum distance between elements.
                broad_phase: Broad phase collision detection.
            )ipc_Qu8mg5v7",
            "mesh"_a, "x"_a, "collisions"_a, "min_distance"_a = 0.0,
            "broad_phase"_a = make_default_broad_phase())
        .def(
            "update_if_needed", &ogc::TrustRegion::update_if_needed,
            R"ipc_Qu8mg5v7(
            Update the trust region if needed based on the current positions.

            Parameters:
                mesh: The collision mesh.
                x: Current vertex positions.
                collisions: Collisions to be updated.
                min_distance: Minimum distance between elements.
                broad_phase: Broad phase collision detection.
            )ipc_Qu8mg5v7",
            "mesh"_a, "x"_a, "collisions"_a, "min_distance"_a = 0.0,
            "broad_phase"_a = make_default_broad_phase())
        .def(
            "filter_step", &ogc::TrustRegion::filter_step,
            R"ipc_Qu8mg5v7(
            Filter the optimization step dx to stay within the trust region.

            Parameters:
                mesh: The collision mesh.
                x: Current vertex positions.
                dx: Proposed vertex displacements.

            Note:
                Sets should_update_trust_region to true if the trust region should be updated on the next iteration.
            )ipc_Qu8mg5v7",
            "mesh"_a, "x"_a, "dx"_a)
        .def_readwrite(
            "trust_region_centers", &ogc::TrustRegion::trust_region_centers,
            "Centers of the trust regions for each vertex.")
        .def_readwrite(
            "trust_region_radii", &ogc::TrustRegion::trust_region_radii,
            "Radii of the trust regions for each vertex.")
        .def_readwrite(
            "relaxed_radius_scaling", &ogc::TrustRegion::relaxed_radius_scaling,
            R"ipc_Qu8mg5v7(
            Scaling factor for the relaxed trust region radius.

            Note:
                - This should be in (0, 1).
                - This is referred to as :math:`2\gamma_p` in the paper.
            )ipc_Qu8mg5v7")
        .def_readwrite(
            "update_threshold", &ogc::TrustRegion::update_threshold,
            R"ipc_Qu8mg5v7(
            Threshold for updating the trust region.

            If more than this fraction of vertices are restricted by the trust region, we update the trust region.

            Note
            ----
                This is referred to as \f\(\gamma_e\f\) in the paper.
            )ipc_Qu8mg5v7")
        .def_readwrite(
            "trust_region_inflation_radius",
            &ogc::TrustRegion::trust_region_inflation_radius,
            R"ipc_Qu8mg5v7(
            Inflation radius for the trust region.

            This is computed each time step based on the predicted motion.
            )ipc_Qu8mg5v7")
        .def_readwrite(
            "should_update_trust_region",
            &ogc::TrustRegion::should_update_trust_region,
            "If true, the trust region will be updated on the next call to ``update_if_needed``.");
}
