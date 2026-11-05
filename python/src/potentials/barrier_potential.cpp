#include <common.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/gcp/gcp_potential.hpp>
#include <ipc/esp/esp_potential.hpp>

using namespace ipc;

void define_barrier_potential(py::module_& m)
{
    py::class_<BarrierPotential, NormalPotential>(m, "BarrierPotential")
        .def(
            py::init([](const double dhat, const double stiffness,
                        const bool use_physical_barrier) {
                assert_positive(dhat, "dhat");
                assert_positive(stiffness, "stiffness");
                return BarrierPotential(dhat, stiffness, use_physical_barrier);
            }),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                dhat: The activation distance of the barrier.
                stiffness: The stiffness of the barrier.
                use_physical_barrier: Whether to use the physical barrier.

            Raises:
                ValueError: If dhat or stiffness is not positive.
            )ipc_Qu8mg5v7",
            "dhat"_a, "stiffness"_a, "use_physical_barrier"_a = false)
        .def(
            py::init([](std::shared_ptr<Barrier> barrier, const double dhat,
                        const double stiffness,
                        const bool use_physical_barrier) {
                assert_not_none(barrier, "barrier");
                assert_positive(dhat, "dhat");
                assert_positive(stiffness, "stiffness");
                return BarrierPotential(
                    std::move(barrier), dhat, stiffness, use_physical_barrier);
            }),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                barrier: The barrier function.
                dhat: The activation distance of the barrier.
                stiffness: The stiffness of the barrier.
                use_physical_barrier: Whether to use the physical barrier.

            Raises:
                ValueError: If barrier is None, or dhat or stiffness is not positive.
            )ipc_Qu8mg5v7",
            "barrier"_a, "dhat"_a, "stiffness"_a,
            "use_physical_barrier"_a = false)
        .def_property(
            "dhat", &BarrierPotential::dhat,
            [](BarrierPotential& self, const double dhat) {
                assert_positive(dhat, "dhat");
                self.set_dhat(dhat);
            },
            "Barrier activation distance. Must be positive.")
        .def_property(
            "stiffness", &BarrierPotential::stiffness,
            [](BarrierPotential& self, const double stiffness) {
                assert_positive(stiffness, "stiffness");
                self.set_stiffness(stiffness);
            },
            "Barrier stiffness. Must be positive.")
        .def_property(
            "use_physical_barrier", &BarrierPotential::use_physical_barrier,
            &BarrierPotential::set_use_physical_barrier,
            "Whether to use the physical barrier.")
        .def_property(
            "barrier",
            py::cpp_function(
                &BarrierPotential::barrier, py::return_value_policy::reference),
            [](BarrierPotential& self,
               const std::shared_ptr<Barrier>& barrier) {
                assert_not_none(barrier, "barrier");
                self.set_barrier(barrier);
            },
            "Barrier function used to compute the potential. Must not be None.");
}

void define_smooth_potential(py::module_& m)
{
    py::class_<GCPParameters>(m, "GCPParameters")
        .def(
            py::init<
                const double, const double, const double, const double,
                const double, const int>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for smooth contact.

            Parameters:
                dhat, alpha_t, beta_t, alpha_n, beta_n, r
            )ipc_Qu8mg5v7",
            "dhat"_a, "alpha_t"_a, "beta_t"_a, "alpha_n"_a, "beta_n"_a, "r"_a)
        .def(
            py::init<const double, const double, const double, const int>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for smooth contact.

            Parameters:
                dhat, alpha_t, beta_t, r
            )ipc_Qu8mg5v7",
            "dhat"_a, "alpha_t"_a, "beta_t"_a, "r"_a)
        .def_readonly("dhat", &GCPParameters::dhat)
        .def_readonly("alpha_t", &GCPParameters::alpha_t)
        .def_readonly("beta_t", &GCPParameters::beta_t)
        .def_readonly("alpha_n", &GCPParameters::alpha_n)
        .def_readonly("beta_n", &GCPParameters::beta_n)
        .def_readonly("r", &GCPParameters::r)
        .def_property(
            "adaptive_dhat_ratio", &GCPParameters::adaptive_dhat_ratio,
            &GCPParameters::set_adaptive_dhat_ratio,
            "Ratio of the distance to the interaction set in the rest "
            "configuration used as the per-element adaptive dhat.");

    py::class_<GCPPotential>(m, "GCPPotential")
        .def(
            py::init<const GCPParameters&>(),
            R"ipc_Qu8mg5v7(
            Construct a smooth barrier potential.

            Parameters:
                param: A set of parameters.
            )ipc_Qu8mg5v7",
            "param"_a)
        .def(
            "__call__",
            py::overload_cast<
                const GCPCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &ipc::GCPPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a)
        .def(
            "gradient",
            py::overload_cast<
                const GCPCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &ipc::GCPPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of \|vertices\|.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a)
        .def(
            "hessian",
            py::overload_cast<
                const GCPCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>, const PSDProjectionMethod>(
                &ipc::GCPPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of \|vertices\|x\|vertices\|.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE)
        .def(
            "__call__",
            py::overload_cast<
                const GCPCollision&, Eigen::ConstRef<Eigen::VectorXd>>(
                &ipc::GCPPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a)
        .def(
            "gradient",
            py::overload_cast<
                const GCPCollision&, Eigen::ConstRef<Eigen::VectorXd>>(
                &GCPPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a)
        .def(
            "hessian",
            py::overload_cast<
                const GCPCollision&, Eigen::ConstRef<Eigen::VectorXd>,
                const PSDProjectionMethod>(&GCPPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE);
}

void define_esp_parameters(py::module& m)
{
    py::enum_<EEMollifierType>(
        m, "EEMollifierType", "Shape of the 3D edge-edge mollifier.")
        .value(
            "SMOOTHSTEP_REL", EEMollifierType::SMOOTHSTEP_REL,
            "Smoothstep over a width relative to the gap, w^2 = rel_eps d^2")
        .value(
            "SMOOTHSTEP", EEMollifierType::SMOOTHSTEP,
            "Smoothstep over a fixed width")
        .value(
            "SMOOTHSTEP_SQRT", EEMollifierType::SMOOTHSTEP_SQRT,
            "Smoothstep over w^2 = eps d length")
        .value("GCP", EEMollifierType::GCP, "Original GCP mollifier")
        .value(
            "GCP_FIXED", EEMollifierType::GCP_FIXED,
            "GCP shape with a fixed width")
        .value(
            "GCP_SQRT", EEMollifierType::GCP_SQRT,
            "GCP shape with width^2 = eps d length");

    py::class_<EEMollifierConfig>(
        m, "EEMollifierConfig", "Options of the 3D edge-edge mollifier.")
        .def(py::init())
        .def_readwrite("type", &EEMollifierConfig::type)
        .def_readwrite(
            "rel_eps", &EEMollifierConfig::rel_eps,
            "SMOOTHSTEP_REL: squared width over squared gap.")
        .def_readwrite(
            "width", &EEMollifierConfig::width,
            "SMOOTHSTEP, GCP_FIXED: fixed width (distance units).")
        .def_readwrite(
            "length", &EEMollifierConfig::length,
            "SMOOTHSTEP_SQRT, GCP_SQRT: length scale (distance units).");

    py::class_<ESPParameters> esp_parameters(m, "ESPParameters");

    py::enum_<ESPParameters::IntegrationType>(esp_parameters, "IntegrationType")
        .value(
            "BRUTE_FORCE", ESPParameters::IntegrationType::BRUTE_FORCE,
            "Integrate all pairs with no obstacle filtering")
        .value(
            "NORMAL", ESPParameters::IntegrationType::NORMAL,
            "Filter obstacle-obstacle pairs; skip primitives with only "
            "obstacle candidates")
        .value(
            "NO_OBST", ESPParameters::IntegrationType::NO_OBST,
            "Skip obstacle sources entirely, may miss collisions!");

    esp_parameters
        .def(
            py::init<
                const double, const double, const int, bool,
                ESPParameters::IntegrationType>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for ESP contact.

            Parameters:
                dhat: Activation distance of the barrier.
                dbar_factor: Near/far split parameter dbar / dhat, in [0, 2]. The far part vanishes identically at 2; the edge-edge weights use the support min(dbar, dhat); 0 disables edge-edge contact.
                quad_order: Gauss-Lobatto edge quadrature order (2D).
                area_weights: Whether to weight the quadrature by area.
                integration_type: Which obstacle primitives to integrate.

            Raises:
                ValueError: If dhat is not positive and finite, dbar_factor is not in [0, 2], or quad_order is negative.
            )ipc_Qu8mg5v7",
            py::arg("dhat"), py::arg("dbar_factor") = 0.2,
            py::arg("quad_order") = 1, py::arg("area_weights") = true,
            py::arg_v(
                "integration_type", ESPParameters::IntegrationType::NORMAL,
                "IntegrationType.NORMAL"))
        .def_readonly("dhat", &ESPParameters::dhat)
        .def_readonly("dbar", &ESPParameters::dbar)
        .def_readonly("quad_order", &ESPParameters::quad_order)
        .def_readonly("area_weights", &ESPParameters::area_weights)
        .def_readwrite(
            "local_quadrature_transfer",
            &ESPParameters::local_quadrature_transfer)
        .def_readwrite(
            "legacy_element_weights", &ESPParameters::legacy_element_weights,
            "Weigh faces by area / 9 and edges in no face by length / 4, as "
            "earlier versions did.")
        .def_readwrite(
            "ee_mollifier", &ESPParameters::ee_mollifier,
            "Shape of the edge-edge mollifier.")
        .def_readwrite(
            "ee_mollifier_power", &ESPParameters::ee_mollifier_power,
            "Power of the edge-edge sample weight; 0 chooses it from the "
            "mollifier (1 for smoothstep; round(p) + 1 for GCP under an "
            "inverse-power barrier of power p).")
        .def_readwrite(
            "ee_weight", &ESPParameters::ee_weight,
            "Weight of an edge-edge sample relative to a vertex sample.")
        .def_readwrite(
            "ee_area_scaling", &ESPParameters::ee_area_scaling,
            "Scale edge-edge samples so each stands for about one r x r "
            "patch.")
        .def_readwrite(
            "ee_area_ref", &ESPParameters::ee_area_ref,
            "Reference length r of ee_area_scaling; 0 uses dbar.")
        .def_property_readonly(
            "ee_mollifier_exponent", &ESPParameters::ee_mollifier_exponent,
            "Power applied to the edge-edge mollifier weight.")
        .def_readwrite("barrier", &ESPParameters::barrier)
        .def_readonly("integration_type", &ESPParameters::integration_type)
        .def(
            "set_lattice_rule",
            [](ESPParameters& self, const int m) {
                self.set_quad_rule(ESPParameters::lattice_rule(m));
            },
            R"ipc_Qu8mg5v7(
            Use the composite vertex rule on m x m sub-triangles for every face (3D). m = 1 is the vertex rule.
            )ipc_Qu8mg5v7",
            py::arg("m"))
        .def(
            "set_auto_lattice", &ESPParameters::set_auto_lattice,
            R"ipc_Qu8mg5v7(
            Choose each face's lattice rule from its rest area so the points are about `spacing` apart (3D).

            Parameters:
                spacing: Target point spacing.
                max_subdivisions: Cap on the subdivisions per face.
            )ipc_Qu8mg5v7",
            py::arg("spacing"), py::arg("max_subdivisions") = 16)
        .def_property_readonly(
            "lattice_spacing", &ESPParameters::lattice_spacing,
            "Target spacing of the automatic lattice rules, or 0 for a fixed "
            "rule.");
}

void define_esp_potential(py::module& m)
{
    py::class_<ESPPotential>(m, "ESPPotential")
        .def(
            py::init<const ESPParameters&, const bool>(),
            R"ipc_Qu8mg5v7(
            Construct an Extremum Sum Potential (ESP).

            Parameters:
                param: A set of parameters.
                use_near_far: Whether to normalize the quadrature weights so they sum to 1, splitting near and far contributions when 0 < dbar_factor <= 2 (the far part vanishes at 2).
            )ipc_Qu8mg5v7",
            py::arg("param"), py::arg("use_near_far") = true)
        .def(
            "element_near_weights", &ESPPotential::element_near_weights,
            "Each integration element's near-weight total at the given "
            "vertices (see set_lagged_weights).",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "set_lagged_weights", &ESPPotential::set_lagged_weights,
            "Hold each element's normalization fixed at these weights.",
            py::arg("weights"))
        .def("clear_lagged_weights", &ESPPotential::clear_lagged_weights)
        .def_property_readonly("lagged_weights", &ESPPotential::lagged_weights)
        .def(
            "__call__",
            py::overload_cast<
                const ESPCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &ipc::ESPPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "gradient",
            py::overload_cast<
                const ESPCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &ipc::ESPPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of \|vertices\|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "hessian",
            py::overload_cast<
                const ESPCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>, const PSDProjectionMethod>(
                &ipc::ESPPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite. With normalized weights and edge-edge terms, each face's block is projected as a whole.

            Returns:
                The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of \|vertices\|x\|vertices\|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"),
            py::arg_v(
                "project_hessian_to_psd", PSDProjectionMethod::NONE,
                "PSDProjectionMethod.NONE"));
}
