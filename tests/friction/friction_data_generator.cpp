#include "friction_data_generator.hpp"

#include <ipc/config.hpp>

#include <igl/edges.h>
#include <iostream>

Eigen::VectorXd LogSpaced(int num, double start, double stop, double base)
{
    return pow(base, Eigen::VectorXd::LinSpaced(num, start, stop).array());
}

Eigen::VectorXd GeomSpaced(int num, double start, double stop)
{
    return LogSpaced(num, log10(start), log10(stop), /*base=*/10);
}

FrictionData friction_data_generator()
{
    FrictionData data;

    auto& [V0, V1, E, F, constraints, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    mu = GENERATE(range(0.0, 1.0, 0.2));
#ifdef NDEBUG
    epsv_times_h = pow(10, GENERATE(range(-6, 0)));
    dhat = pow(10, GENERATE(range(-4, 0)));
    barrier_stiffness = pow(10, GENERATE(range(0, 2)));
#else
    epsv_times_h = pow(10, GENERATE(range(-6, 0, 2)));
    dhat = pow(10, GENERATE(range(-4, 0, 2)));
    barrier_stiffness = 100;
#endif

    SECTION("point-triangle")
    {
        V0.resize(4, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << 0, d, 0;   // point at t=0
        V0.row(1) << -1, 0, 1;  // triangle vertex 0 at t=0
        V0.row(2) << 2, 0, 0;   // triangle vertex 1 at t=0
        V0.row(3) << -1, 0, -1; // triangle vertex 2 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 1, d + dy, 0; // point at t=1

        F.resize(1, 3);
        F << 1, 2, 3;
        igl::edges(F, E);
        REQUIRE(E.rows() == 3);

        constraints.fv_constraints.emplace_back(0, 0);
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        constraints.fv_constraints.back().weight_gradient.resize(V0.size());
#endif
    }
    SECTION("edge-edge")
    {
        V0.resize(4, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -1, d, 0; // edge a vertex 0 at t=0
        V0.row(1) << 1, d, 0;  // edge a vertex 1 at t=0
        V0.row(2) << 0, 0, -1; // edge b vertex 0 at t=0
        V0.row(3) << 0, 0, 1;  // edge b vertex 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d, 0; // edge a vertex 0 at t=1
        V1.row(1) << 2.5, d, 0; // edge a vertex 1 at t=1

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        constraints.ee_constraints.emplace_back(0, 1, 0.0);
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        constraints.ee_constraints.back().weight_gradient.resize(V0.size());
#endif
    }
    SECTION("point-edge")
    {
        V0.resize(3, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -0.5, d, 0; // point at t=0
        V0.row(1) << 0, 0, -1;   // edge vertex 0 at t=0
        V0.row(2) << 0, 0, 1;    // edge vertex 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d, 0; // point at t=1

        E.resize(1, 2);
        E.row(0) << 1, 2;

        constraints.ev_constraints.emplace_back(0, 1);
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        constraints.ev_constraints.back().weight_gradient.resize(V0.size());
#endif
    }
    SECTION("point-point")
    {
        V0.resize(2, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -1, d, 0; // point 0 at t=0
        V0.row(1) << 1, d, 0;  // point 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d, 0;  // edge a vertex 0 at t=1
        V1.row(1) << -0.5, d, 0; // edge a vertex 1 at t=1

        constraints.vv_constraints.emplace_back(0, 1);
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        constraints.vv_constraints.back().weight_gradient.resize(V0.size());
#endif
    }
    SECTION("point-edge")
    {
        V0.resize(3, 2);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -0.5, d; // point at t=0
        V0.row(1) << -1, 0;   // edge vertex 0 at t=0
        V0.row(2) << 1, 0;    // edge vertex 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d; // point at t=1

        E.resize(1, 2);
        E.row(0) << 1, 2;

        constraints.ev_constraints.emplace_back(0, 1);
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        constraints.ev_constraints.back().weight_gradient.resize(V0.size());
#endif
    }
    SECTION("point-point")
    {
        V0.resize(2, 2);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -1, d; // point 0 at t=0
        V0.row(1) << 1, d;  // point 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d;  // edge a vertex 0 at t=1
        V1.row(1) << -0.5, d; // edge a vertex 1 at t=1

        constraints.vv_constraints.emplace_back(0, 1);
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
        constraints.vv_constraints.back().weight_gradient.resize(V0.size());
#endif
    }

    return data;
}