#include "friction_data_generator.hpp"

#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/config.hpp>

#include <igl/edges.h>

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

    auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    collisions.set_are_shape_derivatives_enabled(true);

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

    const double max_d = dhat - 2e-8;
    const double d = GENERATE_COPY(range(0.0, max_d, max_d / 10));
    SECTION("point-triangle")
    {
        V0.resize(4, 3);
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

        collisions.fv_collisions.emplace_back(0, 0);
        collisions.fv_collisions.back().weight_gradient.resize(V0.size());
    }
    SECTION("edge-edge")
    {
        V0.resize(4, 3);
        V0.row(0) << -1, d, 0; // edge a vertex 0 at t=0
        V0.row(1) << 1, d, 0;  // edge a vertex 1 at t=0
        V0.row(2) << 0, 0, -1; // edge b vertex 0 at t=0
        V0.row(3) << 0, 0, 1;  // edge b vertex 1 at t=0

        V1 = V0;
        // double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d, 0; // edge a vertex 0 at t=1
        V1.row(1) << 2.5, d, 0; // edge a vertex 1 at t=1

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        collisions.ee_collisions.emplace_back(0, 1, 0.0);
        collisions.ee_collisions.back().weight_gradient.resize(V0.size());
    }
    SECTION("point-edge")
    {
        V0.resize(3, 3);
        V0.row(0) << -0.5, d, 0; // point at t=0
        V0.row(1) << 0, 0, -1;   // edge vertex 0 at t=0
        V0.row(2) << 0, 0, 1;    // edge vertex 1 at t=0

        V1 = V0;
        // double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d, 0; // point at t=1

        E.resize(1, 2);
        E.row(0) << 1, 2;

        collisions.ev_collisions.emplace_back(0, 1);
        collisions.ev_collisions.back().weight_gradient.resize(V0.size());
    }
    SECTION("point-point")
    {
        V0.resize(2, 3);
        V0.row(0) << -1, d, 0; // point 0 at t=0
        V0.row(1) << 1, d, 0;  // point 1 at t=0

        V1 = V0;
        // double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d, 0;  // edge a vertex 0 at t=1
        V1.row(1) << -0.5, d, 0; // edge a vertex 1 at t=1

        collisions.vv_collisions.emplace_back(0, 1);
        collisions.vv_collisions.back().weight_gradient.resize(V0.size());
    }
    SECTION("point-edge 2D")
    {
        V0.resize(3, 2);
        V0.row(0) << -0.5, d; // point at t=0
        V0.row(1) << -1, 0;   // edge vertex 0 at t=0
        V0.row(2) << 1, 0;    // edge vertex 1 at t=0

        V1 = V0;
        // double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d; // point at t=1

        E.resize(1, 2);
        E.row(0) << 1, 2;

        collisions.ev_collisions.emplace_back(0, 1);
        collisions.ev_collisions.back().weight_gradient.resize(V0.size());
    }
    SECTION("point-point 2D")
    {
        V0.resize(2, 2);
        V0.row(0) << -1, d; // point 0 at t=0
        V0.row(1) << 1, d;  // point 1 at t=0

        V1 = V0;
        // double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 0.5, d;  // edge a vertex 0 at t=1
        V1.row(1) << -0.5, d; // edge a vertex 1 at t=1

        collisions.vv_collisions.emplace_back(0, 1);
        collisions.vv_collisions.back().weight_gradient.resize(V0.size());
    }

    return data;
}