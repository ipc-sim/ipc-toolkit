#include "friction_data_generator.hpp"

#include <catch2/generators/catch_generators_range.hpp>
#include <ipc/config.hpp>
#include <igl/edges.h>

// Function to generate logarithmically spaced values
Eigen::VectorXd LogSpaced(int num, double start, double stop, double base)
{
    return pow(base, Eigen::VectorXd::LinSpaced(num, start, stop).array());
}

// Function to generate geometrically spaced values
Eigen::VectorXd GeomSpaced(int num, double start, double stop)
{
    return LogSpaced(num, log10(start), log10(stop), 10.0);
}

// Generates simple friction test data
FrictionSimpleData friction_data_generator()
{
    FrictionSimpleData data;
    auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] = data;

    collisions.set_enable_shape_derivatives(true);

    // Generate static friction coefficient
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
        V0 << 0, d, 0,    // Point at t=0
              -1, 0, 1,   // Triangle vertex 0 at t=0
              2, 0, 0,    // Triangle vertex 1 at t=0
              -1, 0, -1;  // Triangle vertex 2 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 1, d + dy, 0; // Point at t=1

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
        V0 << -1, d, 0, // Edge A vertex 0 at t=0
              1, d, 0,  // Edge A vertex 1 at t=0
              0, 0, -1, // Edge B vertex 0 at t=0
              0, 0, 1;  // Edge B vertex 1 at t=0

        V1 = V0;
        V1.row(0) << 0.5, d, 0; // Edge A vertex 0 at t=1
        V1.row(1) << 2.5, d, 0; // Edge A vertex 1 at t=1

        E.resize(2, 2);
        E << 0, 1,
             2, 3;

        collisions.ee_collisions.emplace_back(0, 1, 0.0);
        collisions.ee_collisions.back().weight_gradient.resize(V0.size());
    }

    SECTION("point-edge")
    {
        V0.resize(3, 3);
        V0 << -0.5, d, 0, // Point at t=0
              0, 0, -1,   // Edge vertex 0 at t=0
              0, 0, 1;    // Edge vertex 1 at t=0

        V1 = V0;
        V1.row(0) << 0.5, d, 0; // Point at t=1

        E.resize(1, 2);
        E.row(0) << 1, 2;

        collisions.ev_collisions.emplace_back(0, 1);
        collisions.ev_collisions.back().weight_gradient.resize(V0.size());
    }

    SECTION("point-point")
    {
        V0.resize(2, 3);
        V0 << -1, d, 0, // Point 0 at t=0
              1, d, 0;  // Point 1 at t=0

        V1 = V0;
        V1 << 0.5, d, 0,   // Point 0 at t=1
              -0.5, d, 0; // Point 1 at t=1

        collisions.vv_collisions.emplace_back(0, 1);
        collisions.vv_collisions.back().weight_gradient.resize(V0.size());
    }

    SECTION("point-edge 2D")
    {
        V0.resize(3, 2);
        V0 << -0.5, d, // Point at t=0
              -1, 0,   // Edge vertex 0 at t=0
              1, 0;    // Edge vertex 1 at t=0

        V1 = V0;
        V1.row(0) << 0.5, d; // Point at t=1

        E.resize(1, 2);
        E.row(0) << 1, 2;

        collisions.ev_collisions.emplace_back(0, 1);
        collisions.ev_collisions.back().weight_gradient.resize(V0.size());
    }

    SECTION("point-point 2D")
    {
        V0.resize(2, 2);
        V0 << -1, d, // Point 0 at t=0
              1, d;  // Point 1 at t=0

        V1 = V0;
        V1 << 0.5, d,  // Point 0 at t=1
              -0.5, d; // Point 1 at t=1

        collisions.vv_collisions.emplace_back(0, 1);
        collisions.vv_collisions.back().weight_gradient.resize(V0.size());
    }

    return data;
}

// // Generates complex friction test data with pairwise friction handling
// FrictionComplexData friction_data_generator_with_pairwise()
// {
//     FrictionComplexData data;
//     auto& [V0, V1, E, F, collisions, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness, pairwise_friction] = data;

//     collisions.set_enable_shape_derivatives(true);

//     // Set default static and kinetic friction coefficients
//     static_mu = 0.7; // Global default static friction
//     kinetic_mu = 0.5; // Global default kinetic friction

// #ifdef NDEBUG
//     epsv_times_h = pow(10, GENERATE(range(-6, 0)));
//     dhat = pow(10, GENERATE(range(-4, 0)));
//     barrier_stiffness = pow(10, GENERATE(range(0, 2)));
// #else
//     epsv_times_h = pow(10, GENERATE(range(-6, 0, 2)));
//     dhat = pow(10, GENERATE(range(-4, 0, 2)));
//     barrier_stiffness = 100;
// #endif

//     // Define pairwise friction values for different material pairs
//     pairwise_friction[std::make_tuple(1, 2)] = {0.6, 0.4};  // Pairwise static/kinetic friction for material 1-2
//     pairwise_friction[std::make_tuple(2, 3)] = {0.8, 0.6};  // Pairwise static/kinetic friction for material 2-3
//     pairwise_friction[std::make_tuple(1, 3)] = {0.9, 0.7};  // Pairwise static/kinetic friction for material 1-3

//     const double max_d = dhat - 2e-8;
//     const double d = GENERATE_COPY(range(0.0, max_d, max_d / 10));

//     // Utility lambda to configure collisions
//     auto configure_collisions = [&](const std::vector<std::pair<int, int>>& edges, const Eigen::MatrixXi& faces) {
//         if (!edges.empty()) {
//             E.resize(edges.size(), 2);
//             for (size_t i = 0; i < edges.size(); ++i) {
//                 E.row(i) << edges[i].first, edges[i].second;
//             }
//         }
//         F = faces;
//         igl::edges(F, E);
//     };

//     // Different friction collision scenarios
//     SECTION("point-triangle pairwise")
//     {
//         V0.resize(4, 3);
//         V0.row(0) << 0, d, 0;   // Point at t=0
//         V0.row(1) << -1, 0, 1;  // Triangle vertex 0 at t=0
//         V0.row(2) << 2, 0, 0;   // Triangle vertex 1 at t=0
//         V0.row(3) << -1, 0, -1; // Triangle vertex 2 at t=0

//         V1 = V0;
//         double dy = GENERATE(-1, 1, 1e-1);
//         V1.row(0) << 1, d + dy, 0; // Point at t=1

//         // Faces for the point-triangle scenario
//         Eigen::MatrixXi faces(1, 3);
//         faces << 1, 2, 3;
//         configure_collisions({}, faces);

//         REQUIRE(E.rows() == 3); // Ensure edges are correctly generated

//         collisions.fv_collisions.emplace_back(0, 0);
//         collisions.fv_collisions.back().weight_gradient.resize(V0.size());
//     }

//     SECTION("edge-edge pairwise")
//     {
//         V0.resize(4, 3);
//         V0.row(0) << -1, d, 0; // Edge A vertex 0 at t=0
//         V0.row(1) << 1, d, 0;  // Edge A vertex 1 at t=0
//         V0.row(2) << 0, 0, -1; // Edge B vertex 0 at t=0
//         V0.row(3) << 0, 0, 1;  // Edge B vertex 1 at t=0

//         V1 = V0;
//         V1.row(0) << 0.5, d, 0; // Edge A vertex 0 at t=1
//         V1.row(1) << 2.5, d, 0; // Edge A vertex 1 at t=1

//         // Edges for the edge-edge scenario
//         std::vector<std::pair<int, int>> edges = {{0, 1}, {2, 3}};
//         configure_collisions(edges, Eigen::MatrixXi());

//         collisions.ee_collisions.emplace_back(0, 1, 0.0);
//         collisions.ee_collisions.back().weight_gradient.resize(V0.size());
//     }

//     SECTION("point-edge pairwise")
//     {
//         V0.resize(3, 3);
//         V0.row(0) << -0.5, d, 0; // Point at t=0
//         V0.row(1) << 0, 0, -1;   // Edge vertex 0 at t=0
//         V0.row(2) << 0, 0, 1;    // Edge vertex 1 at t=0

//         V1 = V0;
//         V1.row(0) << 0.5, d, 0; // Point at t=1

//         // Edges for point-edge scenario
//         std::vector<std::pair<int, int>> edges = {{1, 2}};
//         configure_collisions(edges, Eigen::MatrixXi());

//         collisions.ev_collisions.emplace_back(0, 1);
//         collisions.ev_collisions.back().weight_gradient.resize(V0.size());
//     }

//     SECTION("point-point pairwise")
//     {
//         V0.resize(2, 3);
//         V0.row(0) << -1, d, 0; // Point 0 at t=0
//         V0.row(1) << 1, d, 0;  // Point 1 at t=0

//         V1 = V0;
//         V1.row(0) << 0.5, d, 0;  // Point 0 at t=1
//         V1.row(1) << -0.5, d, 0; // Point 1 at t=1

//         // Point-point collision
//         collisions.vv_collisions.emplace_back(0, 1);
//         collisions.vv_collisions.back().weight_gradient.resize(V0.size());
//     }

//     return data;
// }
