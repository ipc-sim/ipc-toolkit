#include "friction_data_generator.hpp"

#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/config.hpp>

#include <igl/edges.h>
#include <igl/writeOBJ.h>
#include <ipc/distance/distance_type.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>

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

using namespace ipc;

SmoothFrictionData<3> smooth_friction_data_generator_3d()
{
    SmoothFrictionData<3> data;

    auto& [V0, V1, E, F, collisions, mu, epsv_times_h, param, barrier_stiffness] =
        data;

    double &dhat = param.dhat;
    collisions.set_are_shape_derivatives_enabled(true);

    mu = 1.; // GENERATE(range(0.0, 1.0, 0.2));
#ifdef NDEBUG
    epsv_times_h = pow(10, GENERATE(range(-2, 0)));
    dhat = pow(10, GENERATE(range(-2, 0)));
    barrier_stiffness = 1.;
#else
    epsv_times_h = 1.; // pow(10, GENERATE(range(-6, 0, 2)));
    dhat = 1e-1; // pow(10, GENERATE(range(-4, 0, 2)));
    barrier_stiffness = 1.; // 100;
#endif

    param = ParameterType(dhat, 0.8, 0, 1, 0, 2);
    const double max_d = dhat * 0.9;
    const double min_d = dhat * 0.1;
    const double d = GENERATE_COPY(range(min_d, max_d, max_d / 10));
    SECTION("point-triangle")
    {
        V0.resize(8, 3);
        V0 << 0, d, 0,   // point at t=0
              -1, 0, 1,  // triangle vertex 0 at t=0
              2, 0, 0,   // triangle vertex 1 at t=0
              -1, 0, -1, // triangle vertex 2 at t=0
              -1, d+2, 1,
              2, d+2, 0,
              -1, d+2, -1,
              0, -2, 0;

        V1 = V0;
        Eigen::Vector3d disp;
        disp << 0.5, 0, 0;
        V1.row(0) += disp.transpose();
        V1.row(4) += disp.transpose();
        V1.row(5) += disp.transpose();
        V1.row(6) += disp.transpose();

        F.resize(8, 3);
        F << 1, 2, 3,
             1, 7, 2,
             2, 7, 3,
             3, 7, 1,
             4, 5, 6,
             4, 0, 5,
             5, 0, 6,
             6, 0, 4;
        igl::edges(F, E);

        CollisionMesh mesh(V0, E, F);
        collisions.collisions.push_back(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Face, Point3>>(0, 0, PointTriangleDistanceType::P_T, mesh, param, dhat, V0));
    }
    SECTION("edge-edge")
    {
        V0.resize(10, 3);
        V0.row(0) << -1, d, 0; // edge a vertex 0 at t=0
        V0.row(1) << 1, d, 0;  // edge a vertex 1 at t=0
        V0.row(2) << 0, 0, -1; // edge b vertex 0 at t=0
        V0.row(3) << 0, 0, 1;  // edge b vertex 1 at t=0
        V0.row(4) << 0, d, -1;
        V0.row(5) << 0, d, 1;
        V0.row(6) << 1, 0, 0;
        V0.row(7) << -1, 0, 0;
        V0.row(8) << 0, 2 * d, 0;
        V0.row(9) << 0, -d, 0;

        V1 = V0;
        Eigen::Vector3d disp;
        disp << 0.5, 0, 0;
        V1.row(0) += disp.transpose();
        V1.row(1) += disp.transpose();
        V1.row(4) += disp.transpose();
        V1.row(5) += disp.transpose();
        V1.row(8) += disp.transpose();

        F.resize(12, 3);
        F << 1, 0, 4,
             1, 5, 0,
             1, 4, 8,
             4, 0, 8,
             0, 5, 8,
             5, 1, 8,
             2, 7, 3,
             2, 3, 6,
             2, 6, 9,
             6, 3, 9,
             3, 7, 9,
             7, 2, 9;

        igl::edges(F, E);

        int e0 = E.rows(), e1 = E.rows();
        for (int e = 0; e < E.rows(); e++)
        {
            if (std::min(E(e, 0), E(e, 1)) == 0 &&
                std::max(E(e, 0), E(e, 1)) == 1)
                e0 = e;
            else if (std::min(E(e, 0), E(e, 1)) == 2 &&
                std::max(E(e, 0), E(e, 1)) == 3)
                e1 = e;
        }
        assert(e0 < E.rows());
        assert(e1 < E.rows());

        CollisionMesh mesh(V0, E, F);
        collisions.collisions.push_back(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3>>(e0, e1, EdgeEdgeDistanceType::EA_EB, mesh, param, dhat, V0));
    }
    SECTION("point-edge")
    {
        V0.resize(9, 3);
        V0 << 0, d, 0, // point at t=0
        0, 0, -1,   // edge vertex 0 at t=0
        0, 0, 1,    // edge vertex 1 at t=0
        -1, d * 2, 1,
        2, d * 2, 0,
        -1, d * 2, -1,
        1, 0, 0,
        -1, 0, 0,
        0, -d, 0;

        V1 = V0;
        Eigen::Vector3d disp;
        disp << 0.5, 0, 0;
        V1.row(0) += disp.transpose();
        V1.row(3) += disp.transpose();
        V1.row(4) += disp.transpose();
        V1.row(5) += disp.transpose();

        F.resize(10, 3);
        F << 3, 0, 4,
             4, 0, 5,
             5, 0, 3,
             3, 4, 5,
             1, 7, 2,
             1, 2, 6,
             2, 7, 8,
             7, 1, 8,
             1, 6, 8,
             6, 2, 8;

        igl::edges(F, E);

        int e = 0;
        for (; e < E.rows(); e++)
        {
            if (std::min(E(e, 0), E(e, 1)) == 1 &&
                std::max(E(e, 0), E(e, 1)) == 2)
                break;
        }
        assert(e < E.rows());

        CollisionMesh mesh(V0, E, F);
        collisions.collisions.push_back(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Edge3, Point3>>(e, 0, PointEdgeDistanceType::AUTO, mesh, param, dhat, V0));
    }
    SECTION("point-point")
    {
        V0.resize(8, 3);

        V0 << 0, 0, 0,
            0, d, 0,
            -1, -d, 1,
            2, -d, 0,
            -1, -d, -1,
            -1, 2 * d, 1,
            2, 2 * d, 0,
            -1, 2 * d, -1;

        V1 = V0;
        Eigen::Vector3d disp;
        disp << 0.5, 0, 0;
        V1.row(0) += disp.transpose();
        V1.row(2) += disp.transpose();
        V1.row(3) += disp.transpose();
        V1.row(4) += disp.transpose();

        F.resize(8, 3);
        F << 0, 2, 3,
            0, 3, 4,
            0, 4, 2,
            3, 2, 4,
            1, 6, 5,
            1, 5, 7,
            1, 7, 6,
            5, 6, 7;
        igl::edges(F, E);

        CollisionMesh mesh(V0, E, F);
        collisions.collisions.push_back(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>(0, 1, PointPointDistanceType::AUTO, mesh, param, dhat, V0));
    }

    return data;
}

SmoothFrictionData<2> smooth_friction_data_generator_2d()
{
    SmoothFrictionData<2> data;

    auto& [V0, V1, E, F, collisions, mu, epsv_times_h, param, barrier_stiffness] =
        data;

    double &dhat = param.dhat;
    collisions.set_are_shape_derivatives_enabled(true);

    mu = 1.; // GENERATE(range(0.0, 1.0, 0.2));
#ifdef NDEBUG
    epsv_times_h = pow(10, GENERATE(range(-2, 0)));
    dhat = pow(10, GENERATE(range(-3, 0)));
    barrier_stiffness = 1.;
#else
    epsv_times_h = 1.; // pow(10, GENERATE(range(-6, 0, 2)));
    dhat = 1e-1; // pow(10, GENERATE(range(-4, 0, 2)));
    barrier_stiffness = 1.; // 100;
#endif

    param = ParameterType(dhat, 0.8, 0, 1, 0, 2);
    const double max_d = dhat * 0.9;
    const double min_d = dhat * 0.1;
    const double d = GENERATE_COPY(range(min_d, max_d, max_d / 10));
    SECTION("point-edge 2D")
    {
        V0.resize(6, 2);
        V0 << d, 0, // point at t=0
        0, -1,   // edge vertex 0 at t=0
        0,  1,    // edge vertex 1 at t=0
        1 + d, -1,
        1 + d, 1,
        -1, 0;

        V1 = V0;
        Eigen::Vector2d disp;
        disp << 0, 0.5;
        V1.row(0) += disp.transpose();
        V1.row(3) += disp.transpose();
        V1.row(4) += disp.transpose();

        F.resize(2, 3);
        F << 0, 3, 4,
             1, 2, 5;

        E.resize(6, 2);
        E << 0, 3,
             3, 4,
             4, 0,
             1, 2,
             2, 5,
             5, 1;
        // igl::edges(F, E);

        int e = 3;

        CollisionMesh mesh(V0, E, F);
        collisions.collisions.push_back(std::make_shared<SmoothCollisionTemplate<max_vert_2d, Edge2, Point2>>(e, 0, PointEdgeDistanceType::AUTO, mesh, param, dhat, V0));
    }
    SECTION("point-point 2D")
    {
        V0.resize(6, 2);
        V0 << d, 0,
            0,  0,
            1 + d, -1,
            1 + d, 1,
            0, -1,
            0,  1;

        V1 = V0;
        Eigen::Vector2d disp;
        disp << 0, 0.5;
        V1.row(0) += disp.transpose();
        V1.row(2) += disp.transpose();
        V1.row(3) += disp.transpose();

        F.resize(2, 3);
        F << 0, 2, 3,
             1, 5, 4;

        E.resize(6, 2);
        E << 0, 2,
             2, 3,
             3, 0,
             1, 5,
             5, 4,
             4, 1;

        CollisionMesh mesh(V0, E, F);
        collisions.collisions.push_back(std::make_shared<SmoothCollisionTemplate<max_vert_2d, Point2, Point2>>(0, 1, PointPointDistanceType::AUTO, mesh, param, dhat, V0));
     }

    return data;
}