#include "friction_data_generator.hpp"

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

template <typename Vector> int get_val(const Vector& vals, int& i)
{
    double r = vals[i % vals.size()];
    i /= vals.size();
    return r;
}

FrictionDataGenerator::FrictionDataGenerator()
{
    case_i = 0;
    next();
}

// Attempts to move the generator to the next element.
// Returns true if successful (and thus has another element that can be
// read)
bool FrictionDataGenerator::next()
{
    const int TOTAL_NUM_CASES = MU.size() * EPSV_TIMES_H.size() * DHAT.size()
        * BARRIER_STIFFNESS.size() * DISTANCE_MULTIPLIER.size() * N_CASES;
    if (case_i >= TOTAL_NUM_CASES) {
        return false;
    }

    int i = case_i;
    data.mu = get_val(MU, i);
    data.epsv_times_h = get_val(EPSV_TIMES_H, i);
    data.dhat = get_val(DHAT, i);
    data.barrier_stiffness = get_val(BARRIER_STIFFNESS, i);
    const double d = get_val(DISTANCE_MULTIPLIER, i) * data.dhat;
    const double dy = 0.0; // get_val(DISP_Y, i);
    assert(i >= 0 && i < N_CASES);

    data.V0.resize(0, 0);
    data.V1.resize(0, 0);
    data.E.resize(0, 0);
    data.F.resize(0, 0);
    data.constraints.clear();

    switch (i) {
    case 0: { // point-triangle
        data.V0.resize(4, 3);
        data.V0.row(0) << 0, d, 0;   // point at t=0
        data.V0.row(1) << -1, 0, 1;  // triangle vertex 0 at t=0
        data.V0.row(2) << 2, 0, 0;   // triangle vertex 1 at t=0
        data.V0.row(3) << -1, 0, -1; // triangle vertex 2 at t=0

        data.V1 = data.V0;
        data.V1.row(0) << 1, d + dy, 0; // point at t=1

        data.F.resize(1, 3);
        data.F << 1, 2, 3;
        igl::edges(data.F, data.E);
        REQUIRE(data.E.rows() == 3);

        data.constraints.fv_constraints.emplace_back(0, 0);
        data.constraints.fv_constraints.back().weight_gradient.resize(
            data.V0.size());
    } break;
    case 1: { // edge-edge
        data.V0.resize(4, 3);
        data.V0.row(0) << -1, d, 0; // edge a vertex 0 at t=0
        data.V0.row(1) << 1, d, 0;  // edge a vertex 1 at t=0
        data.V0.row(2) << 0, 0, -1; // edge b vertex 0 at t=0
        data.V0.row(3) << 0, 0, 1;  // edge b vertex 1 at t=0

        data.V1 = data.V0;
        data.V1.row(0) << 0.5, d, 0; // edge a vertex 0 at t=1
        data.V1.row(1) << 2.5, d, 0; // edge a vertex 1 at t=1

        data.E.resize(2, 2);
        data.E.row(0) << 0, 1;
        data.E.row(1) << 2, 3;

        data.constraints.ee_constraints.emplace_back(0, 1, 0.0);
        data.constraints.ee_constraints.back().weight_gradient.resize(
            data.V0.size());
    } break;
    case 2: { // point-edge
        data.V0.resize(3, 3);
        data.V0.row(0) << -0.5, d, 0; // point at t=0
        data.V0.row(1) << 0, 0, -1;   // edge vertex 0 at t=0
        data.V0.row(2) << 0, 0, 1;    // edge vertex 1 at t=0

        data.V1 = data.V0;
        data.V1.row(0) << 0.5, d, 0; // point at t=1

        data.E.resize(1, 2);
        data.E.row(0) << 1, 2;

        data.constraints.ev_constraints.emplace_back(0, 1);
        data.constraints.ev_constraints.back().weight_gradient.resize(
            data.V0.size());
    } break;
    case 3: { // point-point
        data.V0.resize(2, 3);
        data.V0.row(0) << -1, d, 0; // point 0 at t=0
        data.V0.row(1) << 1, d, 0;  // point 1 at t=0

        data.V1 = data.V0;
        data.V1.row(0) << 0.5, d, 0;  // edge a vertex 0 at t=1
        data.V1.row(1) << -0.5, d, 0; // edge a vertex 1 at t=1

        data.constraints.vv_constraints.emplace_back(0, 1);
        data.constraints.vv_constraints.back().weight_gradient.resize(
            data.V0.size());
    } break;
    case 4: { // point-edge
        data.V0.resize(3, 2);
        data.V0.row(0) << -0.5, d; // point at t=0
        data.V0.row(1) << -1, 0;   // edge vertex 0 at t=0
        data.V0.row(2) << 1, 0;    // edge vertex 1 at t=0

        data.V1 = data.V0;
        data.V1.row(0) << 0.5, d; // point at t=1

        data.E.resize(1, 2);
        data.E.row(0) << 1, 2;

        data.constraints.ev_constraints.emplace_back(0, 1);
        data.constraints.ev_constraints.back().weight_gradient.resize(
            data.V0.size());
    } break;
    case 5: { // point-point
        data.V0.resize(2, 2);
        data.V0.row(0) << -1, d; // point 0 at t=0
        data.V0.row(1) << 1, d;  // point 1 at t=0

        data.V1 = data.V0;
        data.V1.row(0) << 0.5, d;  // edge a vertex 0 at t=1
        data.V1.row(1) << -0.5, d; // edge a vertex 1 at t=1

        data.constraints.vv_constraints.emplace_back(0, 1);
        data.constraints.vv_constraints.back().weight_gradient.resize(
            data.V0.size());
    } break;
    default:
        throw "invalid i";
    }

    case_i++;
    return true;
}

// Precondition:
// The generator is either freshly constructed or the last call to next()
// returned true
FrictionData const& FrictionDataGenerator::get() const { return data; }

Catch::Generators::GeneratorWrapper<FrictionData>
FrictionDataGenerator::create()
{
    return Catch::Generators::GeneratorWrapper<FrictionData>(
        std::make_unique<FrictionDataGenerator>());
}
