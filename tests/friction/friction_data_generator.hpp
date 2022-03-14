#pragma once

#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <ipc/friction/friction_constraint.hpp>

struct FrictionData {
    Eigen::MatrixXd V0;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    ipc::Constraints constraints;
    double mu;
    double epsv_times_h;
    double dhat;
    double barrier_stiffness;
};

Eigen::VectorXd LogSpaced(int num, double start, double stop, double base = 10);
Eigen::VectorXd GeomSpaced(int num, double start, double stop);

static const int N_CASES = 6;
static const Eigen::VectorXd MU = Eigen::VectorXd::LinSpaced(11, 0, 1);

#ifdef NDEBUG
static const Eigen::VectorXd EPSV_TIMES_H = GeomSpaced(7, 1e-6, 1);
static const Eigen::VectorXd DHAT = GeomSpaced(4, 1e-4, 1e-1);
static const Eigen::VectorXd BARRIER_STIFFNESS = GeomSpaced(3, 1, 100);
#else
static const Eigen::VectorXd EPSV_TIMES_H = GeomSpaced(4, 1e-6, 1);
static const Eigen::VectorXd DHAT = GeomSpaced(2, 1e-4, 1e-1);
static const Eigen::VectorXd BARRIER_STIFFNESS = GeomSpaced(2, 1, 100);
#endif

static const Eigen::VectorXd DISTANCE_MULTIPLIER = GeomSpaced(6, 2e-5, 2);
static const Eigen::VectorXd DISP_Y = Eigen::VectorXd::LinSpaced(21, -1, 1);

class FrictionDataGenerator
    : public Catch::Generators::IGenerator<FrictionData> {
public:
    FrictionDataGenerator();

    // Attempts to move the generator to the next element.
    // Returns true if successful (and thus has another element that can be
    // read)
    bool next() override;

    // Precondition:
    // The generator is either freshly constructed or the last call to next()
    // returned true
    FrictionData const& get() const override;

    static Catch::Generators::GeneratorWrapper<FrictionData> create();

protected:
    FrictionData data;
    int case_i;
};
