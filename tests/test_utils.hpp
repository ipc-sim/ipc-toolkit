#pragma once

#include <catch2/catch.hpp>

#include <string>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <ipc/collisions/constraints.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/utils/eigen_ext.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA
#define NUM_BROAD_PHASE_METHODS static_cast<int>(BroadPhaseMethod::NUM_METHODS)
#else
#define NUM_BROAD_PHASE_METHODS                                                \
    (static_cast<int>(BroadPhaseMethod::NUM_METHODS) - 1)
#endif

#define GENERATE_BROAD_PHASE_METHODS()                                         \
    static_cast<BroadPhaseMethod>(GENERATE(range(0, NUM_BROAD_PHASE_METHODS)));

///////////////////////////////////////////////////////////////////////////////

static const std::string TEST_DATA_DIR(TEST_DATA_DIR_CSTR);

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F);

///////////////////////////////////////////////////////////////////////////////

void mmcvids_to_constraints(
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    ipc::Constraints& constraints);

///////////////////////////////////////////////////////////////////////////////
// Rotation generator

class RotationGenerator
    : public Catch::Generators::IGenerator<Eigen::Matrix3d> {
public:
    // Attempts to move the generator to the next element.
    // Returns true if successful (and thus has another element that can be
    // read)
    bool next() override;

    // Precondition:
    // The generator is either freshly constructed or the last call to next()
    // returned true
    Eigen::Matrix3d const& get() const override;

    static Catch::Generators::GeneratorWrapper<Eigen::Matrix3d> create();

protected:
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
};

///////////////////////////////////////////////////////////////////////////////
// JSON utils

template <typename T>
inline void from_json(const nlohmann::json& json, ipc::VectorX<T>& vec)
{
    vec = Eigen::Map<ipc::VectorX<T>>(
        json.get<std::vector<T>>().data(), json.size());
}

template <typename T>
void from_json(const nlohmann::json& json, ipc::MatrixX<T>& mat)
{
    typedef std::vector<std::vector<T>> L;
    L list = json.get<L>();

    size_t num_rows = list.size();
    if (num_rows == 0) {
        return;
    }
    size_t num_cols = list[0].size();
    mat.resize(num_rows, num_cols);

    for (size_t i = 0; i < num_rows; ++i) {
        assert(num_cols == list[i].size());
        mat.row(i) = Eigen::Map<ipc::RowVectorX<T>>(list[i].data(), num_cols);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Matrix Market file utils

Eigen::MatrixXd loadMarketXd(const std::string& f);
Eigen::MatrixXi loadMarketXi(const std::string& f);

///////////////////////////////////////////////////////////////////////////////

void print_compare_nonzero(
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& B,
    bool print_only_different = true);

///////////////////////////////////////////////////////////////////////////////

inline Eigen::Vector2d
edge_normal(const Eigen::Vector2d& e0, const Eigen::Vector2d& e1)
{
    Eigen::Vector2d e = e1 - e0;
    Eigen::Vector2d normal(-e.y(), e.x());
    return normal.normalized();
}