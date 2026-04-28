#pragma once

#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <nlohmann/json.hpp>

#include <string>

namespace ipc::tests {

class BroadPhaseGenerator
    : public Catch::Generators::IGenerator<std::shared_ptr<BroadPhase>> {
public:
    BroadPhaseGenerator();

    // Attempts to move the generator to the next element.
    // Returns true if successful (and thus has another element that can be
    // read)
    bool next() override;

    // Precondition:
    // The generator is either freshly constructed or the last call to next()
    // returned true
    std::shared_ptr<BroadPhase> const& get() const override;

    static Catch::Generators::GeneratorWrapper<std::shared_ptr<BroadPhase>>
    create();

private:
    std::vector<std::shared_ptr<BroadPhase>> m_broad_phases;
    int m_current = 0;
};

std::vector<std::shared_ptr<BroadPhase>> broad_phases();

// ============================================================================

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F);

// ============================================================================

void mmcvids_to_collisions(
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    ipc::NormalCollisions& collisions);

// ============================================================================
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

// ============================================================================
// Matrix Market file utils

Eigen::MatrixXd loadMarketXd(const std::string& f);
Eigen::MatrixXi loadMarketXi(const std::string& f);

// ============================================================================

void print_compare_nonzero(
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& B,
    bool print_only_different = true);

// ============================================================================

inline Eigen::Vector2d edge_normal(
    Eigen::ConstRef<Eigen::Vector2d> e0, Eigen::ConstRef<Eigen::Vector2d> e1)
{
    Eigen::Vector2d e = e1 - e0;
    Eigen::Vector2d normal(-e.y(), e.x());
    return normal.normalized();
}

} // namespace ipc::tests

namespace nlohmann {

// Based on https://github.com/nlohmann/json/discussions/3330
template <
    typename Scalar,
    int Rows,
    int Cols,
    int Options,
    int MaxRows,
    int MaxCols>
struct adl_serializer<
    Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
private:
    using Matrix = Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>;
    using Index = Eigen::Index;

public:
    static void to_json(nlohmann::json& j, const Matrix& M)
    {
        for (Index r = 0; r < M.rows(); ++r) {
            if (M.cols() > 1) {
                nlohmann::json jrow = nlohmann::json::array();
                for (Index c = 0; c < M.cols(); ++c) {
                    jrow.push_back(M(r, c));
                }
                j.emplace_back(std::move(jrow));
            } else {
                j.push_back(M(r, 0));
            }
        }
    }

    static void from_json(const nlohmann::json& j, Matrix& M)
    {
        bool resized = false;
        auto resize = [&](size_t nrows, size_t ncols) {
            if (!resized) {
                M.resize(nrows, ncols);
                resized = true;
            } else {
                assert(static_cast<size_t>(M.rows()) == nrows);
                assert(static_cast<size_t>(M.cols()) == ncols);
            }
        };

        for (size_t r = 0, nrows = j.size(); r < nrows; ++r) {
            const auto& jrow = j.at(r);
            if (jrow.is_array()) {
                const size_t ncols = jrow.size();
                resize(nrows, ncols);
                for (size_t c = 0; c < ncols; ++c) {
                    const auto& value = jrow.at(c);
                    M(static_cast<Index>(r), static_cast<Index>(c)) =
                        value.get<Scalar>();
                }
            } else {
                resize(nrows, 1);
                M(static_cast<Index>(r), 0) = jrow.get<Scalar>();
            }
        }
    }
};
} // namespace nlohmann