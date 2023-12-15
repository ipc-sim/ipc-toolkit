#include "utils.hpp"

#include <tests/config.hpp>

#include <catch2/catch_test_macros.hpp>

#include <unsupported/Eigen/SparseExtra>

#include <igl/read_triangle_mesh.h>
#include <igl/edges.h>

namespace ipc::tests {

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    const bool success =
        igl::read_triangle_mesh((tests::DATA_DIR / mesh_name).string(), V, F);
    if (F.size()) {
        igl::edges(F, E);
    }
    return success && V.size() && F.size() && E.size();
}

void mmcvids_to_collisions(
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    ipc::Collisions& collisions)
{
    for (int mmcvid_i = 0; mmcvid_i < mmcvids.rows(); mmcvid_i++) {
        const auto mmcvid = mmcvids.row(mmcvid_i);

        if (mmcvid[0] >= 0) { // Is EE?
            int ei, ej;
            // Find the edge index
            for (ei = 0; ei < E.rows(); ei++) {
                if (E(ei, 0) == mmcvid[0] && E(ei, 1) == mmcvid[1]) {
                    break;
                }
            }
            // Find the edge index
            for (ej = 0; ej < E.rows(); ej++) {
                if (E(ej, 0) == mmcvid[2] && E(ej, 1) == mmcvid[3]) {
                    break;
                }
            }
            assert(ei < E.rows() && ej < E.rows());
            collisions.ee_collisions.emplace_back(ei, ej, 0.0);
        } else {
            if (mmcvid[2] < 0) { // Is VV?
                collisions.vv_collisions.emplace_back(
                    -mmcvid[0] - 1, mmcvid[1]);
                assert(-mmcvid[3] >= 1);
                collisions.vv_collisions.back().weight = -mmcvid[3];

            } else if (mmcvid[3] < 0) { // Is EV?
                int ei;
                for (ei = 0; ei < E.rows(); ei++) {
                    if (E(ei, 0) == mmcvid[1] && E(ei, 1) == mmcvid[2]) {
                        break;
                    }
                }
                assert(ei < E.rows());
                collisions.ev_collisions.emplace_back(ei, -mmcvid[0] - 1);
                collisions.ev_collisions.back().weight = -mmcvid[3];

            } else { // Is FV.
                int fi;
                for (fi = 0; fi < F.rows(); fi++) {
                    if (F(fi, 0) == mmcvid[1] && F(fi, 1) == mmcvid[2]
                        && F(fi, 2) == mmcvid[3]) {
                        break;
                    }
                }
                assert(fi < F.rows());
                collisions.fv_collisions.emplace_back(fi, -mmcvid[0] - 1);
            }
        }
    }
}

// Attempts to move the generator to the next element.
// Returns true if successful (and thus has another element that can be
// read)
bool RotationGenerator::next()
{
    double angle = ipc::Vector1d::Random()[0];
    R = Eigen::AngleAxisd(angle, Eigen::Vector3d::Random().normalized());
    return true;
}

// Precondition:
// The generator is either freshly constructed or the last call to next()
// returned true
Eigen::Matrix3d const& RotationGenerator::get() const { return R; }

Catch::Generators::GeneratorWrapper<Eigen::Matrix3d> RotationGenerator::create()
{
    return Catch::Generators::GeneratorWrapper<Eigen::Matrix3d>(
        Catch::Detail::make_unique<RotationGenerator>());
}

///////////////////////////////////////////////////////////////////////////////
// Matrix Market file utils

Eigen::MatrixXd loadMarketXd(const std::string& f)
{
    Eigen::SparseMatrix<double> tmp;
    REQUIRE(Eigen::loadMarket(tmp, f));
    return Eigen::MatrixXd(tmp);
}

Eigen::MatrixXi loadMarketXi(const std::string& f)
{
    Eigen::SparseMatrix<int> tmp;
    REQUIRE(Eigen::loadMarket(tmp, f));
    return Eigen::MatrixXi(tmp);
}

///////////////////////////////////////////////////////////////////////////////

void print_compare_nonzero(
    const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& B,
    bool print_only_different)
{
    fmt::print(
        "A.norm()={}, B.norm()={} (A-B).norm()={}\n", A.norm(), B.norm(),
        (A - B).norm());
    fmt::print("(i,j): A(i,j), B(i,j), abs_diff, rel_diff\n");
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.rows(); j++) {
            const double abs_diff = std::abs(A(i, j) - B(i, j));
            const double rel_diff =
                abs_diff / std::max(std::abs(A(i, j)), std::abs(B(i, j)));

            const double tol =
                std::max({ std::abs(A(i, j)), std::abs(B(i, j)), double(1.0) })
                * 1e-5;

            if ((A(i, j) != 0 || B(i, j) != 0)
                && (!print_only_different || abs_diff > tol)) {
                fmt::print(
                    "({:d},{:d}): {:g}, {:g}, {:g}, {:g}\n", i, j, A(i, j),
                    B(i, j), abs_diff, rel_diff);
            }
        }
    }
    fmt::print("\n");
}

} // namespace ipc::tests