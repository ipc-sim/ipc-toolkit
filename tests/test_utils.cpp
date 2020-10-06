#include "test_utils.hpp"

#include <catch2/catch.hpp>

#include <igl/dirname.h>
#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

// Flatten the matrix rowwise
Eigen::VectorXd flatten(const Eigen::MatrixXd& X)
{
    Eigen::MatrixXd XT = X.transpose();
    return Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(XT.data(), XT.size()));
}

/// Unflatten rowwise
Eigen::MatrixXd unflatten(const Eigen::VectorXd& x, int dim)
{
    assert(x.size() % dim == 0);
    Eigen::MatrixXd unflat_x(x.size() / dim, dim);
    for (int i = 0; i < x.size(); i++) {
        unflat_x(i / dim, i % dim) = x(i);
    }
    return unflat_x;
}

TEST_CASE("Flatten and unflatten", "[utils]")
{
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(1000, 3);
    Eigen::MatrixXd R = unflatten(flatten(X), X.cols());
    CHECK(X == R);
}

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    bool success = igl::read_triangle_mesh(TEST_DATA_DIR + mesh_name, V, F);
    if (F.size()) {
        igl::edges(F, E);
    }
    return success && V.size() && F.size() && E.size();
}

void mmcvids_to_constraints(
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& mmcvids,
    ipc::Constraints& constraints)
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
            constraints.ee_constraints.emplace_back(ei, ej, 0.0);
        } else {
            if (mmcvid[2] < 0) { // Is VV?
                constraints.vv_constraints.emplace_back(
                    -mmcvid[0] - 1, mmcvid[1]);
                assert(-mmcvid[3] >= 1);
                constraints.vv_constraints.back().multiplicity = -mmcvid[3];

            } else if (mmcvid[3] < 0) { // Is EV?

                for (int i = 0; i < E.rows(); i++) {
                    if (E(i, 0) == mmcvid[1] && E(i, 1) == mmcvid[2]) {
                        constraints.ev_constraints.emplace_back(
                            i, -mmcvid[0] - 1);
                        constraints.ev_constraints.back().multiplicity =
                            -mmcvid[3];
                        break;
                    }
                }

            } else { // Is FV.
                for (int i = 0; i < F.rows(); i++) {
                    if (F(i, 0) == mmcvid[1] && F(i, 1) == mmcvid[2]
                        && F(i, 2) == mmcvid[3]) {
                        constraints.fv_constraints.emplace_back(
                            i, -mmcvid[0] - 1);
                        break;
                    }
                }
            }
        }
    }
}
