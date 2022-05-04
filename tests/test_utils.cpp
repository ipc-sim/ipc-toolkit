#include "test_utils.hpp"

#include <igl/dirname.h>
#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <ipc/utils/eigen_ext.hpp>

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
                constraints.vv_constraints.back().weight = -mmcvid[3];

            } else if (mmcvid[3] < 0) { // Is EV?
                int ei;
                for (ei = 0; ei < E.rows(); ei++) {
                    if (E(ei, 0) == mmcvid[1] && E(ei, 1) == mmcvid[2]) {
                        break;
                    }
                }
                assert(ei < E.rows());
                constraints.ev_constraints.emplace_back(ei, -mmcvid[0] - 1);
                constraints.ev_constraints.back().weight = -mmcvid[3];

            } else { // Is FV.
                int fi;
                for (fi = 0; fi < F.rows(); fi++) {
                    if (F(fi, 0) == mmcvid[1] && F(fi, 1) == mmcvid[2]
                        && F(fi, 2) == mmcvid[3]) {
                        break;
                    }
                }
                assert(fi < F.rows());
                constraints.fv_constraints.emplace_back(fi, -mmcvid[0] - 1);
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
        std::make_unique<RotationGenerator>());
}
