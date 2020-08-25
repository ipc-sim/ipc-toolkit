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

TEST_CASE("Flatten and unflatten", "[utils]")
{
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(1000, 3);
    Eigen::MatrixXd R = unflatten(flatten(X), X.cols());
    CHECK(X == R);
}
