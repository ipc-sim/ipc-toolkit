#include <ipc/collision_filter.hpp>

#include <igl/adjacency_matrix.h>
#include <igl/connected_components.h>

namespace ipc {

CollisionFilter
make_connected_components_filter(Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    Eigen::SparseMatrix<int> A;
    igl::adjacency_matrix(faces, A);
    Eigen::VectorXi C, K;
    igl::connected_components(A, C, K);
    return make_vertex_patches_filter(C);
}

} // namespace ipc
