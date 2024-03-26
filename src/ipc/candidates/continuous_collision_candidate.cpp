#include "continuous_collision_candidate.hpp"

namespace ipc {

std::ostream& ContinuousCollisionCandidate::write_ccd_query(
    std::ostream& out,
    const VectorMax12d& vertices_t0,
    const VectorMax12d& vertices_t1) const
{
    assert(vertices_t0.size() == vertices_t1.size());

    const int dim = vertices_t0.size() / num_vertices();
    assert(vertices_t0.size() % num_vertices() == 0);

    for (int i = 0; i < num_vertices(); i++) {
        out << vertices_t0.segment(dim * i, dim)
                   .transpose()
                   .format(OBJ_VERTEX_FORMAT);
    }

    for (int i = 0; i < num_vertices(); i++) {
        out << vertices_t1.segment(dim * i, dim)
                   .transpose()
                   .format(OBJ_VERTEX_FORMAT);
    }

    return out;
}

} // namespace ipc
