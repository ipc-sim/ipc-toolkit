#include "point2.hpp"

#include <ipc/utils/logger.hpp>

namespace ipc {

namespace {
    std::vector<index_t>
    find_vertex_adjacent_vertices(const CollisionMesh& mesh, const index_t v)
    {
        assert(mesh.dim() == 2);
        std::vector<index_t> neighbors;
        neighbors.assign(2, -1);

        for (index_t i : mesh.vertex_edge_adjacencies()[v]) {
            if (mesh.edges()(i, 0) == v) {
                neighbors[0] = mesh.edges()(i, 1);
            } else if (mesh.edges()(i, 1) == v) {
                neighbors[1] = mesh.edges()(i, 0);
            } else {
                log_and_throw_error("Invalid edge-vertex adjacency!");
            }
        }

        return neighbors;
    }

    bool smooth_point2_term_type(
        Eigen::ConstRef<Eigen::Vector2d> v,
        Eigen::ConstRef<Eigen::Vector2d> direc,
        Eigen::ConstRef<Eigen::Vector2d> e0,
        Eigen::ConstRef<Eigen::Vector2d> e1,
        const GCPParameters& params,
        const bool orientable)
    {
        const Eigen::Vector2d dn = -direc.normalized();
        const Eigen::Vector2d t0 = (e0 - v).normalized(),
                              t1 = (e1 - v).normalized();

        if (dn.dot(t0) <= -params.alpha_t || dn.dot(t1) <= -params.alpha_t) {
            return false;
        }

        if (orientable) {
            const double tmp = Math<double>::smooth_heaviside(
                                   -Math<double>::cross2(dn, t0),
                                   params.alpha_n, params.beta_n)
                + Math<double>::smooth_heaviside(
                                   Math<double>::cross2(dn, t1), params.alpha_n,
                                   params.beta_n);
            if (tmp <= 1. - params.alpha_n) {
                return false;
            }
        }

        return true;
    }

    template <class scalar>
    scalar smooth_point2_term(
        Eigen::ConstRef<Eigen::Vector2<scalar>> v,
        Eigen::ConstRef<Eigen::Vector2<scalar>> direc,
        Eigen::ConstRef<Eigen::Vector2<scalar>> e0,
        Eigen::ConstRef<Eigen::Vector2<scalar>> e1,
        const GCPParameters& params,
        const bool orientable,
        const scalar measure)
    {
        const Eigen::Vector2<scalar> dn = -direc.normalized();
        const Eigen::Vector2<scalar> t0 = (e0 - v).normalized(),
                                     t1 = (v - e1).normalized();

        scalar val = Math<scalar>::smooth_heaviside(
                         dn.dot(t0), params.alpha_t, params.beta_t)
            * Math<scalar>::smooth_heaviside(
                         -dn.dot(t1), params.alpha_t, params.beta_t);

        if (orientable) {
            const scalar tmp = Math<scalar>::smooth_heaviside(
                                   -Math<scalar>::cross2(dn, t0),
                                   params.alpha_n, params.beta_n)
                + Math<scalar>::smooth_heaviside(
                                   -Math<scalar>::cross2(dn, t1),
                                   params.alpha_n, params.beta_n);
            val = val
                * Math<scalar>::smooth_heaviside(tmp - 1., params.alpha_n, 0);
        }

        return val * measure;
    }

    template <class scalar>
    scalar smooth_point2_term_one_side(
        Eigen::ConstRef<Eigen::Vector2<scalar>> v,
        Eigen::ConstRef<Eigen::Vector2<scalar>> direc,
        Eigen::ConstRef<Eigen::Vector2<scalar>> e0,
        const GCPParameters& params,
        const scalar measure)
    {
        const Eigen::Vector2<scalar> dn = -direc.normalized();
        const Eigen::Vector2<scalar> t0 = e0 - v;

        const scalar tangent_term = Math<scalar>::smooth_heaviside(
            dn.dot(t0) / t0.norm(), params.alpha_t, params.beta_t);

        return tangent_term * measure;
    }
    // We select the quadrature weight in one place so the potential, gradient,
    // and Hessian can never disagree about which measure they differentiate.
    // With the rest-shape measure the weight is a constant, so no derivative
    // flows through it. Otherwise it is the (deformed) average half-edge
    // length, which does depend on the vertex positions.
    template <class scalar>
    scalar point2_two_sided_measure(
        const bool use_rest_shape_measure,
        const double rest_measure,
        Eigen::ConstRef<Eigen::Vector2<scalar>> v,
        Eigen::ConstRef<Eigen::Vector2<scalar>> e0,
        Eigen::ConstRef<Eigen::Vector2<scalar>> e1)
    {
        if (use_rest_shape_measure) {
            return scalar(rest_measure);
        }
        return ((e0 - v).norm() + (e1 - v).norm()) / 2.0;
    }

    // Same as above for a vertex with a single neighboring edge, where the
    // weight is the length of that edge.
    template <class scalar>
    scalar point2_one_sided_measure(
        const bool use_rest_shape_measure,
        const double rest_measure,
        Eigen::ConstRef<Eigen::Vector2<scalar>> v,
        Eigen::ConstRef<Eigen::Vector2<scalar>> e0)
    {
        if (use_rest_shape_measure) {
            return scalar(rest_measure);
        }
        return (e0 - v).norm();
    }
} // namespace

Point2::Point2(
    const index_t id,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<VectorMax3d> d,
    const GCPParameters& params)
    : Primitive(id, params)
{
    orientable = mesh.is_orient_vertex(id);
    auto neighbor_verts = find_vertex_adjacent_vertices(mesh, id);
    has_neighbor_1 = neighbor_verts[0] >= 0;
    has_neighbor_2 = neighbor_verts[1] >= 0;

    if (has_neighbor_1 && has_neighbor_2) {
        m_vertex_ids = { { id, neighbor_verts[0], neighbor_verts[1] } };
        m_is_active = smooth_point2_term_type(
            vertices.row(id), d, vertices.row(m_vertex_ids[1]),
            vertices.row(m_vertex_ids[2]), params, orientable);

        // rest-shape measure: average half-edge length from rest positions
        const Eigen::MatrixXd& rp = mesh.rest_positions();
        m_rest_measure = ((rp.row(neighbor_verts[0]) - rp.row(id)).norm()
                          + (rp.row(neighbor_verts[1]) - rp.row(id)).norm())
            / 2.0;
    } else if (has_neighbor_1 || has_neighbor_2) {
        m_vertex_ids = {
            { id, has_neighbor_1 ? neighbor_verts[0] : neighbor_verts[1] }
        };

        const Eigen::Vector2d dn = -d.normalized();
        const Eigen::Vector2d t0 =
            (vertices.row(m_vertex_ids[1]) - vertices.row(m_vertex_ids[0]))
                .normalized();

        m_is_active = dn.dot(t0) > -params.alpha_t;
        // rest-shape measure: half-edge length from rest positions
        const Eigen::MatrixXd& rp = mesh.rest_positions();
        m_rest_measure =
            (rp.row(m_vertex_ids[1]) - rp.row(m_vertex_ids[0])).norm();
    } else {
        m_vertex_ids.resize(1);
        m_vertex_ids[0] = id;
        m_is_active = true;
        m_rest_measure = 1.0; // no neighbors — unit weight
    }
}

int Point2::n_vertices() const { return m_vertex_ids.size(); }

double Point2::potential(
    const Eigen::Vector<double, DIM>& d,
    const VectorMax<double, MAX_SIZE>& x) const
{
    if (has_neighbor_1 && has_neighbor_2) {
        const double measure = point2_two_sided_measure<double>(
            m_params.use_rest_shape_measure, m_rest_measure, x.segment<DIM>(0),
            x.segment<DIM>(DIM), x.segment<DIM>(2 * DIM));
        return smooth_point2_term<double>(
            x.segment<DIM>(0), d, x.segment<DIM>(DIM), x.segment<DIM>(2 * DIM),
            m_params, orientable, measure);
    } else if (has_neighbor_1 || has_neighbor_2) {
        const double measure = point2_one_sided_measure<double>(
            m_params.use_rest_shape_measure, m_rest_measure, x.segment<DIM>(0),
            x.segment<DIM>(DIM));
        return smooth_point2_term_one_side<double>(
            x.segment<DIM>(0), d, x.segment<DIM>(DIM), m_params, measure);
    } else {
        return m_rest_measure;
    }
}
VectorMax<double, Point2::MAX_SIZE + Point2::DIM> Point2::grad(
    const Eigen::Vector<double, DIM>& d,
    const VectorMax<double, MAX_SIZE>& x) const
{
    if (has_neighbor_1 && has_neighbor_2) {
        ScalarBase::setVariableCount(4 * DIM);
        using T = ADGrad<4 * DIM>;
        Eigen::Vector<double, 4 * DIM> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 4, DIM> X = slice_positions<T, 4, DIM>(tmp);
        const T measure = point2_two_sided_measure<T>(
            m_params.use_rest_shape_measure, m_rest_measure, X.row(1), X.row(2),
            X.row(3));
        return smooth_point2_term<T>(
                   X.row(1), X.row(0), X.row(2), X.row(3), m_params, orientable,
                   measure)
            .grad;
    } else if (has_neighbor_1 || has_neighbor_2) {
        ScalarBase::setVariableCount(3 * DIM);
        using T = ADGrad<3 * DIM>;
        Eigen::Vector<double, 3 * DIM> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 3, DIM> X = slice_positions<T, 3, DIM>(tmp);
        const T measure = point2_one_sided_measure<T>(
            m_params.use_rest_shape_measure, m_rest_measure, X.row(1),
            X.row(2));
        return smooth_point2_term_one_side<T>(
                   X.row(1), X.row(0), X.row(2), m_params, measure)
            .grad;
    } else {
        return VectorMax<double, Point2::MAX_SIZE + Point2::DIM>::Zero(
            x.size() + d.size());
    }
}
MatrixMax<
    double,
    Point2::MAX_SIZE + Point2::DIM,
    Point2::MAX_SIZE + Point2::DIM>
Point2::hessian(
    const Eigen::Vector<double, DIM>& d,
    const VectorMax<double, MAX_SIZE>& x) const
{
    if (has_neighbor_1 && has_neighbor_2) {
        ScalarBase::setVariableCount(4 * DIM);
        using T = ADHessian<4 * DIM>;
        Eigen::Vector<double, 4 * DIM> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 4, DIM> X = slice_positions<T, 4, DIM>(tmp);
        const T measure = point2_two_sided_measure<T>(
            m_params.use_rest_shape_measure, m_rest_measure, X.row(1), X.row(2),
            X.row(3));
        return smooth_point2_term<T>(
                   X.row(1), X.row(0), X.row(2), X.row(3), m_params, orientable,
                   measure)
            .Hess;
    } else if (has_neighbor_1 || has_neighbor_2) {
        ScalarBase::setVariableCount(3 * DIM);
        using T = ADHessian<3 * DIM>;
        Eigen::Vector<double, 3 * DIM> tmp;
        tmp << d, x;
        Eigen::Matrix<T, 3, DIM> X = slice_positions<T, 3, DIM>(tmp);
        const T measure = point2_one_sided_measure<T>(
            m_params.use_rest_shape_measure, m_rest_measure, X.row(1),
            X.row(2));
        return smooth_point2_term_one_side<T>(
                   X.row(1), X.row(0), X.row(2), m_params, measure)
            .Hess;
    } else {
        return MatrixMax<double, -1, Point2::MAX_SIZE + Point2::DIM>::Zero(
            x.size() + d.size(), x.size() + d.size());
    }
}
} // namespace ipc