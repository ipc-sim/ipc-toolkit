#include "edge_edge.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/quadrature.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {
    namespace {
        template <class T>
        std::array<VectorMax3<T>, 4> slice_positions(const VectorMax12d &positions, const int &_dim)
        {
            VectorMax3<T> e00, e01, e10, e11;
            e00.setZero(_dim);
            e01.setZero(_dim);
            e10.setZero(_dim);
            e11.setZero(_dim);
            for (int d = 0; d < _dim; d++)
            {
                e00(d) = T(d, positions(d));
                e01(d) = T(_dim + d, positions(_dim + d));
                e10(d) = T(2*_dim + d, positions(2*_dim + d));
                e11(d) = T(3*_dim + d, positions(3*_dim + d));
            }

            return {{e00, e01, e10, e11}};
        }
    }

    Vector12d SmoothEdgeEdgeCollision::positions_to_3d(const VectorMax12d& positions) const
    {
        assert(positions.size() == 8 || positions.size() == 12);
        const int dim = positions.size() / 4;
        Vector12d positions_full;
        positions_full.setZero();
        for (int i = 0; i < 4; i++)
            positions_full.segment(i * 3, dim) = positions.segment(i * dim, dim);
        return positions_full;
    }

    double SmoothEdgeEdgeCollision::compute_distance(const VectorMax12d& positions) const
    {
        return EdgeEdgeCandidate::compute_distance(positions_to_3d(positions));
    }

    VectorMax12d SmoothEdgeEdgeCollision::compute_distance_gradient(
        const VectorMax12d& positions) const
    {
        return EdgeEdgeCandidate::compute_distance_gradient(positions_to_3d(positions));
    }

    MatrixMax12d
    SmoothEdgeEdgeCollision::compute_distance_hessian(const VectorMax12d& positions) const
    {
        return EdgeEdgeCandidate::compute_distance_hessian(positions_to_3d(positions));
    }

    bool
    SmoothEdgeEdgeCollision::ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance,
        const double tmax,
        const double tolerance,
        const long max_iterations,
        const double conservative_rescaling) const
    {
        return EdgeEdgeCandidate::ccd(positions_to_3d(vertices_t0), positions_to_3d(vertices_t1), toi, min_distance, tmax, tolerance, max_iterations, conservative_rescaling);
    }

    double SmoothEdgeEdgeCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const int _dim = positions.size() / num_vertices();
        double val = 0;
        Eigen::VectorXd uv, w;
        line_quadrature(params.n_quadrature, uv, w);
        VectorMax3d p;
        VectorMax3d edge0 = positions.segment(_dim, _dim) - positions.segment(0, _dim);
        for (int i = 0; i < uv.size(); ++i)
        {
            p = positions.segment(0, _dim) + edge0 * uv(i);
            val += w(i) * smooth_point_edge_potential_single_point<double>(p, positions.segment(2*_dim, _dim), positions.segment(3*_dim, _dim), params);
        }
        return val * edge0.norm();

        // tbb::enumerable_thread_specific<VectorMax3d> storage;
        // tbb::enumerable_thread_specific<double> out(std::numeric_limits<double>::infinity());

        // tbb::parallel_for(
        //     tbb::blocked_range<size_t>(size_t(0), uv.size()),
        //     [&](const tbb::blocked_range<size_t>& r) {
        //         double &local_potential = out.local();
        //         local_potential = 0;
        //         for (int i = r.start(); i < r.end(); ++i)
        //         {
        //             storage.local() = positions.segment(0, _dim) + (positions.segment(_dim, _dim) - positions.segment(0, _dim)) * uv(i);
        //             local_potential += w(i) * smooth_point_edge_potential_single_point<double>(p, positions.segment(2*_dim, _dim), positions.segment(3*_dim, _dim), params);
        //         }
        //     });
        
        // return out.combine([](double a, double b) { return a + b; });
    }

    VectorMax12d SmoothEdgeEdgeCollision::gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        const int _dim = positions.size() / num_vertices();
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarGrad<12>;
        auto [e00, e01, e10, e11] = slice_positions<Diff>(positions, _dim);

        VectorMax12d grad;
        grad.setZero(positions.size());
        Eigen::VectorXd uv, w;
        line_quadrature(params.n_quadrature, uv, w);
        VectorMax3<Diff> p;
        VectorMax3<Diff> edge0 = e01 - e00;
        for (int i = 0; i < uv.size(); ++i)
        {
            p = e00 + static_cast<Diff>(uv(i)) * edge0;
            const auto val = edge0.norm() * smooth_point_edge_potential_single_point<Diff>(p, e10, e11, params);
            grad += w(i) * val.getGradient().head(4*_dim);
        }
        return grad;
    }

    MatrixMax12d SmoothEdgeEdgeCollision::hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd) const
    {
        const int _dim = positions.size() / num_vertices();
        DiffScalarBase::setVariableCount(12);
        using Diff=AutodiffScalarHessian<12>;
        auto [e00, e01, e10, e11] = slice_positions<Diff>(positions, _dim);

        MatrixMax12d hess;
        hess.setZero(positions.size(), positions.size());
        Eigen::VectorXd uv, w;
        line_quadrature(params.n_quadrature, uv, w);
        VectorMax3<Diff> p;
        VectorMax3<Diff> edge0 = e01 - e00;
        for (int i = 0; i < uv.size(); ++i)
        {
            p = e00 + static_cast<Diff>(uv(i)) * edge0;
            const auto val = edge0.norm() * smooth_point_edge_potential_single_point<Diff>(p, e10, e11, params);
            hess += w(i) * val.getHessian().topLeftCorner(4*_dim, 4*_dim);
        }

        return hess;
    }
}