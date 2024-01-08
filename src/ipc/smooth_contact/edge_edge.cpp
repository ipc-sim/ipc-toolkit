#include "edge_edge.hpp"
#include "smooth_point_edge.hpp"
#include <ipc/utils/quadrature.hpp>
#include <ipc/utils/AutodiffTypes.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <iostream>


namespace ipc {
    namespace {
        template <class T>
        std::array<VectorMax3<T>, 4> slice_positions(const VectorMax12d &positions, const int &_dim)
        {
            std::array<VectorMax3<T>, 4> points;
            points.fill(VectorMax3<T>::Zero(_dim));
            
            for (int i = 0, id = 0; i < 4; i++)
                for (int d = 0; d < _dim; d++, id++)
                    if constexpr (std::is_same<T, double>::value)
                        points[i](d) = positions(id);
                    else
                        points[i](d) = T(id, positions(id));

            return points;
        }
    }

    Vector12d SmoothEdgeEdgeCollision::positions_to_3d(const VectorMax12d& positions) const
    {
        if (positions.size() == 12)
            return positions;
        
        assert(positions.size() == 8);
        Vector12d positions_full;
        positions_full.setZero();
        for (int i = 0; i < 4; i++)
            positions_full.segment<2>(i * 3) = positions.segment<2>(i * 2);
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

    double SmoothEdgeEdgeCollision::operator()(
        const VectorMax12d& positions, 
        const ParameterType &params) const
    {
        std::cout << "smooth edge edge: " << this->edge0_id << this->edge1_id << std::endl;
        const int _dim = positions.size() / num_vertices();
        auto [e00, e01, e10, e11] = slice_positions<double>(positions, _dim);

        double val = 0;
        Eigen::VectorXd uv, w;
        line_quadrature(params.n_quadrature, uv, w);
        VectorMax3d p;
        const VectorMax3d edge0 = e01 - e00;
        const VectorMax3d edge1 = e11 - e10;
        const double len0 = edge0.norm(), len1 = edge1.norm();
        for (int i = 0; i < uv.size(); ++i)
        {
            p = e00 + edge0 * uv(i);
            val += len0 * w(i) * smooth_point_edge_potential_quadrature<double>(p, e10, e11, params);

            p = e10 + edge1 * uv(i);
            val += len1 * w(i) * smooth_point_edge_potential_quadrature<double>(p, e00, e01, params);
        }
        return val;

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

        Eigen::VectorXd uv, w;
        line_quadrature(params.n_quadrature, uv, w);
        VectorMax3<Diff> p;
        const VectorMax3<Diff> edge0 = e01 - e00;
        const VectorMax3<Diff> edge1 = e11 - e10;
        const Diff len0 = edge0.norm(), len1 = edge1.norm();
        Diff val = Diff(0.);
        for (int i = 0; i < uv.size(); ++i)
        {
            p = e00 + Diff(uv(i)) * edge0;
            val += Diff(w(i)) * len0 * smooth_point_edge_potential_quadrature<Diff>(p, e10, e11, params);

            p = e10 + Diff(uv(i)) * edge1;
            val += Diff(w(i)) * len1 * smooth_point_edge_potential_quadrature<Diff>(p, e00, e01, params);
        }
        return val.getGradient().head(4*_dim);
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

        Eigen::VectorXd uv, w;
        line_quadrature(params.n_quadrature, uv, w);
        VectorMax3<Diff> p;
        const VectorMax3<Diff> edge0 = e01 - e00;
        const VectorMax3<Diff> edge1 = e11 - e10;
        const Diff len0 = edge0.norm(), len1 = edge1.norm();
        Diff val = Diff(0.);
        for (int i = 0; i < uv.size(); ++i)
        {
            p = e00 + Diff(uv(i)) * edge0;
            val += Diff(w(i)) * len0 * smooth_point_edge_potential_quadrature<Diff>(p, e10, e11, params);

            p = e10 + Diff(uv(i)) * edge1;
            val += Diff(w(i)) * len1 * smooth_point_edge_potential_quadrature<Diff>(p, e00, e01, params);
        }

        return val.getHessian().topLeftCorner(4*_dim, 4*_dim);
    }
}