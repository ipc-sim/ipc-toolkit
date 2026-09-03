#include "point_edge.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>

#include <stdexcept> // std::invalid_argument

namespace ipc {

namespace detail {

    template <typename T, int dim>
    Eigen::Matrix<T, 3 * dim, 3 * dim> point_edge_distance_hessian(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1,
        PointEdgeDistanceType dtype)
    {
        static_assert(dim == 2 || dim == 3, "point-edge is only 2D or 3D");

        if constexpr (std::is_floating_point_v<T>) {
            if (dtype == PointEdgeDistanceType::AUTO) {
                dtype = point_edge_distance_type<T, dim>(p, e0, e1);
            }
        } else if (dtype == PointEdgeDistanceType::AUTO) {
            throw_auto_requires_explicit_dtype("point_edge_distance_hessian");
        }

        Eigen::Matrix<T, 3 * dim, 3 * dim> hess =
            Eigen::Matrix<T, 3 * dim, 3 * dim>::Zero();

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            hess.template topLeftCorner<2 * dim, 2 * dim>() =
                point_point_distance_hessian<T, dim>(p, e0);
            break;

        case PointEdgeDistanceType::P_E1: {
            const Eigen::Matrix<T, 2 * dim, 2 * dim> local_hess =
                point_point_distance_hessian<T, dim>(p, e1);
            hess.template topLeftCorner<dim, dim>() =
                local_hess.template topLeftCorner<dim, dim>();
            hess.template topRightCorner<dim, dim>() =
                local_hess.template topRightCorner<dim, dim>();
            hess.template bottomLeftCorner<dim, dim>() =
                local_hess.template bottomLeftCorner<dim, dim>();
            hess.template bottomRightCorner<dim, dim>() =
                local_hess.template bottomRightCorner<dim, dim>();
            break;
        }

        case PointEdgeDistanceType::P_E:
            hess = point_line_distance_hessian<T, dim>(p, e0, e1);
            break;

        default:
            throw_invalid_distance_type("point_edge_distance_hessian");
        }
        return hess;
    }

#define IPC_INSTANTIATE_POINT_EDGE_DISTANCE_HESSIAN(T, dim)                    \
    template Eigen::Matrix<T, 3 * dim, 3 * dim>                                \
    point_edge_distance_hessian<T, dim>(                                       \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>,                                \
        Eigen::ConstRef<Eigen::Vector<T, dim>>, PointEdgeDistanceType)

    IPC_INSTANTIATE_POINT_EDGE_DISTANCE_HESSIAN(float, 2);
    IPC_INSTANTIATE_POINT_EDGE_DISTANCE_HESSIAN(float, 3);
    IPC_INSTANTIATE_POINT_EDGE_DISTANCE_HESSIAN(double, 2);
    IPC_INSTANTIATE_POINT_EDGE_DISTANCE_HESSIAN(double, 3);

#undef IPC_INSTANTIATE_POINT_EDGE_DISTANCE_HESSIAN

} // namespace detail

} // namespace ipc
