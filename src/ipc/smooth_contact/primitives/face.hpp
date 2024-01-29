#pragma once

#include <ipc/smooth_contact/distance/mollifier.hpp>
#include "primitive.hpp"

namespace ipc {
    class Face : public Primitive
    {
    public:
        // d is a vector from closest point on the face to the point outside of the face
        Face(const long &fid,
            const Eigen::Ref<const Eigen::Vector3d>& d,
            const Eigen::Ref<const Eigen::Vector3d>& v0,
            const Eigen::Ref<const Eigen::Vector3d>& v1,
            const Eigen::Ref<const Eigen::Vector3d>& v2);
        
        bool is_active() override;
        int n_vertices() const override;
        
        double potential() const override;
        Eigen::VectorXd grad() const override;
        Eigen::MatrixXd hessian() const override;
    private:
        const Eigen::Vector3d _d, _v0, _v1, _v2;
        Eigen::Vector3d normal;
    };

    /// @brief d points from triangle to the point
    template <typename scalar>
    scalar smooth_face_term(
        const Eigen::Ref<const Vector3<scalar>>& v0,
        const Eigen::Ref<const Vector3<scalar>>& v1,
        const Eigen::Ref<const Vector3<scalar>>& v2)
    {
        return 0.5 * (v1 - v0).cross(v2 - v0).norm(); // area of triangle
    }
}