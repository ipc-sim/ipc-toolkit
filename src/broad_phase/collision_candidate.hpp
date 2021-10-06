#pragma once

#include <vector>
#include <fstream>

#include <Eigen/Core>

#include <ipc/ccd/ccd.hpp>

namespace ipc {

struct CollisionCandidate {
    virtual ~CollisionCandidate() {}

    virtual bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double& toi,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const = 0;
};

struct VertexVertexCandidate {
    VertexVertexCandidate(long vertex0_index, long vertex1_index);

    bool operator==(const VertexVertexCandidate& other) const;
    bool operator!=(const VertexVertexCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const VertexVertexCandidate& other) const;

    long vertex0_index;
    long vertex1_index;
};

struct EdgeVertexCandidate : CollisionCandidate {
    EdgeVertexCandidate(long edge_index, long vertex_index);

    bool operator==(const EdgeVertexCandidate& other) const;
    bool operator!=(const EdgeVertexCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const EdgeVertexCandidate& other) const;

    bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double& toi,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    long edge_index;
    long vertex_index;
};

struct EdgeEdgeCandidate : CollisionCandidate {
    EdgeEdgeCandidate(long edge0_index, long edge1_index);

    bool operator==(const EdgeEdgeCandidate& other) const;
    bool operator!=(const EdgeEdgeCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeEdgeCandidates for sorting.
    bool operator<(const EdgeEdgeCandidate& other) const;

    bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double& toi,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    long edge0_index;
    long edge1_index;
};

/// @brief Candidate for <b>intersection</b> between edge and face.
///
/// Not included in Candidates because it is not a collision candidate.
struct EdgeFaceCandidate {
    EdgeFaceCandidate(long edge_index, long face_index);

    bool operator==(const EdgeFaceCandidate& other) const;
    bool operator!=(const EdgeFaceCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeFaceCandidate for sorting.
    bool operator<(const EdgeFaceCandidate& other) const;

    long edge_index;
    long face_index;
};

struct FaceVertexCandidate : CollisionCandidate {
    FaceVertexCandidate(long face_index, long vertex_index);

    bool operator==(const FaceVertexCandidate& other) const;
    bool operator!=(const FaceVertexCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare FaceVertexCandidate for sorting.
    bool operator<(const FaceVertexCandidate& other) const;

    bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double& toi,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    long face_index;
    long vertex_index;
};

struct Candidates {
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;

    size_t size() const;

    bool empty() const;

    void clear();

    CollisionCandidate& operator[](size_t idx);
    const CollisionCandidate& operator[](size_t idx) const;

    bool save_obj(
        const std::string& filename,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F);
};

////////////////////////////////////////////////////////////////////////////////
// Debugging functions

void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeVertexCandidate>& ev_candidates);
void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeEdgeCandidate>& ee_candidates);
void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<FaceVertexCandidate>& fv_candidates);
void save_obj(
    std::ofstream& out,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<EdgeFaceCandidate>& ef_candidates);

template <typename Candidate>
bool save_obj(
    const std::string& filename,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::vector<Candidate>& candidates)
{
    std::ofstream obj(filename);
    if (!obj.is_open()) {
        return false;
    }
    save_obj(obj, V, E, F, candidates);
    return true;
}

} // namespace ipc
