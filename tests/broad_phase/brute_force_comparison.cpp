#include <catch2/catch.hpp>

#include <tbb/parallel_sort.h>

#include <ipc/broad_phase/brute_force.hpp>

void brute_force_comparison(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& group_ids,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    using namespace ipc;

    Candidates bf_candidates;

    auto can_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids.size() == 0 || group_ids(vi) != group_ids(vj);
    };

    detect_collision_candidates_brute_force(
        V0, V1, E, F, bf_candidates,
        /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true,
        /*perform_aabb_check=*/false, inflation_radius, can_collide);

    auto& ev_candidates = candidates.ev_candidates;
    auto& ee_candidates = candidates.ee_candidates;
    auto& fv_candidates = candidates.fv_candidates;
    auto& bf_ev_candidates = bf_candidates.ev_candidates;
    auto& bf_ee_candidates = bf_candidates.ee_candidates;
    auto& bf_fv_candidates = bf_candidates.fv_candidates;

    CHECK(ev_candidates.size() <= bf_ev_candidates.size());
    CHECK(ee_candidates.size() <= bf_ee_candidates.size());
    CHECK(fv_candidates.size() <= bf_fv_candidates.size());

    tbb::parallel_sort(ev_candidates.begin(), ev_candidates.end());
    tbb::parallel_sort(bf_ev_candidates.begin(), bf_ev_candidates.end());
    int ci = 0;
    for (int bf_ci = 0; bf_ci < bf_ev_candidates.size(); bf_ci++) {
        if (ev_candidates.size() <= ci
            || bf_ev_candidates[bf_ci] != ev_candidates[ci]) {
            // Perform CCD to make sure the candidate is not a collision
            double toi;
            bool hit = bf_ev_candidates[bf_ci].ccd(
                V0, V1, E, F, toi, /*tmax=*/1.0, /*tolerance=*/1e-6,
                /*max_iterations=*/1e7, /*conservative_rescaling=*/1.0);
            CHECK(!hit); // Check for FN

        } else {
            ci++;
        }
    }
    CHECK(ci >= ev_candidates.size());

    tbb::parallel_sort(ee_candidates.begin(), ee_candidates.end());
    tbb::parallel_sort(bf_ee_candidates.begin(), bf_ee_candidates.end());
    ci = 0;
    for (int bf_ci = 0; bf_ci < bf_ee_candidates.size(); bf_ci++) {
        if (ee_candidates.size() <= ci
            || bf_ee_candidates[bf_ci] != ee_candidates[ci]) {
            // Perform CCD to make sure the candidate is not a collision
            double toi;
            bool hit = bf_ee_candidates[bf_ci].ccd(
                V0, V1, E, F, toi, /*tmax=*/1.0, /*tolerance=*/1e-6,
                /*max_iterations=*/1e7, /*conservative_rescaling=*/1.0);
            CHECK(!hit); // Check for FN

        } else {
            ci++;
        }
    }
    CHECK(ci >= ee_candidates.size());

    tbb::parallel_sort(fv_candidates.begin(), fv_candidates.end());
    tbb::parallel_sort(bf_fv_candidates.begin(), bf_fv_candidates.end());
    ci = 0;
    for (int bf_ci = 0; bf_ci < bf_fv_candidates.size(); bf_ci++) {
        if (fv_candidates.size() <= ci
            || bf_fv_candidates[bf_ci] != fv_candidates[ci]) {
            // Perform CCD to make sure the candidate is not a collision
            double toi;
            bool hit = bf_fv_candidates[bf_ci].ccd(
                V0, V1, E, F, toi, /*tmax=*/1.0, /*tolerance=*/1e-6,
                /*max_iterations=*/1e7, /*conservative_rescaling=*/1.0);
            CHECK(!hit); // Check for FN
        } else {
            ci++;
        }
    }
    CHECK(ci >= fv_candidates.size());
}
