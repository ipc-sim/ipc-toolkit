#include <catch2/catch_test_macros.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/candidates/edge_face.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/io/write_candidates_obj.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>

#include <filesystem>
#include <fstream>
#include <sstream>

TEST_CASE("Logger", "[utils][logger]")
{
    const std::shared_ptr<spdlog::logger> custom_logger =
        std::make_shared<spdlog::logger>(
            "custom", std::make_shared<spdlog::sinks::stdout_color_sink_mt>());

    CHECK(&ipc::logger() != custom_logger.get());
    ipc::set_logger(custom_logger);
    CHECK(&ipc::logger() == custom_logger.get());
}

TEST_CASE("Project to PSD", "[utils][project_to_psd]")
{
    Eigen::MatrixXd A, A_psd;

    A.setIdentity(3, 3);
    A_psd = ipc::project_to_psd(A, ipc::PSDProjectionMethod::CLAMP);
    CHECK(A_psd.isApprox(A));
    A_psd = ipc::project_to_psd(A, ipc::PSDProjectionMethod::ABS);
    CHECK(A_psd.isApprox(A));

    A *= -1;
    A_psd = ipc::project_to_psd(A, ipc::PSDProjectionMethod::CLAMP);
    CHECK(A_psd.isZero());

    A.resize(2, 2);
    A.row(0) << 2, 1;
    A.row(1) << 1, 2;
    A_psd = ipc::project_to_psd(A, ipc::PSDProjectionMethod::CLAMP);
    CHECK(A_psd.isApprox(A));
    A_psd = ipc::project_to_psd(A, ipc::PSDProjectionMethod::ABS);
    CHECK(A_psd.isApprox(A));
}

TEST_CASE("Project to PD", "[utils][project_to_pd]")
{
    Eigen::MatrixXd A, A_pd;

    A.setIdentity(3, 3);
    A_pd = ipc::project_to_pd(A);
    CHECK(A_pd.isApprox(A));

    A *= -1;
    A_pd = ipc::project_to_pd(A);
    CHECK(A_pd.isApprox(1e-8 * Eigen::MatrixXd::Identity(3, 3)));

    A.resize(2, 2);
    A.row(0) << 2, 1;
    A.row(1) << 1, 2;
    A_pd = ipc::project_to_pd(A);
    CHECK(A_pd.isApprox(A));
}

TEST_CASE("Write OBJ of candidates", "[utils][write_candidates_obj]")
{
    Eigen::MatrixXd V(4, 3);
    V.row(0) << 0, 0, 0;
    V.row(1) << 1, 0, 0;
    V.row(2) << 0, 1, 0;
    V.row(3) << 0, 0, 1;
    Eigen::MatrixXi E(2, 2);
    E.row(0) << 1, 2;
    E.row(1) << 0, 3;
    Eigen::MatrixXi F(1, 3);
    F.row(0) << 1, 2, 3;
    SECTION("VertexVertexCandidate")
    {
        std::stringstream ss;
        ipc::write_candidates_obj<ipc::VertexVertexCandidate>(
            ss, V, E, F, { { ipc::VertexVertexCandidate(0, 1) } });
        CHECK(ss.str() == "o VV\nv 0 0 0\nv 1 0 0\n");
    }
    SECTION("EdgeVertexCandidate")
    {
        std::stringstream ss;
        ipc::write_candidates_obj<ipc::EdgeVertexCandidate>(
            ss, V, E, F, { { ipc::EdgeVertexCandidate(0, 0) } });
        CHECK(ss.str() == "o EV\nv 1 0 0\nv 0 1 0\nv 0 0 0\nl 1 2\n");
    }
    SECTION("EdgeEdgeCandidate")
    {
        std::stringstream ss;
        ipc::write_candidates_obj<ipc::EdgeEdgeCandidate>(
            ss, V, E, F, { { ipc::EdgeEdgeCandidate(0, 1) } });
        CHECK(
            ss.str()
            == "o EE\nv 1 0 0\nv 0 1 0\nv 0 0 0\nv 0 0 1\nl 1 2\nl 3 4\n");
    }
    SECTION("FaceVertexCandidate")
    {
        std::stringstream ss;
        ipc::write_candidates_obj<ipc::FaceVertexCandidate>(
            ss, V, E, F, { { ipc::FaceVertexCandidate(0, 0) } });
        CHECK(
            ss.str() == "o FV\nv 1 0 0\nv 0 1 0\nv 0 0 1\nv 0 0 0\nf 1 2 3\n");
    }
    SECTION("EdgeFaceCandidate")
    {
        std::stringstream ss;
        ipc::write_candidates_obj<ipc::EdgeFaceCandidate>(
            ss, V, E, F, { { ipc::EdgeFaceCandidate(0, 0) } });
        CHECK(
            ss.str()
            == "o EF\nv 1 0 0\nv 0 1 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nl 1 2\nf 3 4 5\n");
    }
}
TEST_CASE("Write OBJ of candidates to a file", "[utils][write_candidates_obj]")
{
    Eigen::MatrixXd V(3, 3);
    V.row(0) << 0, 0, 0;
    V.row(1) << 1, 0, 0;
    V.row(2) << 0, 1, 0;
    Eigen::MatrixXi E(1, 2);
    E.row(0) << 0, 1;
    const Eigen::MatrixXi F(0, 3);

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "ipctk_candidates_test.obj";
    const auto read_file = [&path]() {
        std::ifstream in(path);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    };

    ipc::Candidates candidates;
    candidates.vv_candidates.emplace_back(0, 1);

    SECTION("Candidates::write_obj")
    {
        CHECK(candidates.write_obj(path.string(), V, E, F));
        // Groups with no candidates still get an (empty) object header.
        CHECK(read_file() == "o VV\nv 0 0 0\nv 1 0 0\no EV\no EE\no FV\n");
    }
    SECTION("write_candidates_obj with a filename")
    {
        CHECK(
            ipc::write_candidates_obj<ipc::VertexVertexCandidate>(
                path.string(), V, E, F, candidates.vv_candidates));
        CHECK(read_file() == "o VV\nv 0 0 0\nv 1 0 0\n");
    }
    SECTION("An unwritable path reports failure")
    {
        const std::string bad = (std::filesystem::temp_directory_path()
                                 / "ipctk_no_such_dir" / "candidates.obj")
                                    .string();
        CHECK_FALSE(candidates.write_obj(bad, V, E, F));
        CHECK_FALSE(
            ipc::write_candidates_obj<ipc::VertexVertexCandidate>(
                bad, V, E, F, candidates.vv_candidates));
    }

    std::filesystem::remove(path);
}
