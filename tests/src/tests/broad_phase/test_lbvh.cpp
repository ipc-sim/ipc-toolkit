#include <ipc/config.hpp>

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <ipc/broad_phase/bvh.hpp>
#include <ipc/broad_phase/lbvh.hpp>
// #include <ipc/broad_phase/vulkan/shaders/single_radixsort.hpp>
// #include <ipc/utils/profiler.hpp>

#include <catch2/catch_test_macros.hpp>

#include <iostream>

#include <random>

using namespace ipc;

namespace {

bool is_aabb_union(LBVH::Node parent, LBVH::Node childA, LBVH::Node childB)
{
    AABB children;
    children.min = childA.aabb_min.min(childB.aabb_min);
    children.max = childA.aabb_max.max(childB.aabb_max);
    constexpr float EPS = 1e-4f;
    return (abs(parent.aabb_max - children.max) < EPS).all()
        && (abs(parent.aabb_min - children.min) < EPS).all();
}

void traverse_lbvh(
    const std::vector<LBVH::Node>& lbvh_nodes,
    const uint32_t index,
    std::vector<bool>& visited)
{
    const LBVH::Node& node = lbvh_nodes[index];
    CHECK(node.is_valid());

    if (node.left == LBVH::Node::INVALID_POINTER) {
        // leaf
        CHECK(!visited[index]);
        visited[index] = true;
    } else {
        // inner node
        CHECK(!visited[index]);
        visited[index] = true;

        uint32_t left_child_index = LBVH::Node::POINTER(index, node.left);
        uint32_t right_child_index = LBVH::Node::POINTER(index, node.right);

        // verify aabbs
        LBVH::Node childA = lbvh_nodes[left_child_index];
        LBVH::Node childB = lbvh_nodes[right_child_index];

        {
            CAPTURE(
                index, left_child_index, right_child_index,
                node.aabb_min.transpose(), childA.aabb_min.transpose(),
                childB.aabb_min.transpose(), node.aabb_max.transpose(),
                childA.aabb_max.transpose(), childB.aabb_max.transpose());
            CHECK(is_aabb_union(node, childA, childB));
        }

        // continue traversal
        traverse_lbvh(lbvh_nodes, left_child_index, visited);
        traverse_lbvh(lbvh_nodes, right_child_index, visited);
    }
}

void check_valid_lbvh_nodes(const std::vector<LBVH::Node>& lbvh_nodes)
{
    std::vector<bool> visited(lbvh_nodes.size(), false);
    traverse_lbvh(lbvh_nodes, 0, visited);
    REQUIRE(
        std::all_of(visited.begin(), visited.end(), [](bool v) { return v; }));
}
} // namespace

TEST_CASE("LBVH::build", "[broad_phase][lbvh]")
{
    constexpr double inflation_radius = 1e-3;

    const std::shared_ptr<LBVH> lbvh = std::make_shared<LBVH>();

    SECTION("Static")
    {
        const std::string mesh = GENERATE("cube.ply", "bunny.ply");

        Eigen::MatrixXd vertices;
        Eigen::MatrixXi edges, faces;
        REQUIRE(tests::load_mesh(mesh, vertices, edges, faces));

        lbvh->build(vertices, edges, faces, inflation_radius);
    }
    SECTION("Dynamic")
    {
        const std::string mesh_t0 = "cloth_ball92.ply";
        const std::string mesh_t1 = "cloth_ball93.ply";

        Eigen::MatrixXd vertices_t0, vertices_t1;
        Eigen::MatrixXi edges, faces;
        REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
        REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

        lbvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
    }

    std::ofstream dot_file("/Users/zachary/Downloads/lbvh.dot");
    dot_file << "digraph G {\nnode [shape = circle;];" << std::endl;

    for (int i = 0; i < lbvh->vertex_nodes().size(); ++i) {
        dot_file << "N" << i << " [label = \""
                 << (lbvh->vertex_nodes()[i].is_inner() ? "I" : "L")
                 << (lbvh->vertex_nodes()[i].is_inner()
                         ? i
                         : lbvh->vertex_nodes()[i].primitive_id)
                 << "\"];" << std::endl;
    }

    for (int i = 0; i < lbvh->vertex_nodes().size(); ++i) {
        const auto& node = lbvh->vertex_nodes()[i];
        if (node.is_inner()) {
            dot_file << "N" << i << " -> N" << node.left << ";" << std::endl;
            dot_file << "N" << i << " -> N" << node.right << ";" << std::endl;
        }
    }

    dot_file << "}" << std::endl;
    dot_file.close();

    // -- TODO: Check the morton codes ----------------------------------------
    // -- TODO: Check the morton codes are sorted -----------------------------

    // -- Check the LBVH nodes are all reachable and contain their children ---

    check_valid_lbvh_nodes(lbvh->vertex_nodes());
    // check_valid_lbvh_nodes(lbvh->edge_nodes());
    // check_valid_lbvh_nodes(lbvh->face_nodes());

    // -- Check clear() works -------------------------------------------------
    lbvh->clear();

    // CHECK(lbvh->vertex_boxes().empty());
    // CHECK(lbvh->edge_boxes().empty());
    // CHECK(lbvh->face_boxes().empty());

    CHECK(lbvh->vertex_nodes().empty());
    CHECK(lbvh->edge_nodes().empty());
    CHECK(lbvh->face_nodes().empty());
}

TEST_CASE("LBVH::detect_*_candidates", "[broad_phase][lbvh]")
{
    constexpr double inflation_radius = 0;

    std::string mesh_t0, mesh_t1;
    SECTION("Two cubes")
    {
        mesh_t0 = "two-cubes-far.ply";
        mesh_t1 = "two-cubes-intersecting.ply";
    }
    SECTION("Cloth-Ball")
    {
        mesh_t0 = "cloth_ball92.ply";
        mesh_t1 = "cloth_ball93.ply";
    }

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
    REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

    const std::shared_ptr<LBVH> lbvh = std::make_shared<LBVH>();
    lbvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    const std::shared_ptr<BVH> bvh = std::make_shared<BVH>();
    bvh->build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    // detect_vertex_vertex_candidates
    {
        std::vector<VertexVertexCandidate> vv_candidates;
        lbvh->detect_vertex_vertex_candidates(vv_candidates);

        std::vector<VertexVertexCandidate> expected_vv_candidates;
        bvh->detect_vertex_vertex_candidates(expected_vv_candidates);

        CHECK(vv_candidates.size() == expected_vv_candidates.size());
        // TODO: Check the candidates are the same
    }

    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        lbvh->detect_edge_vertex_candidates(ev_candidates);

        std::vector<EdgeVertexCandidate> expected_ev_candidates;
        bvh->detect_edge_vertex_candidates(expected_ev_candidates);

        CHECK(ev_candidates.size() == expected_ev_candidates.size());
        // TODO: Check the candidates are the same
    }

    {
        std::vector<EdgeEdgeCandidate> ee_candidates;
        lbvh->detect_edge_edge_candidates(ee_candidates);

        std::vector<EdgeEdgeCandidate> expected_ee_candidates;
        bvh->detect_edge_edge_candidates(expected_ee_candidates);

        CHECK(ee_candidates.size() == expected_ee_candidates.size());
        // TODO: Check the candidates are the same
    }

    {
        std::vector<FaceVertexCandidate> fv_candidates;
        lbvh->detect_face_vertex_candidates(fv_candidates);

        std::vector<FaceVertexCandidate> expected_fv_candidates;
        bvh->detect_face_vertex_candidates(expected_fv_candidates);

        CHECK(fv_candidates.size() == expected_fv_candidates.size());
        // TODO: Check the candidates are the same
    }

    {
        std::vector<EdgeFaceCandidate> ef_candidates;
        lbvh->detect_edge_face_candidates(ef_candidates);

        std::vector<EdgeFaceCandidate> expected_ef_candidates;
        bvh->detect_edge_face_candidates(expected_ef_candidates);

        CHECK(ef_candidates.size() == expected_ef_candidates.size());
        // TODO: Check the candidates are the same
    }

    {
        std::vector<FaceFaceCandidate> ff_candidates;
        lbvh->detect_face_face_candidates(ff_candidates);

        std::vector<FaceFaceCandidate> expected_ff_candidates;
        bvh->detect_face_face_candidates(expected_ff_candidates);

        CHECK(ff_candidates.size() == expected_ff_candidates.size());
        // TODO: Check the candidates are the same
    }

#ifdef IPC_TOOLKIT_WITH_PROFILER
    ipc::profiler().print();
    ipc::profiler().clear();
#endif
}

// TEST_CASE("Radix sort", "[radix_sort]")
// {
//     constexpr size_t size = 100'000;
//     constexpr size_t max_value = 0xFFFFFFFF;

//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_int_distribution<uint32_t> distrib(0, max_value);

//     using MortonCode = std::array<uint32_t, 2>;
//     std::vector<MortonCode> data(size);
//     for (size_t i = 0; i < size; ++i) {
//         data[i][0] = distrib(gen);
//         data[i][1] = 0;
//     }

//     // Radix sort

//     kp::Manager mgr;

//     auto t_data = mgr.tensor(
//         data.data(), data.size(), sizeof(MortonCode),
//         kp::Memory::DataTypes::eCustom);
//     auto t_sorted_data = mgr.tensor(
//         data.size(), sizeof(MortonCode), kp::Memory::DataTypes::eCustom);

//     std::vector<std::shared_ptr<kp::Memory>> params = { {
//         mgr.tensor(1, sizeof(LBVH::AABB), kp::Memory::DataTypes::eCustom),
//         t_data,
//         t_sorted_data,
//     } };

//     auto algorithm = mgr.algorithm(
//         params,
//         SPVShader(SINGLE_RADIXSORT_SPV.begin(), SINGLE_RADIXSORT_SPV.end()),
//         kp::Workgroup({ { 1, 1, 1 } }), {},
//         std::vector<uint32_t>(1, data.size()));

//     mgr.sequence()
//         ->record<kp::OpSyncDevice>(params)
//         ->record<kp::OpAlgoDispatch>(algorithm)
//         ->record<kp::OpSyncLocal>(params)
//         ->eval();

//     std::vector<MortonCode> sorted_data =
//     t_sorted_data->vector<MortonCode>();

//     {
//         CAPTURE(sorted_data);
//         REQUIRE(
//             std::is_sorted(
//                 sorted_data.begin(), sorted_data.end(),
//                 [](const MortonCode& a, const MortonCode& b) {
//                     return a[0] < b[0];
//                 }));
//     }

//     // Compare against CPU sort
//     std::sort(
//         data.begin(), data.end(),
//         [](const MortonCode& a, const MortonCode& b) { return a[0] < b[0];
//         });
//     CHECK(data == sorted_data);
// }
