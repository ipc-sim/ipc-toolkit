#include <catch2/catch_test_macros.hpp>

#include <ipc/collision_filter.hpp>

#include <Eigen/Core>
#include <vector>

using namespace ipc;

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter default", "[collision_filter]")
{
    CollisionFilter f;
    CHECK(f(0, 1));
    CHECK(f(5, 10));
    CHECK(f(0, 0)); // self-pair is accepted by default
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter from lambda", "[collision_filter]")
{
    // Two bodies: verts 0-2 in body 0, verts 3-5 in body 1.
    const std::vector<int> body = { 0, 0, 0, 1, 1, 1 };
    CollisionFilter self_col(
        [&body](size_t vi, size_t vj) { return body[vi] != body[vj]; });

    CHECK_FALSE(self_col(0, 1)); // same body
    CHECK_FALSE(self_col(3, 5)); // same body
    CHECK(self_col(0, 3));       // different bodies
    CHECK(self_col(2, 4));       // different bodies
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter operator!", "[collision_filter]")
{
    CollisionFilter same_body(
        [](size_t vi, size_t vj) { return vi / 3 == vj / 3; });

    CollisionFilter cross_body = !same_body;

    CHECK(same_body(0, 2));
    CHECK_FALSE(cross_body(0, 2));

    CHECK_FALSE(same_body(0, 4));
    CHECK(cross_body(0, 4));
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter operator|", "[collision_filter]")
{
    // Layer tags: 0,0,1,0,1,1
    const std::vector<int> layer = { 0, 0, 1, 0, 1, 1 };

    CollisionFilter a( // pass if vi is in layer 1
        [&layer](size_t vi, size_t /*vj*/) { return layer[vi] == 1; });
    CollisionFilter b( // pass if vj is in layer 0
        [&layer](size_t /*vi*/, size_t vj) { return layer[vj] == 0; });

    CollisionFilter a_or_b = a | b;

    CHECK(a_or_b(2, 5));       // a passes (vi=2, layer 1)
    CHECK(a_or_b(0, 3));       // b passes (vj=3, layer 0)
    CHECK_FALSE(a_or_b(0, 4)); // neither: vi=0 layer 0, vj=4 layer 1
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter operator&", "[collision_filter]")
{
    const std::vector<int> body = { 0, 0, 0, 1, 1, 1 };
    const std::vector<int> layer = { 0, 0, 1, 0, 1, 1 };

    CollisionFilter self_col(
        [&body](size_t vi, size_t vj) { return body[vi] != body[vj]; });
    CollisionFilter layer_col(
        [&layer](size_t vi, size_t vj) { return layer[vi] != layer[vj]; });

    CollisionFilter both = self_col & layer_col;

    // (0,1): same body → self_col blocks → intersection blocks
    CHECK_FALSE(both(0, 1));
    // (0,2): same body → self_col blocks → intersection blocks
    CHECK_FALSE(both(0, 2));
    // (0,3): cross body, same layer (0,0) → layer_col blocks → blocks
    CHECK_FALSE(both(0, 3));
    // (0,4): cross body, cross layer (0,1) → both pass
    CHECK(both(0, 4));
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter union/intersection semantics", "[collision_filter]")
{
    const std::vector<int> body = { 0, 0, 0, 1, 1, 1 };
    const std::vector<int> layer = { 0, 0, 1, 0, 1, 1 };

    CollisionFilter self_col(
        [&body](size_t vi, size_t vj) { return body[vi] != body[vj]; });
    CollisionFilter layer_col(
        [&layer](size_t vi, size_t vj) { return layer[vi] != layer[vj]; });

    CollisionFilter cf_union = self_col | layer_col;
    CollisionFilter cf_inter = self_col & layer_col;

    // (0,1): both block → union also blocks
    CHECK_FALSE(cf_union(0, 1));
    CHECK_FALSE(cf_inter(0, 1));

    // (0,2): self blocks, layer passes → union passes, inter blocks
    CHECK(cf_union(0, 2));
    CHECK_FALSE(cf_inter(0, 2));

    // (0,3): self passes, layer blocks → union passes, inter blocks
    CHECK(cf_union(0, 3));
    CHECK_FALSE(cf_inter(0, 3));

    // (0,4): both pass
    CHECK(cf_union(0, 4));
    CHECK(cf_inter(0, 4));
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter operator|=", "[collision_filter]")
{
    CollisionFilter f([](size_t vi, size_t /*vj*/) { return vi < 3; });
    f |= CollisionFilter([](size_t /*vi*/, size_t vj) { return vj < 3; });

    CHECK(f(5, 2));       // vj < 3
    CHECK(f(2, 5));       // vi < 3
    CHECK_FALSE(f(5, 5)); // neither
}

TEST_CASE("CollisionFilter operator&=", "[collision_filter]")
{
    CollisionFilter f([](size_t vi, size_t /*vj*/) { return vi < 10; });
    f &= CollisionFilter([](size_t /*vi*/, size_t vj) { return vj < 10; });

    CHECK(f(3, 4));        // both
    CHECK_FALSE(f(3, 15)); // vj >= 10
    CHECK_FALSE(f(15, 3)); // vi >= 10
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("make_vertex_patches_filter", "[collision_filter]")
{
    Eigen::VectorXi patches(6);
    patches << 0, 0, 0, 1, 1, 1;

    CollisionFilter f = make_vertex_patches_filter(patches);

    CHECK(f(0, 3));       // different patches
    CHECK(f(2, 4));       // different patches
    CHECK_FALSE(f(0, 1)); // same patch 0
    CHECK_FALSE(f(3, 5)); // same patch 1
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("make_static_obstacle_filter", "[collision_filter]")
{
    // Vertices 0-3 are dynamic; vertices 4-5 are static obstacles.
    const size_t n_dynamic = 4;
    CollisionFilter f = make_static_obstacle_filter(n_dynamic);

    // dynamic–dynamic: allowed
    CHECK(f(0, 3));
    CHECK(f(1, 2));

    // dynamic–static: allowed
    CHECK(f(2, 4));
    CHECK(f(0, 5));
    CHECK(f(4, 1)); // order flipped

    // static–static: blocked
    CHECK_FALSE(f(4, 5));
    CHECK_FALSE(f(5, 4));
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("make_connected_components_filter", "[collision_filter]")
{
    // Two disconnected triangles: {0,1,2} and {3,4,5}
    Eigen::MatrixXi faces(2, 3);
    faces << 0, 1, 2, 3, 4, 5;

    CollisionFilter f = make_connected_components_filter(faces);

    // Cross-component pairs: allowed
    CHECK(f(0, 3));
    CHECK(f(1, 4));
    CHECK(f(2, 5));

    // Intra-component pairs: blocked
    CHECK_FALSE(f(0, 1));
    CHECK_FALSE(f(0, 2));
    CHECK_FALSE(f(1, 2));
    CHECK_FALSE(f(3, 4));
    CHECK_FALSE(f(3, 5));
    CHECK_FALSE(f(4, 5));
}

// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("CollisionFilter composition chain", "[collision_filter]")
{
    // Combine: cross-component AND not-static-static
    Eigen::MatrixXi faces(2, 3);
    faces << 0, 1, 2, 3, 4, 5;

    // Verts 4 and 5 are also "static obstacles"
    const size_t n_dynamic = 4;

    CollisionFilter no_self = make_connected_components_filter(faces);
    CollisionFilter no_static = make_static_obstacle_filter(n_dynamic);

    CollisionFilter active = no_self & no_static;

    // Cross-component, at least one dynamic → allowed
    CHECK(active(0, 3));
    CHECK(active(1, 4)); // vj=4 is static but vi=1 is dynamic
    CHECK(active(2, 5)); // vj=5 is static but vi=2 is dynamic

    // Static–static cross-component pair → blocked by no_static
    CHECK_FALSE(active(4, 5));

    // Intra-component → blocked by no_self
    CHECK_FALSE(active(0, 1));
    CHECK_FALSE(active(3, 5));
}
