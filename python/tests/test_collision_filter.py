import find_ipctk
import numpy as np
from ipctk import (
    CollisionFilter,
    CollisionMesh,
    make_connected_components_filter,
    make_sparse_filter,
    make_static_obstacle_filter,
    make_vertex_patches_filter,
)

# ─────────────────────────────────────────────────────────────────────────────
# CollisionFilter basic construction
# ─────────────────────────────────────────────────────────────────────────────


def test_collision_filter_default():
    """Default-constructed filter accepts all pairs."""
    f = CollisionFilter()
    assert f(0, 1)
    assert f(5, 10)
    assert f(0, 0)  # self-pair accepted by default


def test_collision_filter_from_callable():
    """Filter constructed from a Python callable."""
    # Two bodies: verts 0-2 in body 0, verts 3-5 in body 1.
    body = [0, 0, 0, 1, 1, 1]
    f = CollisionFilter(lambda i, j: body[i] != body[j])

    assert not f(0, 1)  # same body
    assert not f(3, 5)  # same body
    assert f(0, 3)  # different bodies
    assert f(2, 4)  # different bodies


# ─────────────────────────────────────────────────────────────────────────────
# Composition operators
# ─────────────────────────────────────────────────────────────────────────────


def test_collision_filter_negation():
    """~ inverts the filter result."""
    body = [0, 0, 0, 1, 1, 1]
    same_body = CollisionFilter(lambda i, j: body[i] == body[j])
    cross_body = ~same_body

    assert same_body(0, 2)
    assert not cross_body(0, 2)

    assert not same_body(0, 4)
    assert cross_body(0, 4)


def test_collision_filter_union():
    """f | g accepts a pair if EITHER filter passes."""
    # Layer tags: 0,0,1,0,1,1
    layer = [0, 0, 1, 0, 1, 1]

    a = CollisionFilter(lambda i, j: layer[i] == 1)  # pass if vi in layer 1
    b = CollisionFilter(lambda i, j: layer[j] == 0)  # pass if vj in layer 0

    a_or_b = a | b

    assert a_or_b(2, 5)  # a passes (vi=2, layer 1)
    assert a_or_b(0, 3)  # b passes (vj=3, layer 0)
    assert not a_or_b(0, 4)  # neither: vi=0 layer 0, vj=4 layer 1


def test_collision_filter_intersection():
    """f & g accepts a pair only if BOTH filters pass."""
    body = [0, 0, 0, 1, 1, 1]
    layer = [0, 0, 1, 0, 1, 1]

    self_col = CollisionFilter(lambda i, j: body[i] != body[j])
    layer_col = CollisionFilter(lambda i, j: layer[i] != layer[j])

    both = self_col & layer_col

    # (0,1): same body → self_col blocks → intersection blocks
    assert not both(0, 1)
    # (0,2): same body → self_col blocks → intersection blocks
    assert not both(0, 2)
    # (0,3): cross body, same layer (0,0) → layer_col blocks → blocks
    assert not both(0, 3)
    # (0,4): cross body, cross layer (0,1) → both pass
    assert both(0, 4)


def test_collision_filter_union_intersection_semantics():
    """Union passes when either passes; intersection requires both."""
    body = [0, 0, 0, 1, 1, 1]
    layer = [0, 0, 1, 0, 1, 1]

    self_col = CollisionFilter(lambda i, j: body[i] != body[j])
    layer_col = CollisionFilter(lambda i, j: layer[i] != layer[j])

    cf_union = self_col | layer_col
    cf_inter = self_col & layer_col

    # (0,1): both block → union and intersection both block
    assert not cf_union(0, 1)
    assert not cf_inter(0, 1)

    # (0,2): self blocks, layer passes → union passes, inter blocks
    assert cf_union(0, 2)
    assert not cf_inter(0, 2)

    # (0,3): self passes, layer blocks → union passes, inter blocks
    assert cf_union(0, 3)
    assert not cf_inter(0, 3)

    # (0,4): both pass
    assert cf_union(0, 4)
    assert cf_inter(0, 4)


def test_collision_filter_ior():
    """f |= g modifies f to be the union of f and g."""
    f = CollisionFilter(lambda i, j: i < 3)
    f |= CollisionFilter(lambda i, j: j < 3)

    assert f(5, 2)  # j < 3
    assert f(2, 5)  # i < 3
    assert not f(5, 5)  # neither


def test_collision_filter_iand():
    """f &= g modifies f to be the intersection of f and g."""
    f = CollisionFilter(lambda i, j: i < 10)
    f &= CollisionFilter(lambda i, j: j < 10)

    assert f(3, 4)  # both
    assert not f(3, 15)  # j >= 10
    assert not f(15, 3)  # i >= 10


# ─────────────────────────────────────────────────────────────────────────────
# Factory functions
# ─────────────────────────────────────────────────────────────────────────────


def test_make_vertex_patches_filter():
    """Blocks pairs within the same patch, allows cross-patch pairs."""
    patches = np.array([0, 0, 0, 1, 1, 1], dtype=np.int32)
    f = make_vertex_patches_filter(patches)

    assert f(0, 3)  # different patches
    assert f(2, 4)  # different patches
    assert not f(0, 1)  # same patch 0
    assert not f(3, 5)  # same patch 1


def test_make_static_obstacle_filter():
    """Blocks static-static pairs; allows dynamic-* pairs."""
    # Vertices 0-3 are dynamic; vertices 4-5 are static obstacles.
    n_dynamic = 4
    f = make_static_obstacle_filter(n_dynamic)

    # dynamic–dynamic: allowed
    assert f(0, 3)
    assert f(1, 2)

    # dynamic–static: allowed (either order)
    assert f(2, 4)
    assert f(0, 5)
    assert f(4, 1)  # order flipped

    # static–static: blocked
    assert not f(4, 5)
    assert not f(5, 4)


def test_make_connected_components_filter():
    """Blocks intra-component pairs; allows cross-component pairs."""
    # Two disconnected triangles: {0,1,2} and {3,4,5}
    faces = np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int32)
    f = make_connected_components_filter(faces)

    # Cross-component: allowed
    assert f(0, 3)
    assert f(1, 4)
    assert f(2, 5)

    # Intra-component: blocked
    assert not f(0, 1)
    assert not f(0, 2)
    assert not f(1, 2)
    assert not f(3, 4)
    assert not f(3, 5)
    assert not f(4, 5)


# ─────────────────────────────────────────────────────────────────────────────
# Composition chain
# ─────────────────────────────────────────────────────────────────────────────


def test_collision_filter_composition_chain():
    """Combine connected-components and static-obstacle filters."""
    # Two disconnected triangles: {0,1,2} and {3,4,5}
    faces = np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int32)

    # Vertices 4 and 5 are also "static obstacles"
    n_dynamic = 4

    no_self = make_connected_components_filter(faces)
    no_static = make_static_obstacle_filter(n_dynamic)

    active = no_self & no_static

    # Cross-component, at least one dynamic → allowed
    assert active(0, 3)
    assert active(1, 4)  # vj=4 is static but vi=1 is dynamic
    assert active(2, 5)  # vj=5 is static but vi=2 is dynamic

    # Static–static cross-component pair → blocked by no_static
    assert not active(4, 5)

    # Intra-component → blocked by no_self
    assert not active(0, 1)
    assert not active(3, 5)


# ─────────────────────────────────────────────────────────────────────────────
# SparseCanCollide composition
# ─────────────────────────────────────────────────────────────────────────────


def test_make_sparse_filter():
    """make_sparse_filter returns a composable CollisionFilter."""
    # Block pair (0, 1) explicitly; allow everything else.
    sparse = make_sparse_filter({(0, 1): False}, True)

    assert sparse(0, 2)      # not in map → default True
    assert not sparse(0, 1)  # explicit False

    # ~ (invert): (0,1) now passes; everything else is blocked
    inverted = ~sparse
    assert inverted(0, 1)      # was blocked, now passes
    assert not inverted(0, 2)  # was allowed, now blocked

    # | (union): the blocked (0,1) pair is rescued by a second filter
    rescue_01 = CollisionFilter(lambda i, j: i == 0 and j == 1)
    rescued = sparse | rescue_01
    assert rescued(0, 1)  # sparse blocks but rescue_01 passes → union passes
    assert rescued(0, 2)  # sparse passes → union passes

    # & (intersection): add a restriction on top of sparse
    only_small = CollisionFilter(lambda i, j: i < 3 and j < 3)
    restricted = sparse & only_small
    assert restricted(0, 2)      # sparse allows, only_small allows → passes
    assert not restricted(0, 1)  # sparse blocks → intersection fails
    assert not restricted(0, 5)  # only_small fails (j=5 ≥ 3) → intersection fails


def test_collision_filter_on_mesh_can_collide():
    """CollisionFilter assigned to mesh.can_collide works correctly."""
    V = np.array(
        [[0, 0], [1, 0], [0, 1], [1, 1], [2, 0], [3, 0], [2, 1], [3, 1]], dtype=float
    )
    E = np.array(
        [[0, 1], [1, 3], [3, 2], [2, 0], [4, 5], [5, 7], [7, 6], [6, 4]], dtype=int
    )
    mesh = CollisionMesh(V, E)

    patches = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=np.int32)
    f = make_vertex_patches_filter(patches)

    mesh.can_collide = f

    for i in range(V.shape[0]):
        for j in range(V.shape[0]):
            assert mesh.can_collide(i, j) == (patches[i] != patches[j])
