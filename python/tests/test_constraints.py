import find_ipctk
import ipctk


def test_collisions():
    c = ipctk.VertexVertexCollision(0, 1)
    c.weight = 10
    assert c.weight == 10

    cs = ipctk.Candidates()
    num_candidates = 20
    cs.ev_candidates = [
        ipctk.EdgeVertexCandidate(i, i+1) for i in range(num_candidates)]
    assert len(cs) == num_candidates

    assert ipctk.EdgeVertexCandidate(0, 1) == ipctk.EdgeVertexCandidate(0, 1)
