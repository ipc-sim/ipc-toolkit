from mymath import *
import torch
import torch.nn.functional as F

# given some uv sampled on the edge (e0, e1), compute the integrand inside the Psi(x) for point p outside of edge
# alpha > 0, the smaller the less smooth
# eps > 0, the smaller the less smooth
# r >= 1, the larger the less smooth
def point_edge_potential_exact_pointwise(p, e0, e1, eps, r, alpha, uv):
    tangent = F.normalize(e1 - e0, p=2, dim=None)
    diff = (e0 - p).view(-1, 1) + torch.outer(e1 - e0, uv)
    dists = torch.norm(diff, dim=0)
    barrier_val = inv_barrier(dists, eps, r)
    Phi_val = spline3((tangent.view(1, 2) @ diff) / dists * (2 / alpha))
    return barrier_val * Phi_val.view(-1)

# discrete version of point edge potential Psi(x) using the smoothed closest point
# points: N x d
# smooth_factor >= 1, no smoothing at all if = 1
# alpha > 0, the smaller the less smooth
# eps > 0, the smaller the less smooth
# r >= 1, the larger the less smooth
def point_edge_potential_discrete(points, e0, e1, eps, r, alpha, smooth_factor):
    tangent = e1 - e0
    length = torch.norm(tangent)
    normal = torch.tensor([-tangent[1], tangent[0]]) / length

    # uv of point projected onto the line
    s = ((points - e0.view(1, -1)) @ (tangent / length**2)).view(-1)
    L = L_ns(s)**smooth_factor
    # smooth closest point
    y = e0 + torch.outer((s - L).view(-1), tangent)

    Phi = (F.normalize(y - points, p=2, dim=1) @ tangent.view(-1,1)).view(-1) ** 2 / length**2
    dist_sqr = ((points - e0) @ normal.view(-1,1)).view(-1) **2 + length**2 * L**2
    
    return torch.norm(e1 - e0) * spline3(Phi * (2 / alpha)) * inv_barrier(dist_sqr, eps, r)
