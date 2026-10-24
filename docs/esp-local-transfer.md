# Experimental local quadrature transfer

This branch adds `ESPParameters::local_quadrature_transfer` (Python: the same
attribute). It is off by default. Enable it before constructing `ESPPotential`:

```python
params = ipctk.ESPParameters(dhat=0.03, dbar_factor=0.05, area_weights=True)
params.local_quadrature_transfer = True
potential = ipctk.ESPPotential(params, use_near_far=True)
```

The first implementation supports 3D, normalized near/far integration, the
existing vertex quadrature rule, and frictionless contact. Custom face rules and
unnormalized integration are rejected. Building tangential collisions with a
nonzero friction coefficient is rejected to avoid using inconsistent legacy
friction weights. In 2D, or with `dbar_factor=0`, evaluation follows the existing
formulation. Collision detection and the endpoint mollifier are unchanged.

## Formulation

For a source triangle, rod, or particle, let `w` be its existing rest-geometry
weight and `n` its number of vertex samples. Let `N_i` and `F_i` be the near and
far point potentials at vertex `i`. Each interior edge-edge sample `j` has the
existing mollifier weight `m_j`, closest point `q_j`, and source-edge parameter
`u_j`. Its reference potential is

```
R_j = (1-u_j) N_a + u_j N_b
M = sum_j m_j
E = (w/n) [sum_i (N_i + F_i) + sum_j m_j (N(q_j) - R_j)/(1+M)]
```

Thus each EE sample receives mass `m_j/[n(1+M)]`, taken only from its source
edge's endpoints with barycentric factors `(1-u_j,u_j)`. Since `0<=u_j<=1` for
retained interior EE samples and `m_j>=0`, the total donated mass is
`M/[n(1+M)] < 1/n`: every vertex retains a nonnegative coefficient, even with
multiple EE samples. Total near-field quadrature mass stays fixed. The far-field
vertex rule is unchanged.

For a single EE sample, the weight-gradient contribution is
`(w/n)*(N(q)-R)/(1+m)^2 * grad(m)`. It multiplies a matching local residual
instead of a difference against an element-wide normalized average. If the
point potentials have compatible endpoint limits, that residual vanishes as
the EE point approaches its source endpoint.

With several samples, define `S=sum_j m_j*(N(q_j)-R_j)` and `D=1+sum_j m_j`.
The weight-gradient contribution is
`(w/n)*sum_j ((N(q_j)-R_j)-S/D)/D * grad(m_j)`. The common denominator still
couples their residuals, so cancellation of one residual does not guarantee
cancellation of its complete weight-gradient term. Transitions at a target
edge endpoint can also retain a large residual. A symmetric feature handoff
would require further quadrature and collision-builder changes.

The transfer preserves the first moment as well as total mass, because each
sample position is the same barycentric interpolation used for its donation.
A common affine point-potential field therefore integrates exactly as in the
vertex-only rule. This is an algebraic property, not an assertion that actual
nonlinear collision potentials or different contact stencils are affine.

Value, gradient, and Hessian all differentiate this expression, including `u_j`,
the reference endpoint potentials, and the denominator. Pointwise derivatives
use the existing helpers; geometry uses TinyAD; element jets retain the full
cross-stencil Hessian. Optional PSD projection is applied to the assembled
source-element Hessian before global scatter.

## Validation and limits

`[esp_local_transfer]` tests cover the reported edge-over-box values, complete
finite-difference gradients and Hessians (log and inverse-power barriers, with
and without area weighting), multiple EE samples, translation invariance, PSD
projection, empty contacts, disabled EE support, and unsupported modes. Hessian
checks use configurations away from feature boundaries: the existing formulation
is C1, and its Hessian need not be unique at an exact feature switch.

`python/examples/esp_local_transfer.py` produces an analytic old/new force sweep,
CSV, plot, and independently rebuilt-contact finite-difference checks on issue
#280's original edge-over-box geometry. In a 631-sample sweep over
`A.x=[49.700,50.015] mm`, maximum horizontal force drops from 103.1273 to 8.8900
in the reproduction's units (91.38%). A smaller target-side peak remains. The
energy and vertical force change too; this is a changed quadrature formulation.

The interior-crossing unit test checks increasing energy as the gap shrinks over
four decades. This finite test and the coefficient argument are not a proof of
global barrier behavior or mesh convergence. Keep this mode experimental while
broader scene and refinement studies are performed.
