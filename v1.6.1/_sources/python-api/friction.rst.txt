Friction
========

.. seealso::

    :doc:`/tutorials/advanced_friction` describes the friction model and
    anisotropic usage. Additional anisotropic helpers are in :doc:`/cpp-api/friction`.

Smooth Mollifier
----------------

.. autofunction:: ipctk.smooth_friction_f0
.. autofunction:: ipctk.smooth_friction_f1
.. autofunction:: ipctk.smooth_friction_f2
.. autofunction:: ipctk.smooth_friction_f1_over_x
.. autofunction:: ipctk.smooth_friction_f2_x_minus_f1_over_x3

Smooth :math:`\mu`
------------------

.. autofunction:: ipctk.smooth_mu
.. autofunction:: ipctk.smooth_mu_derivative
.. autofunction:: ipctk.smooth_mu_f0
.. autofunction:: ipctk.smooth_mu_f1
.. autofunction:: ipctk.smooth_mu_f2
.. autofunction:: ipctk.smooth_mu_f1_over_x
.. autofunction:: ipctk.smooth_mu_f2_x_minus_mu_f1_over_x3

Anisotropic Friction Helpers
-----------------------------

``anisotropic_mu_eff_f`` implements the elliptical :math:`L^2` (matchstick)
effective-μ formula (:cite:t:`Erleben2019Matchstick`). The C++ API provides
``anisotropic_mu_eff_from_tau_aniso`` and related helpers. For
direction-dependent ellipse axes on tangential collisions, call
``TangentialCollisions.update_lagged_anisotropic_friction_coefficients`` so
effective μ matches the lagged slip direction; the built-in friction paths then
use those lagged scalars (they do not differentiate μ with respect to slip).

.. autofunction:: ipctk.anisotropic_mu_eff_f