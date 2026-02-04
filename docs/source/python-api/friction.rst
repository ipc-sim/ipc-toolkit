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

``anisotropic_mu_eff_f`` and ``anisotropic_mu_eff_f_dtau`` implement the
elliptical L2 (matchstick) model (:cite:t:`Erleben2019Matchstick`). The C++
API provides ``anisotropic_mu_eff_from_tau_aniso``, ``anisotropic_mu_eff_f_grad``,
and related helpers; the solver uses them when you set anisotropic coefficients
on tangential collisions.

.. autofunction:: ipctk.anisotropic_mu_eff_f
.. autofunction:: ipctk.anisotropic_mu_eff_f_dtau