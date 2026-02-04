Friction
========

.. seealso::

    :doc:`/tutorials/advanced_friction` explains the friction model, the
    static/kinetic transition, and anisotropic usage. Full derivation and
    plots are in ``notebooks/anisotropic_friction_math.ipynb``.

Smooth Mollifier
----------------

.. doxygenfunction:: smooth_friction_f0
.. doxygenfunction:: smooth_friction_f1
.. doxygenfunction:: smooth_friction_f2
.. doxygenfunction:: smooth_friction_f1_over_x
.. doxygenfunction:: smooth_friction_f2_x_minus_f1_over_x3

Smooth :math:`\mu`
------------------

.. doxygenfunction:: smooth_mu
.. doxygenfunction:: smooth_mu_derivative

.. doxygenfunction:: smooth_mu_f0
.. doxygenfunction:: smooth_mu_f1
.. doxygenfunction:: smooth_mu_f2
.. doxygenfunction:: smooth_mu_f1_over_x
.. doxygenfunction:: smooth_mu_f2_x_minus_mu_f1_over_x3

Anisotropic Friction Helpers
-----------------------------

Effective friction follows an elliptical L2 projection (matchstick cone):
:math:`\mu_{\text{eff}} = \sqrt{(\mu_0 t_0)^2 + (\mu_1 t_1)^2}` with
:math:`t = \tau / \|\tau\|`. Use ``anisotropic_mu_eff_from_tau_aniso`` when you
have :math:`\tau_{\text{aniso}}` and need :math:`\mu_s`, :math:`\mu_k` for the
smooth transition; use ``anisotropic_mu_eff_f`` when you have the unit
direction. Zero ``mu_s_aniso`` and ``mu_k_aniso`` falls back to scalar
:math:`\mu_s`, :math:`\mu_k`. See :cite:t:`Erleben2019Matchstick` for the
Matchstick model; code: `erleben/matchstick <https://github.com/erleben/matchstick>`_.

.. doxygenfunction:: anisotropic_mu_eff_f
.. doxygenfunction:: anisotropic_mu_eff_f_dtau
.. doxygenfunction:: anisotropic_x_from_tau_aniso
.. doxygenfunction:: anisotropic_mu_eff_from_tau_aniso
.. doxygenfunction:: anisotropic_mu_eff_f_grad