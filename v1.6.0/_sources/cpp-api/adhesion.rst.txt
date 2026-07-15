Adhesion
========

Functions
---------

.. list-table::
    :header-rows: 1

    * - Name
      - Description
    * - :func:`normal_adhesion_potential`
      - The normal adhesion potential.
    * - :func:`normal_adhesion_potential_first_derivative`
      - The first derivative of the normal adhesion potential wrt d.
    * - :func:`normal_adhesion_potential_second_derivative`
      - The second derivative of the normal adhesion potential wrt d.
    * - :func:`max_normal_adhesion_force_magnitude`
      - The maximum normal adhesion force magnitude.
    * - :func:`tangential_adhesion_f0`
      - The tangential adhesion mollifier function.
    * - :func:`tangential_adhesion_f1`
      - The first derivative of the tangential adhesion mollifier function.
    * - :func:`tangential_adhesion_f2`
      - The second derivative of the tangential adhesion mollifier function.
    * - :func:`tangential_adhesion_f1_over_x`
      - The first derivative of the tangential adhesion mollifier function divided by y.
    * - :func:`tangential_adhesion_f2_x_minus_f1_over_x3`
      - The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.
    * - :func:`smooth_mu_a0`
      - Compute the value of the ∫ μ(y) a₁(y) dy, where a₁ is the first derivative of the smooth tangential adhesion mollifier.
    * - :func:`smooth_mu_a1`
      - Compute the value of the μ(y) a₁(y), where a₁ is the first derivative of the smooth tangential adhesion mollifier.
    * - :func:`smooth_mu_a2`
      - Compute the value of d/dy (μ(y) a₁(y)), where a₁ is the first derivative of the smooth tangential adhesion mollifier.
    * - :func:`smooth_mu_a1_over_x`
      - Compute the value of the μ(y) a₁(y) / y, where a₁ is the first derivative of the smooth tangential adhesion mollifier.
    * - :func:`smooth_mu_a2_x_minus_mu_a1_over_x3`
      - Compute the value of the [(d/dy μ(y) a₁(y)) ⋅ y - μ(y) a₁(y)] / y³, where a₁ and a₂ are the first and second derivatives of the smooth tangential adhesion mollifier.


Normal Adhesion Potential
-------------------------

.. doxygengroup:: normal_adhesion
   :content-only:

Tangential Adhesion Potential
-----------------------------

.. doxygengroup:: tangential_adhesion
   :content-only: