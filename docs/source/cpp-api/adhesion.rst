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


Normal Adhesion Potential
-------------------------

.. doxygengroup:: normal_adhesion
   :content-only:

Tangential Adhesion Potential
-----------------------------

.. doxygengroup:: tangential_adhesion
   :content-only: