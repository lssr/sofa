2D Cross Section Optimization for Frame
=======================================
------------
Introduction
------------

Suppose we want to create a 2D profile which can be extruded to form members of a frame. These members will take axial compressive loads and a single transverse load, meaning that the profile will be optimized for both maximum buckling force and for least compliance in a single direction, since the beam could be rotated axially to align the single expected shear force direction with the "strong" axis of the profile.
For a beam of length :math:`\text{L}` and constant cross section A, our objective function will be based on the maximum bending load :math:`P_{\tau}(L, A)` and our maximum buckling load :math:`P_{N}(L, A)`.
Euler-Bernoulli theory tells us that the stress in a bending beam is equal to:

.. math::
   \sigma = \frac{Mz}{I}

Where :math:`M` is the bending moment, :math:`z` is the distance from the central axis of the profile, and :math:`I` is the moment of area of the profile. Since the beam is fixed on both ends and assuming the worst case of the load being applied at the center of the beam, we know:

.. math::
   M = \frac{PL}{8}

Where :math:`P` is the applied bending load and :math:`L` is the length of the beam. Since :math:`z` is maximized at the furthest point from the central axis of the profile, we know:

.. math::
   \sigma_{max} = \frac{\frac{PL}{8}\frac{T}{2}}{I} = \frac{PLT}{16I}

Where T is the thickness of the beam in the direction of loading. If we set :math:`\sigma_{max} = :math:`\sigma_y`, the yield stress of the beam material, we can solve for the maximum bending load as:

.. math::
   P_\tau = \frac{16\sigma_y I}{LT}

Euler-Bernoulli beam theory also tells us that the max buckling load of a beam is:

.. math::
   \P_N = \frac{4\pi^2EI}{L^2}

When the bending load is applied in the y direction, we get :math:`J = P_\tau + P_N` as:

.. math::
   J = \frac{16\sigma_y}{LT}I_x + \frac{4\pi^2E}{L^2}\text{min}(I_x, I_y)