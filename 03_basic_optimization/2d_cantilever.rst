2D Cantilever Optimization
==========================

Introduction
------------

Suppose we want to make an aluminum shelf no taller than 10cm to best support an :math:`4\text{cm}` wide and :math:`80\text{cm}` long rod weighing :math:`500\text{N}` :math:`23\text{cm}` away from a wall while keeping the mass of the shelf under :math:`16.2\text{kg}`. This shelf will simply be a 2D profile extruded out of the page :math:`80\text{cm}` for ease of manufacturing, so we can just focus on optimizing the 2D profile. In order to ensure the rod is supported, we require that :math:`1\text{cm}` of material is present under it through its whole width. What is the best shape for the shelf to have? The design space can be summarized in a sketch:

.. image:: 2d_cantilever/design_domain.png

This is known as the design domain of the problem. We want to figure out how to best support the rod while keeping the total mass of the shelf under :math:`16.2\text{kg}`. How do we express this mathematically?
First, let's find what fraction of our design domain can be filled with material. At any location, we can either have aluminum or not - there is no partial density. This means that :math:`16.2\text{kg} * \frac{1 \text{cm}^3}{2.7 g} = 6000 \text{cm}^3` can be filled with aluminum, which works out to :math:`\frac{6000 \text{cm}^3}{25\text{cm} * 10\text{cm} * 80 \text{cm} = .3`, so up to 30% of the design domain can be filled with aluminum.
We also need to define "best support" as an objective function. We can do this by minimizing the compliance which can be defined as the sum of the strain energy in the part. Since we are operating in the elastic region we know that the strain energy in the part is equal to the input energy, or

.. math::
   \text{J} = \int_\Omega \langle \boldsymbol{F}, \boldsymbol{v} \rangle \text{dA}

Where :math:`\boldsymbol{F}` represents the applied forces at the subsurfaces :math:`\text{dA}` on the applied load region :math:`\Omega` and :math:`\boldsymbol{v}` represent the displacements on :math:`\text{dA}`. Importantly, :math:`\text{J}` is a single scalar value so it is simple to optimize. :math:`\text{J}` is a function of :math:`\boldsymbol{\rho}`, the density distribution in the region, which is a continuous value between 0 and 1. In the end we ideally want :math:`\boldsymbol{\rho}` to converge to 0 or 1 on each element.
In summary, we want to vary :math:`\boldsymbol{\rho}` to minimize :math:`\text{J}(\boldsymbol{\rho})` while keeping :math:`\int_V \rho \text{dV} < 0.3\int_V 1 \text{dV}`, where V is the domain space. Further, :math:`\forall (x,y | x \in [21, 25] \cap y \in [9, 10]) \rightarrow \rho(x,y) = 1`.