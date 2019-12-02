Intro to FEniCS Solver
======================

At its core, FEniCS is a package that can be used to solve numerically solve PDEs. There are many examples for how to use FEniCS to solve PDEs modelling thermal conductivity, solving the Navier-Stokes equation, etc. but at Sofa we focus on structural optimization and for that we first need to be able to model the behavior of structures under loads using the laws of continuum mechanics. We recommend reading the `Numerical tours of continuum mechanics using FEniCS <https://comet-fenics.readthedocs.io/en/latest/intro.html>`_ by Jeremy Bleyer for more advanced models using FEniCS, but the following examples cover what we will use in Sofa.

Examples:

.. toctree::
   :maxdepth: 1
   
   02_fenics_intro/2d_cantilever_beam
   02_fenics_intro/3d_cantilever_beam
   02_fenics_intro/mesh_1