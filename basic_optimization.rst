Intro to dolfin-adjoint
===========================
dolfin-adjoint is a tool used alongside dolfin (FEniCS) for optimization. The basic idea is that a control variable is defined by the user and dolfin can be used to solve some PDE based on the control variable. The user defines an objective function (generally one which should be minimized) based on the results of the dolfin simulation and dolfin-adjoint optimizes the control variable to best fit the objective while staying within user defined constraints.
Examples:

.. toctree::
   :maxdepth: 1

   03_basic_optimization/2d_cantilever