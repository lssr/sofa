Optimization for Rating
==========================
-------------
Introduction
-------------
So far, our optimization has been to maximize performance while staying under a certain weight. This makes sense from a cost-minimizing perspective, but sometimes we may want to meet some load rating no matter how much material it would take, though we still want to minimize the weight (and cost) of the part. We can achieve this in two ways: we can either have the objective function for performance be to minimize the distance to a certain rating rather than to globally maximize performance and have a secondary objective of minimizing weight, or we can have a single objective for weight-minimization and set a constraint on the performance that it must meet some threshold. We will examine both of these options in the following examples.

.. toctree::
   :maxdepth; 1
   
   04_rating_optimization/2d_cross_section_frame