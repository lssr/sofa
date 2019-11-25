Fenics Region Selector Modules
==============================

Fenics has several methods for defining regions for boundary condition control
and we use a method where we define a class for each domain. These classes have
a method which takes a point and returns true or false based on whether the point
qualifies for the domain. Domains can either be a region or a boundary. For
example, to define a domain for all points in 2D between :math:`x=0` and :math:`x=2`
and :math:`y=1` and :math:`y=4`, the following class would be written::

 class GetRegion(SubDomain):
     def __init__(self):
         super().__init__()
 
     def inside(self, x, on_boundary):
         return between(x[0], (0.0,2.0)) and between(x[1], (1.0,4.0))

The vector ``x`` passed into the ``inside`` function is of the form :math:`\vec{x}=
\begin{pmatrix} x \\ y \\ z \end{pmatrix}`. If we wanted to define several rectangles
for the boundary conditions, we could either write a separate class for each one with
hardcoded coordinate values, or we could write one class that users can pass parameters
into to generate multiple unique rectangular boundaries. Here is a simple example::

 class GetRectangularRegion(SubDomain):
     def __init__(self, x1, y1, x2, y2):
         super().__init__()
         self.x1 = x1
         self.y1 = y1
         self.x2 = x2
         self.y2 = y2
 
     def inside(self, x, on_boundary):
         return between(x[0], (self.x1, self.x2)) and \
             between(x[1], (self.y1, self.y2))

This class can now be used to select any rectangular region whose sides are parallel to
the x and y axes. This is obviously somewhat limited, particularly since there is only one
set of accepted inputs, so we have developed some more robust modules to select regions.
They will be used throughout the rest of the project, and their documentation can be found
here:

.. toctree::
   :maxdepth: 1

   region_selectors/region_selector_2d
   region_selectors/region_selector_3d