2D Region Selector
==================

This file contains definitions for different 2D regions we may want to select.

------------------
Rectangular Region
------------------
This class selects a rectangular region whose sides are parallel to the :math:`x` and 
:math:`y` axes. It has two classmethods to allow it to be defined with either a pair of 
points or four floats. In general, it is recommended that the classes be created through 
the classmethods.

::

    class GetRectangularRegion(SubDomain):

The internal object init function takes four floats representing the two opposing corners 
of the rectangle. It sets ``x1`` and ``y1`` to the min of the two given values for x and y 
and ``x2`` and ``y2`` to the max so that the ``between()`` function can be used simply.

::

    def __init__(self, x1:float, y1:float, x2:float, y2:float):
        super().__init__()
        self.x1 = np.minimum(x1, x2)
        self.y1 = np.minimum(y1, y2)
        self.x2 = np.maximum(x1, x2)
        self.y2 = np.maximum(y1, y2)

This classmethod allows for region declaration using a pair of points. Their coordinates 
are used to call the main init function.

::

    @classmethod
    def from_points(cls, p1:Point, p2:Point) -> 'GetRectangularRegion':
        return cls(p1.x(), p1.y(), p2.x(), p2.y())

This classmethod gives an alias classmethod for the object init function.

::

    @classmethod
    def from_floats(cls, x1:float, y1:float, x2:float, y2:float) -> 'GetRectangularRegion':
        return cls(x1, y1, x2, y2)

The ``inside()`` method for this class simply returns whether a given :math:`\vec{x}` is between ``x1`` and ``x2`` and between ``y1`` and ``y2``.

::

    def inside(self, x, on_boundary):
        return(between(x[0], (self.x1, self.x2)) and \
                between(x[1], (self.y1, self.y2))

------------------
Circular Region
------------------

This class selects circular regions. Its internal init method takes three floats describing 
the center of the circle and its radius. It also has two classmethods which can be used to 
instantiate the class, one which takes three floats and one which takes a point and a float.

::

    class GetCircularRegion(SubDomain):
        def __init__(self, cx:float, cy:float, r:float):
            self.cx = cx
            self.cy = cy
            self.r = r
       
        @classmethod
        def from_points(cls, c:Point, r:float) -> 'GetCircularRegion':
            return cls(c.x(), c.y(), r)
       
        @classmethod
        def from_floats(cls, cx:float, cy:float, r:float) -> 'GetCircularRegion':
            return cls(cx, cy, r)
    
        def inside(self, x, on_boundary):
            return (x[0] - self.cx)**2.0 + (x[1] - self.cy)**2.0 <= r**2.0

------------------
Linear Boundary
------------------

------------------
Complete Code
------------------
The complete code follows and can also be downloaded :download:`here </../code/region_selector_2d.py>`.

.. literalinclude:: /../code/region_selector_3d.py
   :language: python