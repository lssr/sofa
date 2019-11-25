3D Region Selector
==================

This file contains definitions for different 3D regions we may want to select.

This file is going to require some more vector operations, so we define two helper 
functions to simplify our code later on. We import dolfin and numpy and redefine a larger 
tolerance since there's greater room for error in these calculations. First, we write a 
function to find the magnitude of the cross product of two vectors::

    from dolfin import *
    import numpy as np
    
    tol = 0.000001
    
    def mag_cross(v1, v2):
        v1xv2 = np.cross(v1, v2)
        return np.sqrt(np.dot(v1xv2,v1xv2))

We also want a function that can check if two vectors are equal without testing for exact 
equality, essentially a version of dolfin's ``near()`` but for vectors. For each element in 
the vector, check if it matches the corresponding element in the other. If not, return 
False. If none of them don't match, then the vectors are equal so we return True::

    def vec_near(v1, v2, tol = DOLFIN_EPS):
        for i in range(0, v1.shape[0]):
            if not near(v1[i], v2[i], tol):
                return False
        return True

----------------------
Cuboid Region Selector
----------------------
This class selects rectangular cuboid region in any orientation. The cuboids are defined 
using 4 points, where one point is a corner of the cuboid and the others are the three 
points adjacent to it.

::

    class GetCuboidRegion(SubDomain):
        def __init__(self, p1:Point, p2:Point, p3:Point, p4:Point):
            super().__init__()
            self.p1 = p1
            self.p2 = p2
            self.p3 = p3
            self.p4 = p4
        
We now define the three edges from ``p1`` to ``p2``, ``p3``, ``p4`` as the vectors 
:math:`\vec{u}, \vec{v}, \vec{w}`, respectively. A point with position vector 
:math:`\vec{x}` is inside the region if and only if its projection onto :math:`\vec{u}, 
\vec{v}, \vec{w}` is greater than :math:`0` and less than the lengths of :math:`\vec{u}, 
\vec{v}, \vec{w}`.

::

        def inside(self, x, on_boundary):
            p1c = self.p1.array()
            p2c = self.p2.array()
            p3c = self.p3.array()
            p4c = self.p4.array()
            
            u = p2c - p1c
            v = p3c - p1c
            w = p4c - p1c
            
            x_vec = x - p1c
            
            ux = np.dot(u, x_vec)
            vx = np.dot(v, x_vec)
            wx = np.dot(w, x_vec)
            
            um = np.dot(u,u)
            vm = np.dot(v,v)
            wm = np.dot(w,w)
            
            return between(ux, (0, um)) and between(vx, (0, vm)) and between(wx, (0, wm))

---------------------------
Cylindrical Region Selector
---------------------------
This class selects a cylindrical region as defined by two points representing the centers of the circles at the end of the cylinder and a float representing the radius of the circle. A point :math:`\vec{x}` is inside the region if and only if the projection of the vector from one point on the cylinder to :math:`\vec{x}` onto the central axis of the cylinder :math:`\vec{a}` is between :math:`0` and :math:`|\vec{a}|` and its distance from :math:`\vec{a}` is between :math:`0` and the radius.

::

    class GetCylindricalRegion(SubDomain):
        def __init__(self, p1:Point, p2:Point, r:float):
            super.__init__()
            self.p1 = p1
            self.p2 = p2
            self.r = r
        
        def inside(self, x, on_boundary):
            p1c = self.p1.array()
            p2c = self.p2.array()
            
            a = p2c - p1c
            xv = x - p1c
            
            xva = np.dot(xv, a)
            
            xvm = np.dot(xv, xv)
            am = np.dot(a, a)
            xam = (xva**2.0)/am
            dsq = xvm - xam
            
            return between(dsq, (0.0, self.r**2.0)) and between(xva, (0.0, am))

------------------------
Planar Boundary Selector
------------------------
This class selects a boundary in the shape of a parallelogram. It can be defined in two ways. One way is to select a coordinate value, say :math:`x=5`, and it selects the points that lay in that plane. The other way is to define a parallelogram region using its four corners as points and the class selects the points which lie in-plane to the parallelogram and are within its boundaries. This class can be used to select a fixed end (e.g. by selecting the plane :math:`x=0`) or to select a loading region which is only a subsection of a face, for example the last 2mm of the top face of a rectangular beam.

::

    class GetPlanarBoundary(SubDomain):
        def __init__(self, p1:Point, p2:Point, p3:Point, p4:Point):
            super().__init__()
            v12 = p2.array() - p1.array()
            v13 = p3.array() - p1.array()
            v14 = p4.array() - p1.array()
            
            self.p1 = p1
            

Here, the code figures out which of the points is across from ``p1`` and stores it as ``p3`` using the parallelogram law of vector addition. If none of the points are registered as being the one that's diagonal, the four points do not form a parallelogram.

::

            if vec_near(v13 + v14, v12, tol): # p2 is across from p1
                self.p2 = p3
                self.p3 = p2
                self.p4 = p4
            elif vec_near(v12 + v14, v13, tol): # p3 is across from p1
                self.p2 = p2
                self.p3 = p3
                self.p4 = p4
            elif vec_near(v12 + v13, v14, tol): # p4 is across from p1
                self.p2 = p2
                self.p3 = p4
                self.p4 = p3
            else:
                raise ValueError("Points must form a parallelogram")

Here we define the ``inside()`` function. If our test point ``x`` is out of plane, we immediately know it isn't inside the region. We then break the parallelogram into four triangles by connecting the vertices with ``x``. ``x`` is inside the region if and only if the areas of the triangles sum up to equal the area of the parallelogram.

::

        def inside(self, x, on_boundary):
            p1c = self.p1.array()
            p2c = self.p2.array()
            p3c = self.p3.array()
            p4c = self.p4.array()
            
            v12 = p2c - p1c
            v13 = p3c - p1c
            v14 = p4c - p1c
            v1x - x - p1c
            
            # If point is close to p1, return true
            if near(np.dot(v1x,v1x), 0.0, tol):
                return True
            
            N = np.cross(v12, v13)
            
            d = np.dot(v1x, N)/np.sqrt(np.dot(v1x, v1x))
            
            # If point is not near plane, return false
            if not near(d, 0.0, tol):
                return False
            
            A = np.array([x-p1c, x-p2c, x-p3c, x-p4c])
            s = 0.0
            s = s + mag_cross(A[0,:], A[1,:])
            s = s + mag_cross(A[1,:], A[2,:])
            s = s + mag_cross(A[2,:], A[3,:])
            s = s + mag_cross(A[3,:], A[0,:])
            s = s / 2.0
            
            return near(s, mag_cross(v12, v14), tol) and on_boundary

Next we define the methods for defining the region from four points and from a coordinate. The four points method validates the inputs by checking the points are coplanar and non-colinear.

::

        @classmethod
        def from_points(cls, p1:Point, p2:Point, p3:Point, p4:Point)
            colinear = False
            coords = np.array([p1.array(), p2.array(), p3.array(), p4.array()])
            for i in range(0,4):
                A = np.delete(coords, i, axis=0)
                v1 = A[1,:] - A[0,:]
                v2 = A[2,:] - A[0,:]
                
                if near(mag_cross(v1, v2), 0.0, tol):
                    colinear = True
            
            if colinear:
                raise ValueError("Points must be non-colinear")
            
            v12 = coords[1,:] - coords[0,:]
            v13 = coords[2,:] - coords[0,:]
            v14 = coords[3,:] - coords[0,:]
            
            coplanar = near(np.dot(v12, np.cross(v13, v14)), 0.0)
            
            if not coplanar:
                raise ValueError("Points must be coplanar")
            
            return cls(p1, p2, p3, p4)
        
        @classmethod
        def from_coord(cls, dimension:str, val:float):
            if dimension == 'x' or dimension == 'X':
                return cls(Point(val,-10000.0,-10000.0), Point(val,10000.0,-10000.0),\
                        Point(val,10000.0,10000.0), Point(val,-10000.0,10000.0))
            elif dimension == 'y' or dimension == 'Y':
                return cls(Point(-10000.0,val,-10000.0), Point(10000.0,val,-10000.0),\
                        Point(10000.0,val,10000.0), Point(-10000.0,val,10000.0))
            elif dimension == 'z' or dimension == 'Z':
                return cls(Point(-10000.0,-10000.0,val), Point(10000.0,-10000.0,val),\
                        Point(10000.0,10000.0,val), Point(-10000.0,10000.0,val))
            else:
                raise ValueError("Dimension must be 'x', 'y', or 'z'")
        
        def area(self):
            return mag_cross(self.p2.array()-self.p1.array(), self.p4.array()-self.p1.array())

-------------
Complete Code
-------------
The complete code follows and can also be downloaded :download:`here </../code/region_selector_3d.py>`.

.. literalinclude:: /../code/region_selector_3d.py
   :language: python