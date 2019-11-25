from dolfin import *
import numpy as np

tol = .000001

# Finds magnitude of cross product of two vectors
def mag_cross(v1, v2):
    v1xv2 = np.cross(v1, v2)
    return np.sqrt(np.dot(v1xv2,v1xv2))

# Like near() function, but for vectors
def vec_near(v1, v2, tol = DOLFIN_EPS):
    for i in range(0, v1.shape[0]):
        if not near(v1[i], v2[i], tol):
            return False
    return True

# Defines a rectangular cuboid region
class GetCuboidRegion(SubDomain):
    # Define using 4 points where p2, p3, p4 are distinct and share edges with p1
    def __init__(self, p1:Point, p2:Point, p3:Point, p4:Point):
        super().__init__()
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

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

# Defines a cylindrical region
class GetCylindricalRegion(SubDomain):
    # Define using 2 points and a radius
    # p1: Center of circle at one end of cylinder
    # p2: Center of circle at other end of cylinder
    # r: Radius of cylinder
    def __init__(self, p1:Point, p2:Point, r:float):
        super().__init__()
        self.p1 = p1
        self.p2 = p2
        self.r = r

    def inside(self, x, on_boundary):
        p1c = self.p1.array()
        p2c = self.p2.array()

        a = p2c - p1c # primary axis of cylinder
        xv = x - p1c

        xva = np.dot(xv, a) # x*a

        xvm = np.dot(xv,xv) # |x|^2
        am = np.dot(a,a) # |a|^2
        xam = (xva**2.0)/am # (x*a)^2/|a|^2
        dsq = xvm - xam # d^2 = |x|^2 - proj(x onto a)^2

        # d^2 <= r^2 && 0 <= proj(x onto a)/|a| <= |a| -> 0 <= x dot a <= |a|^2
        return between(dsq, (0.0,self.r**2.0)) and between(xva, (0.0,am))

# Defines parallelogram boundary in 3-space
class GetPlanarBoundary(SubDomain):
    def __init__(self, p1:Point, p2:Point, p3:Point, p4:Point):
        super().__init__()
        v12 = p2.array() - p1.array()
        v13 = p3.array() - p1.array()
        v14 = p4.array() - p1.array()

        self.p1 = p1
        
        # Use parallelogram law of vector addition to check what's the diagonal
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
        # End up with p1, p2, p3, p4 with p1 across from p3 and p2 across from p4

        

    def inside(self, x, on_boundary):
        p1c = self.p1.array()
        p2c = self.p2.array()
        p3c = self.p3.array()
        p4c = self.p4.array()

        v12 = p2c - p1c
        v13 = p3c - p1c
        v14 = p4c - p1c
        v1x = x - p1c
        
        # If point is close to p1, return true
        if near(np.dot(v1x,v1x),0.0, tol):
            return True

        N = np.cross(v12, v13)
        
        # Point-plane distance
        d = np.dot(v1x, N)/np.sqrt(np.dot(v1x,v1x))
        
        # If point is not near plane, return false
        if not near(d, 0.0, tol):
            return False

        # Check if point is in bounds.
        # Sum of areas of triangles between x and each pair of adjacent corners
        # should equal the area of the parallelogram.
        A = np.array([x-p1c, x-p2c, x-p3c, x-p4c])
        s = 0.0
        s = s + mag_cross(A[0,:], A[1,:])
        s = s + mag_cross(A[1,:], A[2,:])
        s = s + mag_cross(A[2,:], A[3,:])
        s = s + mag_cross(A[3,:], A[0,:])
        s = s / 2.0
 
        return near(s, mag_cross(v12, v14), tol) and on_boundary
    
    # Define a region using the four corners of the parallelogram.
    # These points must be:
    #   Non-colinear
    #   Coplanar
    #   Form a parallelogram
    @classmethod
    def from_points(cls, p1:Point, p2:Point, p3:Point, p4:Point):
        colinear = False
        coords = np.array([p1.array(),p2.array(),p3.array(),p4.array()]);
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
    
    # Define a region that is parallel to the x, y, or z planes
    # at a specific x, y, or z coordinate.
    # Assumes the shape is contained within the +/-10000 square.
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
    
    # Gets area of boundary of parallelogram.
    # DOES NOT necessarily get area of selected region on mesh
    def area(self):
        return mag_cross(self.p2.array()-self.p1.array(), self.p4.array()-self.p1.array())

