from dolfin import *
import numpy as np

# Select rectangular region
class GetRectangularRegion(SubDomain):

    # x1: x coordinate of one corner
    # y1: y coordinate of one corner
    # x2: x coordinate of opposite corner
    # y2: y coordinate of opposite corner
    def __init__(self, x1:float, y1:float, x2:float, y2:float):
        super().__init__()
        self.x1 = np.minimum(x1, x2)
        self.y1 = np.minimum(y1, y2)
        self.x2 = np.maximum(x1, x2)
        self.y2 = np.maximum(y1, y2)

    # p1: Point of one corner of rectangle
    # p2: Point of opposite corner of rectangle
    @classmethod
    def from_points(cls, p1:Point, p2:Point) -> 'GetRectangularRegion':
        return cls(p1.x(), p1.y(), p2.x(), p2.y())

    # x1: x coordinate of one corner
    # y1: y coordinate of one corner
    # x2: x coordinate of opposite corner
    # y2: y coordiante of opposite corner
    @classmethod
    def from_floats(cls, x1:float, y1:float, x2:float, y2:float) -> 'GetRectangularRegion':
        return cls(x1, y1, x2, y2)

    def inside(self,x,on_boundary):
        return between(x[0], (self.x1,self.x2)) and between(x[1], (self.y1,self.y2))

# Select circular region
class GetCircularRegion(SubDomain):
    # cx: float x coord of circle center
    # cy: float y coord of circle center
    # r: float radius of circle
    def __init__(self, cx:float, cy:float, r:float):
        super().__init__()
        self.cx = cx
        self.cy = cy
        self.r = r
    
    # c: Point center of circle
    # r: float radius of circle
    @classmethod
    def from_points(cls, c:Point, r:float) -> 'GetCircularRegion':
        return cls(c.x(), c.y(), r)

    # cx: float x coord of circle center
    # cy: float y coord of circle center
    # r: float radius of circle
    @classmethod
    def from_floats(cls, cx:float, cy:float, r:float) -> 'GetCircularRegion':
        return cls(cx, cy, r)

    def inside(self,x,on_boundary):
        return (x[0] - self.cx)**2.0 + (x[1] - self.cy)**2.0 <= self.r**2.0

# Select linear horizonta/vertical boundary region
class GetLinearBoundary(SubDomain):
    # coord: Constant coordinate of the boundary line
    # range1: Min value along line to select
    # range2: Max value along line to select
    # horizontal: True/False whether the line is horizontal or vertical
    def __init__(self, coord:float, range1:float, range2:float, horizontal:bool):
        super().__init__()
        self.coord = coord
        self.range1 = np.minimum(range1,range2)
        self.range2 = np.maximum(range1,range2)
        self.ishorizontal = horizontal

    # p1: Point on one end of line
    # p2: Point on other end of line
    @classmethod
    def from_points(cls, p1:Point, p2:Point) -> 'GetLinearBoundary':
        x1 = p1.x()
        y1 = p1.y()
        x2 = p2.x()
        y2 = p2.y()

        # If x doesn't change, line is vertical
        if near(x1, x2):
            return cls(x1, y1, y2, False)
        # If y doesn't change, line is horizontal
        elif near(y1, y2):
            return cls(y1, x1, x2, True)
        else:
            raise ValueError("Linear boundaries must be horizontal or vertical")
        
    # coord: Constant coordinate of the boundary line
    # range1: Min value along line to select
    # range2: Max value along line to select
    # horizontal: True/False whether the line is horizontal or vertical
    @classmethod
    def from_floats(cls, coord:float, range1:float, range2:float, horizontal:bool) -> 'GetLinearBoundary':
        return cls(coord, range1, range2, horizontal)

    def inside(self,x,on_boundary):
        if self.ishorizontal:
            return near(x[1], self.coord) and between(x[0], (self.range1, self.range2)) and on_boundary
        else:
            return near(x[0], self.coord) and between(x[1], (self.range1, self.range2)) and on_boundary

# Selects a single point
class SelectPoint(SubDomain):
    # p: Point to select
    def __init__(self, p:Point):
        super().__init__()
        self.p = p

    def inside(self, x, on_boundary):
        return near(x[0], self.p.x()) and near(x[1], self.p.y())
