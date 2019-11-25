from __future__ import print_function
import region_selector_3d as rs
from dolfin import *

length = 20.0 # [mm]
width = 2.0 # [mm]
thickness = 1.0 # [mm]
resolution = 2 # [Nodes/mm]

resX = int(resolution * length) # Num nodes in x axis
resY = int(resolution * width) # Num nodes in y axis
resZ = int(resolution * thickness) # Num nodes in z axis

mesh = BoxMesh(Point(0.0,0.0,0.0), Point(length, width, thickness), resX, resY, resZ)

fixedRegion = rs.GetPlanarBoundary.from_coord('x', 0.0)
loadRegion = rs.GetPlanarBoundary.from_points(\
    Point(length, 0.0, thickness),\
    Point(length, width, thickness),\
    Point(length - 1.0, width, thickness),\
    Point(length - 1.0, 0.0, thickness))

boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
fixedRegion.mark(boundaries, 1)
loadRegion.mark(boundaries, 2)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

folder_name = './3d_cantilever_results'
    
boundaryfile = File('%s/boundaries.pvd' % folder_name)
boundaryfile << boundaries

E = 68900.0 # [N/mm^2]
nu = 0.33 # Poisson's ratio
mu = E / (2.0 * (1.0 + nu))
lmbda = E*nu / ((1.0 + nu) * (1.0-2.0*nu))

load = Constant((0.0,0.0,-5.0)) # [N]
# Convert to N/mm
loadArea = loadRegion.area()
scaledLoad = load / loadArea

def eps(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(eps(u)) * Identity(mesh.topology().dim()) + 2.0*mu*eps(u)

V = VectorFunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)
u = Function(V, name="Displacement")
v = TestFunction(V)

a = inner(sigma(du), eps(v))*dx
L = dot(scaledLoad,v)*ds(2)

bc = DirichletBC(V, Constant((0.0,0.0,0.0)), fixedRegion)

solve(a == L, u, bc)

u.rename("Displacement", "Displacement")
xdmf_file = XDMFFile('%s/results.xdmf' % folder_name)
xdmf_file.write(u,1.0)
