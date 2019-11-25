from __future__ import print_function
import region_selector_2d as rs
from dolfin import *

length = 16.0 # [mm]
thickness = 3.0 # [mm]
resolution = 2 # [Nodes/mm]

resX = int(resolution*length) # Num nodes in x axis
resY = int(resolution * thickness) # Num nodes in y axis

mesh = RectangleMesh(Point(0.0,0.0), Point(length, thickness), resX, resY)

fixedRegion = rs.GetLinearBoundary.from_points(Point(0.0,0.0), Point(0.0,thickness))
loadRegion = rs.GetLinearBoundary.from_points(Point(length-1.0, thickness), Point(length, thickness))
   
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
fixedRegion.mark(boundaries, 1)
loadRegion.mark(boundaries, 2)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

folder_name = './2d_cantilever_results'
    
boundaryFile = File('%s/boundaries.pvd' % folder_name)
boundaryfile << boundaries

E = 70000.0
nu = 0.3
mu = E / (2.0 * (1.0 + nu))
lmbda = E*nu / ((1.0 + nu) * (1.0-2.0*nu))

load = Constant((0.0,-1000.0)) # [N/mm]

def eps(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(eps(u)) * Identity(mesh.topology.dim()) + 2.0*mu*eps

V = VectorFunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)
u = Function(V, name="Displacement")
v = TestFunction(V)

a = inner(sigma(du), eps(v))*dx
L = dot(load,v)*ds(2)

bc = DirichletBC(V, Constant((0.0,0.0)), fixedRegion)

solve(a == L, u, bc)

u.rename("Displacement", "Displacement")
xdmf_file = XDMFFile('%s/results.xdmf' % folder_name)
xdmf_file.write(u,1.0)
