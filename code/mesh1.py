from dolfin import *
import numpy as np
import region_selector_3d as rs

E = 69800.0 # [N/mm^2]
nu = 0.33
mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu))

appliedLoad = Constant((0.0,-300.0,0.0)) # [N]

mesh = Mesh('mesh.xml')
folder_name = './mesh1_results'

fixedPin1 = rs.GetCylindricalRegion(Point(0.0,0.0,5.0), Point(0.0,0.0,-5.0), 1.6)
fixedPin2 = rs.GetCylindricalRegion(Point(25.0,0.0,5.0), Point(25.0,0.0,-5.0), 1.1)
loadFace = rs.GetPlanarBoundary.from_coord('y', 31.0)

regions = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
regions.set_all(0)
fixedPin1.mark(regions, 1)
fixedPin2.mark(regions, 1)
loadFace.mark(regions, 2)

regionfile = File('%s/regions.pvd' % folder_name)
regionfile << regions

ds = Measure('ds', domain = mesh, subdomain_data = regions)

loadArea = assemble(Constant(1.0)*ds(2)) # Note: CAN'T use loadFace.area()
scaledLoad = appliedLoad / loadArea # [N/mm]

def eps(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(eps(u)) * Identity(mesh.topology().dim()) + 2.0*mu*eps(u)

V = VectorFunctionSpace(mesh, "Lagrange", 2)
du = TrialFunction(V)
u = Function(V, name = "Displacement")
v = TestFunction(V)

a = inner(sigma(du), eps(v))*dx
L = dot(scaledLoad, v)*ds(2)

bc1 = DirichletBC(V, Constant((0.0,0.0,0.0)), fixedPin1)
bc2 = DirichletBC(V, Constant((0.0,0.0,0.0)), fixedPin2)

solve(a == L, u, [bc1, bc2])

u.rename("Displacement", "Displacement")
xdmf_file = XDMFFile('%s/results.xdmf' % folder_name)
xdmf_file.write(u,1.0)
