Imported Mesh
=============

In this example we will import a mesh from a CAD file and use Fenics to model its displacement. This file was created in Autodesk Inventor and exported as an .igs file, which can be downloaded :download:`here <../mesh_1.igs>`. To convert it to a a Gmsh .msh file, run this line in the command line::

	gmsh -3 -clmax 1 -o mesh.msh mesh_1.igs

The ``-3`` flag creates a 3D mesh and the ``-clmax 1`` flag sets the max size of an element to 1 to help with the smoothing of the holes in the model. We then use another command to convert the .msh file to a .xml file for Dolfin::

	dolfin-convert mesh.msh mesh.xml

We can now import the mesh into a Fenics script and write it out to a .pvd file with no marked regions to preview the mesh geometry::

	from dolfin import *
	
	mesh = Mesh('mesh.xml')
	
	folder_name = './mesh1_results'
	
	regions = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
	regions.set_all(0)
	
	regionfile = File('%s/regions.pvd' % folder_name)
	regionfile << regions

We get this when we open the .pvd file in ParaView.

.. image:: mesh_1/generated_mesh.png

This piece is designed to be made of aluminum and hold up a load of up to :math:`300\text{N}` on the top platform while deflecting less than :math:`1\text{mm}`. Let's use Fenics to see if it meets the specifications.

First, let's define the fixed regions and loaded regions in the mesh::

	from dolfin import *
	import numpy as np
	import region_selector_3d as rs
	
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

We can now view the marked faces in ParaView:

.. image:: mesh_1/marked_mesh.png