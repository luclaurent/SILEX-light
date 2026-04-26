Complete Examples (All Libraries)
=================================

This page gives practical examples for every SILEXlight library:

- ``SILEXlight.silex_lib_gmsh``
- ``SILEXlight.silex_lib_tri3_python``
- ``SILEXlight.silex_lib_tet4_python``
- ``SILEXlight.silex_lib_tri3_fortran``
- ``SILEXlight.silex_lib_tet4_fortran``

Conventions used in all examples
--------------------------------

- Node numbering in element connectivity starts at ``1``.
- ``nodes`` shapes:
  - 2D: ``(n_nodes, 2)``
  - 3D: ``(n_nodes, 3)``
- ``elements`` shapes:
  - TRI3: ``(n_elem, 3)``
  - TET4: ``(n_elem, 4)``
- Material arrays:
  - TRI3 plane stress: ``[Young, nu, thickness]``
  - TET4 3D isotropic: ``[Young, nu]``

Gmsh IO Library
---------------

Write a minimal TRI3 mesh and nodal displacement field to ``.msh``:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_gmsh as gmsh

   nodes = np.array([
       [0.0, 0.0],
       [1.0, 0.0],
       [0.0, 1.0],
   ])
   elements = np.array([[1, 2, 3]], dtype=int)

   # Example nodal displacement Ux, Uy for each node
   U = np.array([
       [0.0, 0.0],
       [1e-5, 0.0],
       [0.0, 1e-5],
   ])

   fields = [
       [U, 'nodal', 2, 'Displacement'],
   ]

   gmsh.WriteResults('tri3_demo', nodes, elements, 2, fields)

Read back nodes and elements from a gmsh file:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_gmsh as gmsh

   nodes2d = gmsh.ReadGmshNodes('tri3_demo.msh', nbcoord=2)
   tri3, used_node_ids = gmsh.ReadGmshElements('tri3_demo.msh', elmttype=2, prop=1)

   print(nodes2d.shape)       # (n_nodes, 2)
   print(tri3.shape)          # (n_triangles, 3)
   print(used_node_ids[:10])

TRI3 Python Library
-------------------

Compute element assembly triplets (I, J, V) for stiffness matrix:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tri3_python as tri3

   nodes = np.array([
       [0.0, 0.0],
       [1.0, 0.0],
       [0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3]], dtype=int)
   material = np.array([210e9, 0.30, 0.01], dtype=float)

   Ik, Jk, Vk = tri3.stiffnessmatrix(nodes, elements, material)
   print(Ik.shape, Jk.shape, Vk.shape)  # (36,) each for one TRI3 element

Compute stress, strain, and error indicator from a displacement vector:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tri3_python as tri3

   nodes = np.array([
       [0.0, 0.0],
       [1.0, 0.0],
       [0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3]], dtype=int)
   material = np.array([210e9, 0.30, 0.01], dtype=float)

   # 2 dof per node -> 6 dof total
   QQ = np.array([0.0, 0.0, 1e-5, 0.0, 0.0, 1e-5], dtype=float)

   Sigma, sig_smooth, Epsilon, EpsilonNodes, ErrElem, ErrGlob = tri3.compute_stress_strain_error(
       nodes, elements, material, QQ
   )

   print(Sigma.shape)         # (n_elem, 4)  : [sx, sy, txy, von_mises]
   print(sig_smooth.shape)    # (n_nodes, 4)
   print(Epsilon.shape)       # (n_elem, 3)
   print(EpsilonNodes.shape)  # (n_nodes, 3)
   print(float(ErrGlob))

Apply distributed force on line elements:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tri3_python as tri3

   nodes = np.array([
       [0.0, 0.0],
       [1.0, 0.0],
       [0.0, 1.0],
   ], dtype=float)

   # Load applied on edge (node 1 -> node 2)
   edge_elements = np.array([[1, 2]], dtype=int)

   # [fx(pt1), fy(pt1), fx(pt2), fy(pt2)] in N/m
   fs = np.array([0.0, -1000.0, 0.0, -1000.0], dtype=float)

   # Coordinates of pt1 and pt2
   pts = np.array([0.0, 0.0, 1.0, 0.0], dtype=float)

   F = tri3.forceonline(nodes, edge_elements, fs, pts)
   print(F.shape)  # (2*n_nodes,)

TET4 Python Library
-------------------

Compute 3D stiffness assembly triplets:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tet4_python as tet4

   nodes = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.0, 1.0, 0.0],
       [0.0, 0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3, 4]], dtype=int)
   material = np.array([210e9, 0.30], dtype=float)

   Ik, Jk, Vk = tet4.stiffnessmatrix(nodes, elements, material)
   print(Ik.shape, Jk.shape, Vk.shape)  # (144,) each for one TET4 element

Compute consistent mass assembly triplets:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tet4_python as tet4

   nodes = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.0, 1.0, 0.0],
       [0.0, 0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3, 4]], dtype=int)
   rho = 7800.0

   Ik, Jk, Vk = tet4.massmatrix(nodes, elements, rho)
   print(Ik.shape, Jk.shape, Vk.shape)  # (144,) each

Apply pressure (or directional traction) on triangular surface elements:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tet4_python as tet4

   nodes = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.0, 1.0, 0.0],
       [0.0, 0.0, 1.0],
   ], dtype=float)

   # Surface triangle on nodes (1,2,3)
   surface_tri = np.array([[1, 2, 3]], dtype=int)

   press = 2.0e6

   # direction=[0,0,0] -> pressure normal to face
   Fp = tet4.forceonsurface(nodes, surface_tri, press, direction=[0.0, 0.0, 0.0])
   print(Fp.shape)  # (3*n_nodes,)

Compute 3D stress/strain and global error estimator:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tet4_python as tet4

   nodes = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.0, 1.0, 0.0],
       [0.0, 0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3, 4]], dtype=int)
   material = np.array([210e9, 0.30], dtype=float)

   # 3 dof per node -> 12 dof total
   QQ = np.array([
       0.0, 0.0, 0.0,
       1e-5, 0.0, 0.0,
       0.0, 1e-5, 0.0,
       0.0, 0.0, 1e-5,
   ], dtype=float)

   sigma, sig_smooth, eps_elem, eps_nodes, err_elem, err_glob = tet4.compute_stress_strain_error(
       nodes, elements, material, QQ
   )

   print(sigma.shape)       # (n_elem, 7)
   print(sig_smooth.shape)  # (n_nodes, 7)
   print(eps_elem.shape)    # (n_elem, 6)
   print(eps_nodes.shape)   # (n_nodes, 6)
   print(float(err_glob))

Fortran Libraries (Optional, Compiled)
--------------------------------------

The Fortran modules provide the same core finite-element routines with a compiled backend:

- ``SILEXlight.silex_lib_tri3_fortran``
- ``SILEXlight.silex_lib_tet4_fortran``

If the package is installed with the Fortran extension built, usage is analogous.

TRI3 Fortran example:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tri3_fortran as tri3f

   nodes = np.array([
       [0.0, 0.0],
       [1.0, 0.0],
       [0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3]], dtype=int)
   material = np.array([210e9, 0.30, 0.01], dtype=float)

   Ik, Jk, Vk = tri3f.stiffnessmatrix(nodes, elements, material)
   print(Ik.shape, Jk.shape, Vk.shape)

TET4 Fortran example:

.. code-block:: python

   import numpy as np
   from SILEXlight import silex_lib_tet4_fortran as tet4f

   nodes = np.array([
       [0.0, 0.0, 0.0],
       [1.0, 0.0, 0.0],
       [0.0, 1.0, 0.0],
       [0.0, 0.0, 1.0],
   ], dtype=float)
   elements = np.array([[1, 2, 3, 4]], dtype=int)
   material = np.array([210e9, 0.30], dtype=float)

   Ik, Jk, Vk = tet4f.stiffnessmatrix(nodes, elements, material)
   print(Ik.shape, Jk.shape, Vk.shape)

Notes on Solver Integration
---------------------------

All ``stiffnessmatrix`` and ``massmatrix`` functions return assembly triplets ``(Ik, Jk, Vk)``.
A common sparse assembly pattern with SciPy is:

.. code-block:: python

   import scipy.sparse as sp

   ndof = 2 * n_nodes   # TRI3
   # ndof = 3 * n_nodes # TET4

   K = sp.coo_matrix((Vk, (Ik, Jk)), shape=(ndof, ndof)).tocsr()

This assembled matrix can then be used with your preferred linear solver.
