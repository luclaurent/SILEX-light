#############################################################################
#      Import libraries
#############################################################################
import string
import time
import numpy
import scipy.sparse
import scipy.sparse.linalg


# Choose between the Fortran or the Python librairie:

# import silex_lib_tet4_fortran as silex_lib_elt
from SILEXlight import silex_lib_tet4_python as silex_lib_elt
from SILEXlight import silex_lib_gmsh

# install mumps to solve the system more quickly
# import mumps

#############################################################################
print("SILEX CODE - computation of a with tet4 elements")
#############################################################################

tic = time.perf_counter()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName = 'piston-tet4'

# Output result file: define the name of the result file (*.msh)
ResultsFileName = 'Results_piston_tet4'

# choose the element type
eltype = 4

# choose geometry dimension
ndim = 3

# choose the results in the results file
flag_write_fields = 0

# read the mesh from gmsh
nodes = silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh', ndim)
elements, Idnodes = silex_lib_gmsh.ReadGmshElements(
    MeshFileName+'.msh', eltype, 5)

# read surfaces where to impose boundary conditions
elementsS1, IdnodeS1 = silex_lib_gmsh.ReadGmshElements(
    MeshFileName+'.msh', 2, 1)
elementsS2, IdnodeS2 = silex_lib_gmsh.ReadGmshElements(
    MeshFileName+'.msh', 2, 2)
elementsS3, IdnodeS3 = silex_lib_gmsh.ReadGmshElements(
    MeshFileName+'.msh', 2, 3)
elementsS4, IdnodeS4 = silex_lib_gmsh.ReadGmshElements(
    MeshFileName+'.msh', 2, 4)

# write the surface mesh in a gmsh-format file to verify if its correct
# silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,4)
# silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,2)
# silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
# silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)
# silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,2)

# Define material
Young = 200000.0
nu = 0.3

# Boundary conditions
IdNodesFixed_x = IdnodeS3
IdNodesFixed_y = IdnodeS1
IdNodesFixed_z = IdnodeS2

# compute external forces from pressure
press = -10  # MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [0.0, 0.0, 0.0]
F = silex_lib_elt.forceonsurface(nodes, elementsS4, press, direction)


toc = time.perf_counter()
print("Time for the user part:", toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof = nnodes*ndim
nelem = elements.shape[0]
print("Number of nodes:", nnodes)
print("Number of elements:", nelem)

# define fixed dof
Fixed_Dofs = numpy.hstack(
    [(IdNodesFixed_x-1)*3, (IdNodesFixed_y-1)*3+1, (IdNodesFixed_z-1)*3+2])

# define free dof
SolvedDofs = numpy.setdiff1d(range(ndof), Fixed_Dofs)

# initialize displacement vector
Q = numpy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.perf_counter()
tic = time.perf_counter()
Ik, Jk, Vk = silex_lib_elt.stiffnessmatrix(nodes, elements, [Young, nu])
toc = time.perf_counter()

K = scipy.sparse.csr_matrix((Vk, (Ik, Jk)), shape=(ndof, ndof), dtype=float)
print("Time to compute the stiffness matrix:", toc-tic)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.perf_counter()
Q[SolvedDofs] = scipy.sparse.linalg.spsolve(
    K[SolvedDofs, :][:, SolvedDofs], F[SolvedDofs])
# install mumps to improve the computational time
# import mumps
# Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.perf_counter()
print("Time to solve the problem:", toc-tic)

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
tic = time.perf_counter()

SigmaElem, SigmaNodes, EpsilonElem, EpsilonNodes, ErrorElem, ErrorGlobal = silex_lib_elt.compute_stress_strain_error(
    nodes, elements, [Young, nu], Q)

toc = time.perf_counter()
print("Time to compute stress and error:", toc-tic)
print("The global error is:", ErrorGlobal)
print("Total time for the computational part:", toc-tic0)

#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.perf_counter()

# displacement written on 3 columns:
disp = numpy.zeros((nnodes, ndim))
disp[range(nnodes), 0] = Q[list(range(0, ndof, 3))]
disp[range(nnodes), 1] = Q[list(range(1, ndof, 3))]
disp[range(nnodes), 2] = Q[list(range(2, ndof, 3))]

# external load written on 3 columns:
load = numpy.zeros((nnodes, ndim))
load[range(nnodes), 0] = F[list(range(0, ndof, 3))]
load[range(nnodes), 1] = F[list(range(1, ndof, 3))]
load[range(nnodes), 2] = F[list(range(2, ndof, 3))]

if flag_write_fields == 0:
    fields_to_write = [[disp, 'nodal', ndim, 'displacement'],
                       [SigmaElem[range(nelem), [6]],
                        'elemental', 1, 'Sigma V.M.'],
                       [SigmaNodes[range(nnodes), [6]], 'nodal',
                        1, 'Smooth Sigma V.M.'],
                       [ErrorElem, 'elemental', 1, 'error'],
                       ]

if flag_write_fields == 1:
    fields_to_write = [[disp, 'nodal', ndim, 'displacement'],
                       [load, 'nodal', ndim, 'Force'],
                       [SigmaElem[range(nelem), [0]],
                        'elemental', 1, 'Sigma 11'],
                       [SigmaElem[range(nelem), [1]],
                        'elemental', 1, 'Sigma 22'],
                       [SigmaElem[range(nelem), [2]],
                        'elemental', 1, 'Sigma 33'],
                       [SigmaElem[range(nelem), [3]],
                        'elemental', 1, 'Sigma 23'],
                       [SigmaElem[range(nelem), [4]],
                        'elemental', 1, 'Sigma 13'],
                       [SigmaElem[range(nelem), [5]],
                        'elemental', 1, 'Sigma 12'],
                       [SigmaElem[range(nelem), [6]],
                        'elemental', 1, 'Sigma V.M.'],
                       [SigmaNodes[range(nnodes), [0]], 'nodal',
                        1, 'Smooth Sigma 11'],
                       [SigmaNodes[range(nnodes), [1]], 'nodal',
                        1, 'Smooth Sigma 22'],
                       [SigmaNodes[range(nnodes), [2]], 'nodal',
                        1, 'Smooth Sigma 33'],
                       [SigmaNodes[range(nnodes), [3]], 'nodal',
                        1, 'Smooth Sigma 23'],
                       [SigmaNodes[range(nnodes), [4]], 'nodal',
                        1, 'Smooth Sigma 13'],
                       [SigmaNodes[range(nnodes), [5]], 'nodal',
                        1, 'Smooth Sigma 12'],
                       [SigmaNodes[range(nnodes), [6]], 'nodal',
                        1, 'Smooth Sigma V.M.'],
                       [EpsilonElem[range(nelem), [0]],
                        'elemental', 1, 'Epsilon 11'],
                       [EpsilonElem[range(nelem), [1]],
                        'elemental', 1, 'Epsilon 22'],
                       [EpsilonElem[range(nelem), [2]],
                        'elemental', 1, 'Epsilon 33'],
                       [EpsilonElem[range(nelem), [3]]/2.0,
                        'elemental', 1, 'Epsilon 23'],
                       [EpsilonElem[range(nelem), [4]]/2.0,
                        'elemental', 1, 'Epsilon 13'],
                       [EpsilonElem[range(nelem), [5]]/2.0,
                        'elemental', 1, 'Epsilon 12'],
                       [EpsilonNodes[range(nnodes), [0]],
                        'nodal', 1, 'Smooth Epsilon 11'],
                       [EpsilonNodes[range(nnodes), [1]],
                        'nodal', 1, 'Smooth Epsilon 22'],
                       [EpsilonNodes[range(nnodes), [2]],
                        'nodal', 1, 'Smooth Epsilon 33'],
                       [EpsilonNodes[range(nnodes), [3]]/2.0, 'nodal',
                        1, 'Smooth Epsilon 23'],
                       [EpsilonNodes[range(nnodes), [4]]/2.0, 'nodal',
                        1, 'Smooth Epsilon 13'],
                       [EpsilonNodes[range(nnodes), [5]]/2.0, 'nodal',
                        1, 'Smooth Epsilon 12'],
                       [ErrorElem, 'elemental', 1, 'error'],
                       ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName, nodes,
                            elements, eltype, fields_to_write)

toc = time.perf_counter()
print("Time to write results:", toc-tic)
print("----- END -----")
