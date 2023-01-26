#############################################################################
#      Import libraries
#############################################################################
import string; import time; import numpy
import scipy.sparse
import scipy.sparse.linalg


## Choose between the Fortran or the Python librairie:

#import silex_lib_tri3 as silex_lib_elt
from SILEXlight import silex_lib_tri3_python as silex_lib_elt
from SILEXlight import silex_lib_gmsh

#install mumps to solve the system more quickly
#import mumps

#############################################################################
print("SILEX CODE - compute a plate with a hole with tri3 elements")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.perf_counter()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='plate-hole-tri3'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_plate-hole-tri3'

# choose the element type
eltype=2

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,3)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S1_mesh',nodes,elementsS1,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S3_mesh',nodes,elementsS3,1)

# Define material
Young  = 200000.0
nu     = 0.3
thickness = 1.0

# Boundary conditions

# define fixed dof
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS2

F=silex_lib_elt.forceonline(nodes,elementsS3,[0.0,10.0,0.0,10.0],[0.0,200.0,100.0,200.0])

toc = time.perf_counter()
print("Time for the reading data part:",toc-tic)

tic0 = time.perf_counter()
#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*ndim
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

# define fixed dof
Fixed_Dofs = numpy.hstack([(IdNodesFixed_x-1)*2,(IdNodesFixed_y-1)*2+1])

# define free dof
SolvedDofs = numpy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=numpy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic = time.perf_counter()

Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness])

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

toc = time.perf_counter()
print("Time to compute the stiffness matrix:",toc-tic)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.perf_counter()
Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
# install mumps to improve the computational time
#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.perf_counter()
print("Time to solve the problem:",toc-tic)

#############################################################################
#       compute stress, smooth stress, strain and error
#############################################################################
tic = time.perf_counter()

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness],Q)

toc = time.perf_counter()
print("Time to compute stress and error:",toc-tic)
print("The global error is:",ErrorGlobal)


#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.perf_counter()

# displacement written on 2 columns:
disp=numpy.zeros((nnodes,2))
disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
disp[range(nnodes),1]=Q[list(range(1,ndof,2))]

load=numpy.zeros((nnodes,ndim))
load[range(nnodes),0]=F[list(range(0,ndof,2))]
load[range(nnodes),1]=F[list(range(1,ndof,2))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'Displacement'],
                  [load,'nodal',2,'Force'],
                  [SigmaElem[range(nelem),[0]],'elemental',1,'Sigma xx'],
                  [SigmaElem[range(nelem),[1]],'elemental',1,'Sigma yy'],
                  [SigmaElem[range(nelem),[2]],'elemental',1,'Sigma xy'],
                  [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                  [SigmaNodes[range(nnodes),[0]],'nodal',1,'Sigma xx Smooth'],
                  [SigmaNodes[range(nnodes),[1]],'nodal',1,'Sigma yy Smooth'],
                  [SigmaNodes[range(nnodes),[2]],'nodal',1,'Sigma xy Smooth'],
                  [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                  [EpsilonElem[range(nelem),[0]],'elemental',1,'Epsilon xx'],
                  [EpsilonElem[range(nelem),[1]],'elemental',1,'Epsilon yy'],
                  [EpsilonElem[range(nelem),[2]]/2.0,'elemental',1,'Epsilon xy'],
                  [EpsilonNodes[range(nnodes),[0]],'nodal',1,'Epsilon xx Smooth'],
                  [EpsilonNodes[range(nnodes),[1]],'nodal',1,'Epsilon yy Smooth'],
                  [EpsilonNodes[range(nnodes),[2]]/2.0,'nodal',1,'Epsilon xy Smooth'],
                  [ErrorElem,'elemental',1,'Error']
                  ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)


toc = time.perf_counter()
print("Time to write results:",toc-tic)
print("----- END -----")
