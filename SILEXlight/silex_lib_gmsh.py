"""Utilities for reading and writing Gmsh meshes and result fields."""

def WriteResults(filename, nodes, elements, elttype, fields=[]):
    """Write a mesh and optional fields to a Gmsh ``.msh`` file.

    Parameters
    ----------
    filename : str or path-like
        Output file name without or with ``.msh`` extension.
    nodes : numpy.ndarray
        Node coordinates, shape ``(n_nodes, 2)`` or ``(n_nodes, 3)``.
    elements : numpy.ndarray
        Element connectivity array.
    elttype : int
        Gmsh element type id.

        Supported ids in this writer:
        - ``1``: 2-node line
        - ``2``: 3-node triangle
        - ``3``: 4-node quadrangle
        - ``4``: 4-node tetrahedron
        - ``5``: 8-node hexahedron
        - ``8``: 3-node second-order line
        - ``9``: 6-node second-order triangle
        - ``11``: 10-node second-order tetrahedron

        Common Gmsh ids reference:
        - ``1``: 2-node line
        - ``2``: 3-node triangle
        - ``3``: 4-node quadrangle
        - ``4``: 4-node tetrahedron
        - ``5``: 8-node hexahedron
        - ``6``: 6-node prism
        - ``7``: 5-node pyramid
        - ``8``: 3-node second-order line
        - ``9``: 6-node second-order triangle
        - ``10``: 9-node second-order quadrangle
        - ``11``: 10-node second-order tetrahedron
        - ``12``: 27-node second-order hexahedron
        - ``13``: 18-node second-order prism
        - ``14``: 14-node second-order pyramid
        - ``15``: 1-node point
        - ``16``: 8-node second-order quadrangle
        - ``17``: 20-node second-order hexahedron
        - ``18``: 15-node second-order prism
        - ``19``: 13-node second-order pyramid
    fields : list, optional
        Field descriptors as ``[values, 'nodal'|'elemental', ncomp, name]``.
    """
    
    if not isinstance(filename, str):
        filename = str(filename)
    
    if not filename.endswith('.msh'):
        filename = filename + '.msh'

    f = open(filename, 'w')

    nnodes = nodes.shape[0]
    ndof = nnodes*2
    nelem = elements.shape[0]

    f.write('$MeshFormat\n')
    f.write('2.2 0 8\n')
    f.write('$EndMeshFormat\n')
    f.write('$Nodes\n')
    f.write('%d\n' % (nnodes))

    ###############
    # write nodes #
    ###############
    # (2d)
    if nodes.shape[1] == 2:
        for i in range(nnodes):
            f.write('%d %s %s 0.0\n' % (i+1,
                                        nodes[i, 0],
                                        nodes[i, 1]))

    # (3d)
    if nodes.shape[1] == 3:
        for i in range(nnodes):
            f.write('%d %s %s %s\n' % (i+1,
                                       nodes[i, 0],
                                       nodes[i, 1],
                                       nodes[i, 2]))

    f.write('$EndNodes\n')

    ##################
    # write elements #
    ##################
    f.write('$Elements\n')
    f.write('%d\n' % (nelem))
    # 2-node truss elements
    if elttype == 1:
        for e in range(nelem):
            f.write('%d 1 3 1 1 0 %d %d\n' % (e+1,
                                              elements[e, 0],
                                              elements[e, 1]))

    # 3-node triangle
    if elttype == 2:
        for e in range(nelem):
            f.write('%d 2 3 1 1 0 %d %d %d\n' % (e+1,
                                                 elements[e, 0],
                                                 elements[e, 1],
                                                 elements[e, 2]))

    # 4-node quadrangle
    if elttype == 3:
        for e in range(nelem):
            f.write('%d 3 3 1 1 0 %d %d %d %d\n' % (e+1,
                                                    elements[e, 0],
                                                    elements[e, 1],
                                                    elements[e, 2],
                                                    elements[e, 3]))
    # 4-node tetrahedron
    if elttype == 4:
        for e in range(nelem):
            f.write('%d 4 3 1 1 0 %d %d %d %d\n' % (e+1,
                                                    elements[e, 0],
                                                    elements[e, 1],
                                                    elements[e, 2],
                                                    elements[e, 3]))

    # 8-node hexahedron
    if elttype == 5:
        for e in range(nelem):
            f.write('%d 5 3 1 1 0 %d %d %d %d %d %d %d %d\n' % (e+1,
                                                                elements[e, 0],
                                                                elements[e, 1],
                                                                elements[e, 2],
                                                                elements[e, 3],
                                                                elements[e, 4],
                                                                elements[e, 5],
                                                                elements[e, 6],
                                                                elements[e, 7]))
    # 3-node line
    if elttype == 8:
        for e in range(nelem):
            f.write('%d 8 3 1 1 0 %d %d %d \n' % (e+1,
                                                  elements[e, 0],
                                                  elements[e, 1],
                                                  elements[e, 2]))
    # 6-node triangle
    if elttype == 9:
        for e in range(nelem):
            f.write('%d 9 3 1 1 0 %d %d %d %d %d %d\n' % (e+1,
                                                          elements[e, 0],
                                                          elements[e, 1],
                                                          elements[e, 2],
                                                          elements[e, 3],
                                                          elements[e, 4],
                                                          elements[e, 5]))

    # 10-node tetrahedron
    if elttype == 11:
        for e in range(nelem):
            f.write('%d 11 3 1 1 0 %d %d %d %d %d %d %d %d %d %d\n' % (e+1,
                                                                       elements[e, 0],
                                                                       elements[e, 1],
                                                                       elements[e, 2],
                                                                       elements[e, 3],
                                                                       elements[e, 4],
                                                                       elements[e, 5],
                                                                       elements[e, 6],
                                                                       elements[e, 7],
                                                                       elements[e, 8],
                                                                       elements[e, 9]))

    f.write('$EndElements\n')

    #################
    # write results #
    #################
    for p in range(len(fields)):
        values = fields[p][0]
        valuestype = fields[p][1]
        nbpernode = fields[p][2]
        name = fields[p][3]
        # f.write('View[%d].Visible = 0;\n' % p)
        if valuestype == 'nodal':
            f.write('$NodeData\n')
            f.write('1\n')
            f.write('"%s"\n' % name)
            f.write('1\n')
            f.write('0.0\n')
            f.write('3\n')
            f.write('0\n')
            if nbpernode == 1:
                f.write('1\n')
            if nbpernode == 2:
                f.write('3\n')
            if nbpernode == 3:
                f.write('3\n')
            f.write('%d\n' % (nnodes))
            if nbpernode == 1:
                for i in range(nnodes):
                    f.write('%d %s\n' % (i+1, values[i]))

            if nbpernode == 2:
                for i in range(nnodes):
                    f.write('%d %s %s 0.0\n' %
                            (i+1, values[i, 0], values[i, 1]))

            if nbpernode == 3:
                for i in range(nnodes):
                    f.write('%d %s %s %s\n' %
                            (i+1, values[i, 0], values[i, 1], values[i, 2]))

            f.write('$EndNodeData\n')

        if valuestype == 'elemental':
            f.write('$ElementData\n')
            f.write('1\n')
            f.write('"%s"\n' % name)
            f.write('1\n')
            f.write('0.0\n')
            f.write('3\n')
            f.write('0\n')
            if nbpernode == 1:
                f.write('1\n')
            if nbpernode == 2:
                f.write('3\n')
            if nbpernode == 3:
                f.write('3\n')
            f.write('%d\n' % (nelem))

            if nbpernode == 1:
                for e in range(nelem):
                    f.write('%d %s\n' % (e+1, values[e]))

            if nbpernode == 2:
                for e in range(nelem):
                    f.write('%d %s %s\n' % (e+1, values[e, 0], values[e, 1]))

            if nbpernode == 3:
                for e in range(nelem):
                    f.write('%d %s %s %s\n' %
                            (e+1, values[e, 0], values[e, 1], values[e, 2]))

            f.write('$EndElementData\n')

    f.close()
    return

##############################################################
# FUNCTION : read nodes from GMSH output file
##############################################################


def ReadGmshNodes(file, nbcoord):
    """Read node coordinates from a Gmsh ``.msh`` file.

    Parameters
    ----------
    file : str or path-like
        Input mesh file.
    nbcoord : int
        Number of coordinates per node (2 or 3).
    """
    import string
    import numpy
    
    if not isinstance(file, str):
        file = str(file)
    
    f = open(file, 'r')
    i = 0
    j = 0
    # nodes=[]
    goodline = 0
    nbnodes = 0
    for l in f:
        # s=string.split(l)
        s = l.split()
        if (goodline == 1):
            if (nbnodes == 0):
                nbnodes = int(s[0])
                nodes = numpy.zeros((nbnodes, nbcoord))
            else:
                if (i < nbnodes):
                    if nbcoord == 3:
                        nodes[i, 0] = float(s[1])
                        nodes[i, 1] = float(s[2])
                        nodes[i, 2] = float(s[3])
                        i = i+1

                    if nbcoord == 2:
                        nodes[i, 0] = float(s[1])
                        nodes[i, 1] = float(s[2])
                        i = i+1

        if (s[0] == '$Nodes'):
            goodline = 1

    f.close()
    return nodes

##############################################################
# FUNCTION : read elements from GMSH output file
##############################################################


def ReadGmshElements(file, elmttype, prop):
    """Read elements of a given type and physical property from Gmsh file.

    Parameters
    ----------
    file : str or path-like
        Input mesh file.
    elmttype : int
        Gmsh element type id to extract.

        Currently handled by this reader:
        - ``1``: 2-node line
        - ``2``: 3-node triangle
        - ``3``: 4-node quadrangle
        - ``4``: 4-node tetrahedron
        - ``5``: 8-node hexahedron
        - ``8``: 3-node second-order line
        - ``9``: 6-node second-order triangle
        - ``11``: 10-node second-order tetrahedron
    prop : int
        Physical tag to filter elements.

    Returns
    -------
    tuple
        ``(elements, unique_node_ids)``.
    """
    import string
    import numpy
    f = open(file, 'r')
    i = -1
    j = 0
    elements = []
    goodline = 0
    if (elmttype == 11):
        nbnodes = 10
    if (elmttype == 9):
        nbnodes = 6
    if (elmttype == 8):
        nbnodes = 3
    if (elmttype == 5):
        nbnodes = 8
    if (elmttype == 4):
        nbnodes = 4
    if (elmttype == 2):
        nbnodes = 3
    if (elmttype == 3):
        nbnodes = 4
    if (elmttype == 1):
        nbnodes = 2

    for l in f:
        # s=string.split(l)
        s = l.split()
        if (goodline == 1):
            if (i == -1):
                nbelem = int(s[0])
                i = 0
                elements = numpy.zeros((nbelem, nbnodes), dtype=int)
            else:
                if (i < nbelem):
                    if (int(s[3]) == prop):
                        if (int(s[1]) == elmttype):
                            for k in range(nbnodes):
                                elements[j, k] = int(s[5+k])
                            j = j+1
                    i = i+1

        if (s[0] == '$Elements'):
            goodline = 1

    f.close()
    return elements[list(range(j)), :][:, list(range(nbnodes))], numpy.unique(elements[list(range(j)), :][:, list(range(nbnodes))])
