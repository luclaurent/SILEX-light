
import numpy as np


def det33_ligne_de_un(a):
    return a[0, 1]*a[1, 2]-a[0, 2]*a[1, 1]-a[0, 0]*a[1, 2]+a[0, 0]*a[1, 1]+a[1, 0]*a[0, 2]-a[1, 0]*a[0, 1]


def stiffnessmatrix(nodes, elements, material):
    nbnodes = nodes.shape[0]
    nbelem = elements.shape[0]

    Young = material[0]
    nu = material[1]
    thickness = material[2]

    Ik = np.zeros(6*6*nbelem, dtype=int)
    Jk = np.zeros(6*6*nbelem, dtype=int)
    Vk = np.zeros(6*6*nbelem, dtype=float)
    a23 = np.zeros((2, 3), dtype=float)
    b = np.zeros(3, dtype=float)
    c = np.zeros(3, dtype=float)
    BB = np.zeros((3, 6), dtype=float)
    dofelem = np.zeros(6, dtype=int)
    dofx = np.zeros(3, dtype=int)
    dofy = np.zeros(3, dtype=int)
    idnodes = np.zeros(3, dtype=int)

    CC = np.array([[Young/(1-nu**2), nu*Young/(1-nu**2), 0],
                      [nu*Young/(1-nu**2), Young/(1-nu**2), 0],
                      [0, 0, Young/(2*(1+nu))]
                      ])

    p = 0
    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1
        dofx[:] = (idnodes)*2
        dofy[:] = (idnodes)*2+1

        dofelem[0] = dofx[0]
        dofelem[1] = dofx[1]
        dofelem[2] = dofx[2]
        dofelem[3] = dofy[0]
        dofelem[4] = dofy[1]
        dofelem[5] = dofy[2]

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]

        a23[0, 0] = X[0]
        a23[0, 1] = X[1]
        a23[0, 2] = X[2]
        a23[1, 0] = Y[0]
        a23[1, 1] = Y[1]
        a23[1, 2] = Y[2]
        det_of_sys = det33_ligne_de_un(a23)
        Area = abs(0.5*det_of_sys)

        b[0] = -(Y[2]-Y[1])/det_of_sys
        c[0] = +(X[2]-X[1])/det_of_sys

        b[1] = +(Y[2]-Y[0])/det_of_sys
        c[1] = -(X[2]-X[0])/det_of_sys

        b[2] = -(Y[1]-Y[0])/det_of_sys
        c[2] = +(X[1]-X[0])/det_of_sys

        BB[0, 0] = b[0]
        BB[0, 1] = b[1]
        BB[0, 2] = b[2]
        BB[0, 3] = 0.0
        BB[0, 4] = 0.0
        BB[0, 5] = 0.0
        BB[1, 0] = 0.0
        BB[1, 1] = 0.0
        BB[1, 2] = 0.0
        BB[1, 3] = c[0]
        BB[1, 4] = c[1]
        BB[1, 5] = c[2]
        BB[2, 0] = c[0]
        BB[2, 1] = c[1]
        BB[2, 2] = c[2]
        BB[2, 3] = b[0]
        BB[2, 4] = b[1]
        BB[2, 5] = b[2]

        ke = np.dot(BB.T, np.dot(CC, BB))*Area*thickness

        for i in range(6):
            #for j in range(12):
            Ik[p:(p+6)] = dofelem[i]
            Jk[p:(p+6)] = dofelem[:]
            Vk[p:(p+6)] = ke[i, :]
            p = p+6
            #Ik.append(dofelem[i])
            #Jk.append(dofelem[j])
            #Vk.append(ke[i,j])

    return Ik, Jk, Vk

############################################################################


def compute_stress_strain_error(nodes, elements, material, QQ):
    nbnodes = nodes.shape[0]
    nbelem = elements.shape[0]

    Young = material[0]
    nu = material[1]
    thickness = material[2]

    CC = np.array([[Young/(1-nu**2), nu*Young/(1-nu**2), 0],
                      [nu*Young/(1-nu**2), Young/(1-nu**2), 0],
                      [0, 0, Young/(2*(1+nu))]
                      ])
    CCinv = np.array([[1.0/Young, -nu/Young, 0.0],
                         [-nu/Young, 1.0/Young, 0.0],
                         [0.0, 0.0, (2*(1+nu))/Young]
                         ])

    rg = np.zeros(3, dtype=float)
    sg = np.zeros(3, dtype=float)
    tg = np.zeros(3, dtype=float)
    wg = np.zeros(3, dtype=float)

    nodalweight = np.zeros(nbnodes, dtype=float)
    sig_smooth = np.zeros((nbnodes, 4), dtype=float)
    a23 = np.zeros((2, 3), dtype=float)
    b = np.zeros(3, dtype=float)
    c = np.zeros(3, dtype=float)
    BB = np.zeros((3, 6), dtype=float)
    dofelem = np.zeros(6, dtype=int)
    dofx = np.zeros(3, dtype=int)
    dofy = np.zeros(3, dtype=int)
    idnodes = np.zeros(3, dtype=int)

    Sigma = np.zeros((nbelem, 4), dtype=float)
    Epsilon = np.zeros((nbelem, 3), dtype=float)
    EpsilonNodes = np.zeros((nbnodes, 3), dtype=float)
    ErrElem = np.zeros(nbelem, dtype=float)
    NormSigElt = np.zeros(nbelem, dtype=float)
    ErrGlob = 0.0
    NormSig = 0.0

    NN = np.zeros(3, dtype=float)
    sig_smooth_Gauss = np.zeros(3, dtype=float)
    sigdiff = np.zeros(3, dtype=float)

    # Quadrature points
    rg[0] = 2.0/3.0
    sg[0] = 1.0/6.0
    tg[0] = 1.0/6.0
    wg[0] = 0.5*1.0/3.0

    rg[1] = 1.0/6.0
    sg[1] = 2.0/3.0
    tg[1] = 1.0/6.0
    wg[1] = 0.5*1.0/3.0

    rg[2] = 1.0/6.0
    sg[2] = 1.0/6.0
    tg[2] = 2.0/3.0
    wg[2] = 0.5*1.0/3.0

    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1
        dofx[:] = (idnodes)*2
        dofy[:] = (idnodes)*2+1

        dofelem[0] = dofx[0]
        dofelem[1] = dofx[1]
        dofelem[2] = dofx[2]
        dofelem[3] = dofy[0]
        dofelem[4] = dofy[1]
        dofelem[5] = dofy[2]

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]
        Uelem = QQ[dofelem]

        a23[0, 0] = X[0]
        a23[0, 1] = X[1]
        a23[0, 2] = X[2]
        a23[1, 0] = Y[0]
        a23[1, 1] = Y[1]
        a23[1, 2] = Y[2]
        det_of_sys = det33_ligne_de_un(a23)
        Area = abs(0.5*det_of_sys)

        b[0] = -(Y[2]-Y[1])/det_of_sys
        c[0] = +(X[2]-X[1])/det_of_sys

        b[1] = +(Y[2]-Y[0])/det_of_sys
        c[1] = -(X[2]-X[0])/det_of_sys

        b[2] = -(Y[1]-Y[0])/det_of_sys
        c[2] = +(X[1]-X[0])/det_of_sys

        BB[0, 0] = b[0]
        BB[0, 1] = b[1]
        BB[0, 2] = b[2]
        BB[0, 3] = 0.0
        BB[0, 4] = 0.0
        BB[0, 5] = 0.0
        BB[1, 0] = 0.0
        BB[1, 1] = 0.0
        BB[1, 2] = 0.0
        BB[1, 3] = c[0]
        BB[1, 4] = c[1]
        BB[1, 5] = c[2]
        BB[2, 0] = c[0]
        BB[2, 1] = c[1]
        BB[2, 2] = c[2]
        BB[2, 3] = b[0]
        BB[2, 4] = b[1]
        BB[2, 5] = b[2]

        EpsilonElem = np.dot(BB, Uelem)
        SigmaElem = np.dot(CC, EpsilonElem)

        Epsilon[e, 0] = EpsilonElem[0]
        Epsilon[e, 1] = EpsilonElem[1]
        Epsilon[e, 2] = EpsilonElem[2]
        Sigma[e, 0] = SigmaElem[0]
        Sigma[e, 1] = SigmaElem[1]
        Sigma[e, 2] = SigmaElem[2]
        # Von Mises
        Sigma[e, 3] = np.sqrt(1.5*(SigmaElem[0]**2+SigmaElem[1] **
                                 2+2.0*SigmaElem[2]**2)-0.5*(SigmaElem[0]+SigmaElem[1])**2)

        for i in range(3):
            nodalweight[idnodes[i]] = nodalweight[idnodes[i]]+Area
            sig_smooth[idnodes[i], :] = sig_smooth[idnodes[i], :] + \
                Sigma[e, :]*Area

    # norm the smooth stress
    for i in range(nbnodes):
        for j in range(4):
            sig_smooth[i, j] = sig_smooth[i, j]/nodalweight[i]

        for j in range(3):
            EpsilonNodes[i, j] = CCinv[j, 0]*sig_smooth[i, 0] + \
                CCinv[j, 1]*sig_smooth[i, 1]+CCinv[j, 2]*sig_smooth[i, 2]

    # compute error
    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]

        a23[0, 0] = X[0]
        a23[0, 1] = X[1]
        a23[0, 2] = X[2]
        a23[1, 0] = Y[0]
        a23[1, 1] = Y[1]
        a23[1, 2] = Y[2]
        det_of_sys = det33_ligne_de_un(a23)
        Area = abs(0.5*det_of_sys)

        for g in range(3):
            for j in range(3):
                sig_smooth_Gauss[j] = sig_smooth[idnodes[0], j]*rg[g] + \
                    sig_smooth[idnodes[1], j]*sg[g] + \
                    sig_smooth[idnodes[2], j]*tg[g]

            sigdiff[:] = sig_smooth_Gauss[:]-Sigma[e, 0:3]
            ErrElem[e] = ErrElem[e] + \
                np.dot(sigdiff, np.dot(CCinv, sigdiff))*Area/0.5*wg[g]
            NormSigElt[e] = NormSigElt[e]+np.dot(
                sig_smooth_Gauss, np.dot(CCinv, sig_smooth_Gauss))*Area/0.5*wg[g]

    ErrElem = ErrElem/sum(NormSigElt)
    ErrGlob = np.sqrt(sum(ErrElem))

    return Sigma, sig_smooth, Epsilon, EpsilonNodes, ErrElem, ErrGlob

############################################################################


def forceonline(nodes, elements, fs, pts):
    # nodes: node coordinates
    # elements: 2-node line elements on which the force is applied
    # fs=[ surf. load on pt1 x-direc  , surf. load on pt1 y-direc , surf. load on pt2 x-direc  , surf. load on pt2 y-direc ]
    # fs : units = Newton per length
    # pts=[ pt 1 x , pt 1 y ,pt 2 x ,pt 2 y]
    nelem = elements.shape[0]
    rg = [-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)]
    wg = [1.0, 1.0]

    F = np.zeros(nodes.shape[0]*2)

    for e in range(nelem):

        idnodes = elements[e, :]

        dofx = (idnodes-1)*2
        dofy = (idnodes-1)*2+1

        xnodes = nodes[idnodes-1, 0]
        ynodes = nodes[idnodes-1, 1]
        lelem = np.sqrt((xnodes[1]-xnodes[0])**2+(ynodes[1]-ynodes[0])**2)

        forceelem = np.zeros((4, 1), dtype=float)

        for g in range(2):
            N = [(1-rg[g])*0.5, (1+rg[g])*0.5]
            x = N[0]*xnodes[0]+N[1]*xnodes[1]
            y = N[0]*ynodes[0]+N[1]*ynodes[1]

            l0 = np.sqrt((pts[0]-x)**2+(pts[1]-y)**2)
            l1 = np.sqrt((pts[2]-x)**2+(pts[3]-y)**2)

            forcex = (l1*fs[0]+l0*fs[2])/(l0+l1)
            forcey = (l1*fs[1]+l0*fs[3])/(l0+l1)

            Phi = np.array([[N[0], 0.0, N[1], 0.0],
                               [0.0, N[0], 0.0, N[1]]
                               ])
            forcegausspt = np.array([[forcex], [forcey]])

            forceelem = forceelem + \
                np.dot(Phi.T, forcegausspt)*wg[g]*0.5*lelem

        F[dofx[0]] = F[dofx[0]]+forceelem[0]
        F[dofy[0]] = F[dofy[0]]+forceelem[1]
        F[dofx[1]] = F[dofx[1]]+forceelem[2]
        F[dofy[1]] = F[dofy[1]]+forceelem[3]

    return F
