#############################################################################
#                           SILEX CODE
#                    4-node-tetrahedral element
#
#                  Antoine Legay - CNAM - Paris
#############################################################################
#      Define a few functions
#############################################################################

import numpy
import scipy.sparse


def det33_ligne_de_un(a):
    return a[0, 1]*a[1, 2]-a[0, 2]*a[1, 1]+a[0, 0]*(a[1, 1]-a[1, 2])+a[1, 0]*(a[0, 2]-a[0, 1])


def det44_ligne_de_un(a):
    return a[0, 1]*(a[1, 2]*a[2, 3]-a[1, 3]*a[2, 2])+a[1, 1]*(-a[0, 2]*a[2, 3]+a[0, 3]*a[2, 2])+a[2, 1]*(a[0, 2]*a[1, 3]-a[0, 3]*a[1, 2])+a[0, 0]*(-a[1, 2]*a[2, 3]+a[1, 3]*a[2, 2]+a[1, 1]*(a[2, 3]-a[2, 2])+a[2, 1]*(-a[1, 3]+a[1, 2]))+a[1, 0]*(a[0, 2]*a[2, 3]-a[0, 3]*a[2, 2]+a[0, 1]*(-a[2, 3]+a[2, 2])+a[2, 1]*(a[0, 3]-a[0, 2]))+a[2, 0]*(-a[0, 2]*a[1, 3]+a[0, 1]*(a[1, 3]-a[1, 2])+a[0, 3]*(-a[1, 1]+a[1, 2])+a[0, 2]*a[1, 1])


#############################################################################
#      compute stiffness matrix
#############################################################################

def stiffnessmatrix(nodes, elements, material):

    nelem = elements.shape[0]
    Ik = numpy.zeros(12*12*nelem, dtype=int)
    Jk = numpy.zeros(12*12*nelem, dtype=int)
    Vk = numpy.zeros(12*12*nelem, dtype=float)
    beta = numpy.zeros((4), dtype=float)
    gamm = numpy.zeros((4), dtype=float)
    delt = numpy.zeros((4), dtype=float)
    A34 = numpy.zeros((3, 4), dtype=float)
    A23 = numpy.zeros((3, 4), dtype=float)
    B = numpy.zeros((6, 12), dtype=float)
    dofelem = numpy.zeros(12, dtype=int)
    dofx = numpy.zeros(4, dtype=int)
    dofy = numpy.zeros(4, dtype=int)
    dofz = numpy.zeros(4, dtype=int)
    idnodes = numpy.zeros(4, dtype=int)
    young = material[0]
    nu = material[1]
    lamb = nu*young/((1+nu)*(1-2*nu))
    mu = young/(2*(1+nu))
    C = numpy.array([[lamb+2*mu, lamb, lamb, 0.0, 0.0, 0.0],
                     [lamb, lamb+2*mu, lamb, 0.0, 0.0, 0.0],
                     [lamb, lamb, lamb+2*mu, 0.0, 0.0, 0.0],
                     [0.0, 0.0, 0.0, mu, 0.0, 0.0],
                     [0.0, 0.0, 0.0, 0.0, mu, 0.0],
                     [0.0, 0.0, 0.0, 0.0, 0.0, mu]])

    p = 0

    for e in range(nelem):
        idnodes[:] = elements[e, :]-1
        dofx[:] = (idnodes)*3
        dofy[:] = (idnodes)*3+1
        dofz[:] = (idnodes)*3+2
        dofelem[0] = dofx[0]
        dofelem[1] = dofy[0]
        dofelem[2] = dofz[0]
        dofelem[3] = dofx[1]
        dofelem[4] = dofy[1]
        dofelem[5] = dofz[1]
        dofelem[6] = dofx[2]
        dofelem[7] = dofy[2]
        dofelem[8] = dofz[2]
        dofelem[9] = dofx[3]
        dofelem[10] = dofy[3]
        dofelem[11] = dofz[3]

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]
        Z = nodes[idnodes, 2]
        A34[0, :] = X
        A34[1, :] = Y
        A34[2, :] = Z
        det_of_sys = det44_ligne_de_un(A34)
        Vol = abs(det_of_sys/6)

        A23[0, 0] = Y[1]
        A23[0, 1] = Y[2]
        A23[0, 2] = Y[3]
        A23[1, 0] = Z[1]
        A23[1, 1] = Z[2]
        A23[1, 2] = Z[3]
        beta[0] = -det33_ligne_de_un(A23)
        A23[0, 0] = X[1]
        A23[0, 1] = X[2]
        A23[0, 2] = X[3]
        gamm[0] = det33_ligne_de_un(A23)
        A23[1, 0] = Y[1]
        A23[1, 1] = Y[2]
        A23[1, 2] = Y[3]
        delt[0] = -det33_ligne_de_un(A23)

        A23[0, 0] = Y[0]
        A23[0, 1] = Y[2]
        A23[0, 2] = Y[3]
        A23[1, 0] = Z[0]
        A23[1, 1] = Z[2]
        A23[1, 2] = Z[3]
        beta[1] = det33_ligne_de_un(A23)
        A23[0, 0] = X[0]
        A23[0, 1] = X[2]
        A23[0, 2] = X[3]
        gamm[1] = -det33_ligne_de_un(A23)
        A23[1, 0] = Y[0]
        A23[1, 1] = Y[2]
        A23[1, 2] = Y[3]
        delt[1] = det33_ligne_de_un(A23)

        A23[0, 0] = Y[0]
        A23[0, 1] = Y[1]
        A23[0, 2] = Y[3]
        A23[1, 0] = Z[0]
        A23[1, 1] = Z[1]
        A23[1, 2] = Z[3]
        beta[2] = -det33_ligne_de_un(A23)
        A23[0, 0] = X[0]
        A23[0, 1] = X[1]
        A23[0, 2] = X[3]
        gamm[2] = det33_ligne_de_un(A23)
        A23[1, 0] = Y[0]
        A23[1, 1] = Y[1]
        A23[1, 2] = Y[3]
        delt[2] = -det33_ligne_de_un(A23)

        A23[0, 0] = Y[0]
        A23[0, 1] = Y[1]
        A23[0, 2] = Y[2]
        A23[1, 0] = Z[0]
        A23[1, 1] = Z[1]
        A23[1, 2] = Z[2]
        beta[3] = det33_ligne_de_un(A23)
        A23[0, 0] = X[0]
        A23[0, 1] = X[1]
        A23[0, 2] = X[2]
        gamm[3] = -det33_ligne_de_un(A23)
        A23[1, 0] = Y[0]
        A23[1, 1] = Y[1]
        A23[1, 2] = Y[2]
        delt[3] = det33_ligne_de_un(A23)

        B[0, 0] = beta[0]
        B[0, 3] = beta[1]
        B[0, 6] = beta[2]
        B[0, 9] = beta[3]
        B[1, 1] = gamm[0]
        B[1, 4] = gamm[1]
        B[1, 7] = gamm[2]
        B[1, 10] = gamm[3]
        B[2, 2] = delt[0]
        B[2, 5] = delt[1]
        B[2, 8] = delt[2]
        B[2, 11] = delt[3]
        B[3, 2] = gamm[0]
        B[3, 5] = gamm[1]
        B[3, 8] = gamm[2]
        B[3, 11] = gamm[3]
        B[3, 1] = delt[0]
        B[3, 4] = delt[1]
        B[3, 7] = delt[2]
        B[3, 10] = delt[3]
        B[4, 0] = delt[0]
        B[4, 3] = delt[1]
        B[4, 6] = delt[2]
        B[4, 9] = delt[3]
        B[4, 2] = beta[0]
        B[4, 5] = beta[1]
        B[4, 8] = beta[2]
        B[4, 11] = beta[3]
        B[5, 0] = gamm[0]
        B[5, 3] = gamm[1]
        B[5, 6] = gamm[2]
        B[5, 9] = gamm[3]
        B[5, 1] = beta[0]
        B[5, 4] = beta[1]
        B[5, 7] = beta[2]
        B[5, 10] = beta[3]

        ke = numpy.dot(B.T, numpy.dot(C, B))*Vol/(det_of_sys*det_of_sys)

        for i in range(12):
            Ik[p:(p+12)] = dofelem[i]
            Jk[p:(p+12)] = dofelem[:]
            Vk[p:(p+12)] = ke[i, :]
            p = p+12

    return Ik, Jk, Vk

#############################################################


def forceonsurface(nodes, elements, press, direction):
    nbnodes = nodes.shape[0]
    nbelem = elements.shape[0]
    idnodes = numpy.zeros(3, dtype=int)
    dofx = numpy.zeros(3, dtype=int)
    dofy = numpy.zeros(3, dtype=int)
    dofz = numpy.zeros(3, dtype=int)
    Fp = numpy.zeros(3*nbnodes, dtype=float)
    v12 = numpy.zeros(3, dtype=float)
    v13 = numpy.zeros(3, dtype=float)

    tmp = numpy.linalg.norm(numpy.array(direction))
    if (tmp < 1.0e-5):
        pressure_flag = 1
    else:
        pressure_flag = 0
        direction = numpy.array(direction)/tmp

    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1
        dofx[:] = (idnodes)*3
        dofy[:] = (idnodes)*3+1
        dofz[:] = (idnodes)*3+2

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]
        Z = nodes[idnodes, 2]

        v12[0] = X[1]-X[0]
        v12[1] = Y[1]-Y[0]
        v12[2] = Z[1]-Z[0]

        v13[0] = X[2]-X[0]
        v13[1] = Y[2]-Y[0]
        v13[2] = Z[2]-Z[0]

        VecN = numpy.cross(v12, v13)
        NormVecN = numpy.linalg.norm(VecN)
        Ae = 0.5*NormVecN  # area of the triangle
        VecN = VecN/NormVecN

        if pressure_flag == 1:
            FpPerNode = VecN*press*Ae/3.0
        else:
            FpPerNode = direction*press*Ae/3.0

        for i in range(3):
            Fp[dofx[i]] = Fp[dofx[i]]+FpPerNode[0]
            Fp[dofy[i]] = Fp[dofy[i]]+FpPerNode[1]
            Fp[dofz[i]] = Fp[dofz[i]]+FpPerNode[2]

    return Fp

############################################################


def compute_stress_strain_error(nodes, elements, material, QQ):
    nbnodes = nodes.shape[0]
    nbelem = elements.shape[0]
    sigma = numpy.zeros((nbelem, 7), dtype=float)
    sig_smooth = numpy.zeros((nbnodes, 7), dtype=float)
    ErrElem = numpy.zeros(nbelem, dtype=float)
    nodalweight = numpy.zeros(nbnodes, dtype=float)
    EpsilonElem = numpy.zeros((nbelem, 6), dtype=float)
    EpsilonNodes = numpy.zeros((nbnodes, 6), dtype=float)
    tmpeps = numpy.zeros(6, dtype=float)
    tmpsig = numpy.zeros(6, dtype=float)
    X = numpy.zeros(4, dtype=float)
    Y = numpy.zeros(4, dtype=float)
    Z = numpy.zeros(4, dtype=float)
    idnodes = numpy.zeros(4, dtype=int)
    dofx = numpy.zeros(4, dtype=int)
    dofy = numpy.zeros(4, dtype=int)
    dofz = numpy.zeros(4, dtype=int)
    dofelem = numpy.zeros(12, dtype=int)
    A34 = numpy.zeros((3, 4), dtype=float)
    A23 = numpy.zeros((2, 3), dtype=float)
    Cm1S = numpy.zeros(6, dtype=float)
    sigdiff = numpy.zeros(6, dtype=float)
    NormSigElt = numpy.zeros(nbelem, dtype=float)
    RG = numpy.zeros(5, dtype=float)
    SG = numpy.zeros(5, dtype=float)
    TG = numpy.zeros(5, dtype=float)
    WG = numpy.zeros(5, dtype=float)
    NN = numpy.zeros(4, dtype=float)
    sig = numpy.zeros(7, dtype=float)
    Q = numpy.zeros(12, dtype=float)
    beta = numpy.zeros((4), dtype=float)
    gamm = numpy.zeros((4), dtype=float)
    d = numpy.zeros((4), dtype=float)
    B = numpy.zeros((6, 12), dtype=float)

    young = material[0]
    nu = material[1]
    lamb = nu*young/((1+nu)*(1-2*nu))
    mu = young/(2*(1+nu))
    CC = numpy.array([[lamb+2*mu, lamb, lamb, 0.0, 0.0, 0.0],
                      [lamb, lamb+2*mu, lamb, 0.0, 0.0, 0.0],
                      [lamb, lamb, lamb+2*mu, 0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, mu, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, mu, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 0.0, mu]])

    CCinv = numpy.array([[1.0/young, -nu/young, -nu/young, 0.0, 0.0, 0.0],
                         [-nu/young, 1.0/young, -nu/young, 0.0, 0.0, 0.0],
                         [-nu/young, -nu/young, 1.0/young, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 2.0*(1+nu)/young, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 2.0*(1+nu)/young, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 2.0*(1+nu)/young]])

    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1
        dofx[:] = (idnodes)*3
        dofy[:] = (idnodes)*3+1
        dofz[:] = (idnodes)*3+2
        dofelem[0] = dofx[0]
        dofelem[1] = dofy[0]
        dofelem[2] = dofz[0]
        dofelem[3] = dofx[1]
        dofelem[4] = dofy[1]
        dofelem[5] = dofz[1]
        dofelem[6] = dofx[2]
        dofelem[7] = dofy[2]
        dofelem[8] = dofz[2]
        dofelem[9] = dofx[3]
        dofelem[10] = dofy[3]
        dofelem[11] = dofz[3]

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]
        Z = nodes[idnodes, 2]

        Q[:] = QQ[dofelem]

        A34[0, :] = X
        A34[1, :] = Y
        A34[2, :] = Z
        det_of_sys = det44_ligne_de_un(A34)
        Vol = abs(det_of_sys/6)
        tmp = Vol/(det_of_sys*det_of_sys)

        A23[0, 0] = Y[1]
        A23[0, 1] = Y[2]
        A23[0, 2] = Y[3]
        A23[1, 0] = Z[1]
        A23[1, 1] = Z[2]
        A23[1, 2] = Z[3]
        beta[0] = -det33_ligne_de_un(A23)
        A23[0, 0] = X[1]
        A23[0, 1] = X[2]
        A23[0, 2] = X[3]
        gamm[0] = det33_ligne_de_un(A23)
        A23[1, 0] = Y[1]
        A23[1, 1] = Y[2]
        A23[1, 2] = Y[3]
        d[0] = -det33_ligne_de_un(A23)

        A23[0, 0] = Y[0]
        A23[0, 1] = Y[2]
        A23[0, 2] = Y[3]
        A23[1, 0] = Z[0]
        A23[1, 1] = Z[2]
        A23[1, 2] = Z[3]
        beta[1] = det33_ligne_de_un(A23)
        A23[0, 0] = X[0]
        A23[0, 1] = X[2]
        A23[0, 2] = X[3]
        gamm[1] = -det33_ligne_de_un(A23)
        A23[1, 0] = Y[0]
        A23[1, 1] = Y[2]
        A23[1, 2] = Y[3]
        d[1] = det33_ligne_de_un(A23)

        A23[0, 0] = Y[0]
        A23[0, 1] = Y[1]
        A23[0, 2] = Y[3]
        A23[1, 0] = Z[0]
        A23[1, 1] = Z[1]
        A23[1, 2] = Z[3]
        beta[2] = -det33_ligne_de_un(A23)
        A23[0, 0] = X[0]
        A23[0, 1] = X[1]
        A23[0, 2] = X[3]
        gamm[2] = det33_ligne_de_un(A23)
        A23[1, 0] = Y[0]
        A23[1, 1] = Y[1]
        A23[1, 2] = Y[3]
        d[2] = -det33_ligne_de_un(A23)

        A23[0, 0] = Y[0]
        A23[0, 1] = Y[1]
        A23[0, 2] = Y[2]
        A23[1, 0] = Z[0]
        A23[1, 1] = Z[1]
        A23[1, 2] = Z[2]
        beta[3] = det33_ligne_de_un(A23)
        A23[0, 0] = X[0]
        A23[0, 1] = X[1]
        A23[0, 2] = X[2]
        gamm[3] = -det33_ligne_de_un(A23)
        A23[1, 0] = Y[0]
        A23[1, 1] = Y[1]
        A23[1, 2] = Y[2]
        d[3] = det33_ligne_de_un(A23)

        B[0, 0] = beta[0]
        B[0, 3] = beta[1]
        B[0, 6] = beta[2]
        B[0, 9] = beta[3]
        B[1, 1] = gamm[0]
        B[1, 4] = gamm[1]
        B[1, 7] = gamm[2]
        B[1, 10] = gamm[3]
        B[2, 2] = d[0]
        B[2, 5] = d[1]
        B[2, 8] = d[2]
        B[2, 11] = d[3]
        B[3, 2] = gamm[0]
        B[3, 5] = gamm[1]
        B[3, 8] = gamm[2]
        B[3, 11] = gamm[3]
        B[3, 1] = d[0]
        B[3, 4] = d[1]
        B[3, 7] = d[2]
        B[3, 10] = d[3]
        B[4, 0] = d[0]
        B[4, 3] = d[1]
        B[4, 6] = d[2]
        B[4, 9] = d[3]
        B[4, 2] = beta[0]
        B[4, 5] = beta[1]
        B[4, 8] = beta[2]
        B[4, 11] = beta[3]
        B[5, 0] = gamm[0]
        B[5, 3] = gamm[1]
        B[5, 6] = gamm[2]
        B[5, 9] = gamm[3]
        B[5, 1] = beta[0]
        B[5, 4] = beta[1]
        B[5, 7] = beta[2]
        B[5, 10] = beta[3]

        Sig = numpy.dot(numpy.dot(CC, B), Q)/det_of_sys
        sigma[e, 0:6] = Sig[:]
        sigma[e, 6] = numpy.sqrt(1.5*(Sig[0]**2+Sig[1]**2+Sig[2]**2+2.0*Sig[3]
                                 ** 2+2.0*Sig[4]**2+2.0*Sig[5]**2)-0.5*(Sig[0]+Sig[1]+Sig[2])**2)
        EpsilonElem[e, :] = numpy.dot(CCinv, Sig)

    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]
        Z = nodes[idnodes, 2]

        A34[0, :] = X
        A34[1, :] = Y
        A34[2, :] = Z
        det_of_sys = det44_ligne_de_un(A34)
        Vol = abs(det_of_sys/6)

        nodalweight[idnodes] = nodalweight[idnodes]+Vol
        sig_smooth[idnodes, :] = sig_smooth[idnodes, :]+sigma[e, :]*Vol

# norm the smooth stress
    for i in range(nbnodes):
        sig_smooth[i, :] = sig_smooth[i, :]/nodalweight[i]
        tmpsig[:] = sig_smooth[i, 0:6]
        tmpeps = numpy.dot(CCinv, tmpsig)
        EpsilonNodes[i, :] = tmpeps[:]

# Compute Error

# Define Gauss points in reference tetrahedral
    npgt = 5
    RG[0] = 1.0/4.0
    SG[0] = 1.0/4.0
    TG[0] = 1.0/4.0
    WG[0] = -4.0/(5.0*6.0)
    RG[1] = 1.0/2.0
    SG[1] = 1.0/6.0
    TG[1] = 1.0/6.0
    WG[1] = 9.0/(20.0*6.0)
    RG[2] = 1.0/6.0
    SG[2] = 1.0/2.0
    TG[2] = 1.0/6.0
    WG[2] = 9.0/(20.0*6.0)
    RG[3] = 1.0/6.0
    SG[3] = 1.0/6.0
    TG[3] = 1.0/2.0
    WG[3] = 9.0/(20.0*6.0)
    RG[4] = 1.0/6.0
    SG[4] = 1.0/6.0
    TG[4] = 1.0/6.0
    WG[4] = 9.0/(20.0*6.0)

# compute energy norm of uh and sigma-smooth
    NormSig = 0.0
    for e in range(nbelem):
        idnodes[:] = elements[e, :]-1

        X = nodes[idnodes, 0]
        Y = nodes[idnodes, 1]
        Z = nodes[idnodes, 2]

        A34[0, :] = X
        A34[1, :] = Y
        A34[2, :] = Z
        det_of_sys = det44_ligne_de_un(A34)
        Vol = abs(det_of_sys/6)

# loop over Gauss Points
        for g in range(npgt):
            NN[0] = RG[g]
            NN[1] = SG[g]
            NN[2] = TG[g]
            NN[3] = 1-RG[g]-SG[g]-TG[g]
            sig = sig_smooth[idnodes[0], 0:6]*NN[0]+sig_smooth[idnodes[1], 0:6] * \
                NN[1]+sig_smooth[idnodes[2], 0:6]*NN[2] + \
                sig_smooth[idnodes[3], 0:6]*NN[3]

            sigdiff = sig-sigma[e, 0:6]

            NormSigElt[e] = NormSigElt[e] + \
                numpy.dot(sig, numpy.dot(CCinv, sig))*WG[g]*Vol*6.0

            ErrElem[e] = ErrElem[e] + \
                numpy.dot(sigdiff, numpy.dot(CCinv, sigdiff))*WG[g]*Vol*6.0

        NormSig = NormSig+NormSigElt[e]

    ErrElem = ErrElem/NormSig
    ErrGlob = scipy.sqrt(scipy.sum(ErrElem))

    return sigma, sig_smooth, EpsilonElem, EpsilonNodes, ErrElem, ErrGlob
