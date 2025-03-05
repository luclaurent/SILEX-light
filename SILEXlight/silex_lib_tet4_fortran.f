cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           SILEX CODE
c                    4-node-tetrahedral element
c
c                  Antoine Legay - CNAM - Paris
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c to compile this fortran routines to a python library :
c     f2py3 -c -m silex_lib_tet4_fortran silex_lib_tet4_fortran.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine CrossProduct(a,b,c)

      double precision,intent(in) :: a(3),b(3)
      double precision,intent(out) :: c(3)

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NormVector(a,n)

      integer n,i
      double precision a(n),tmp

      tmp = 0.0d0
      do i=1,n
        tmp=tmp+a(i)*a(i)
      enddo
      NormVector=sqrt(tmp)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function det33_ligne_de_un(a)
      implicit none 
      double precision a(2,3)
      
      det33_ligne_de_un= a(1,2)*a(2,3)
     &                  -a(1,3)*a(2,2)
     &                  -a(1,1)*a(2,3)
     &                  +a(1,1)*a(2,2)
     &                  +a(2,1)*a(1,3)
     &                  -a(2,1)*a(1,2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function det33(a)
      implicit none 
  
      double precision a(3,3)
      
      det33=a(1,1)*a(2,2)*a(3,3)
     &     -a(1,1)*a(2,3)*a(3,2)
     &     -a(2,1)*a(1,2)*a(3,3)
     &     +a(2,1)*a(1,3)*a(3,2)
     &     +a(3,1)*a(1,2)*a(2,3)
     &     -a(3,1)*a(1,3)*a(2,2)
      
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function det44_ligne_de_un(a)
      implicit none 
 
      double precision a(3,4)
      
      det44_ligne_de_un= a(1,2)*a(2,3)*a(3,4)
     &                   -a(1,2)*a(2,4)*a(3,3)
     &                   -a(2,2)*a(1,3)*a(3,4)
     &                   +a(2,2)*a(1,4)*a(3,3)
     &                   +a(3,2)*a(1,3)*a(2,4)
     &                   -a(3,2)*a(1,4)*a(2,3)
     &                   -a(1,1)*a(2,3)*a(3,4)
     &                   +a(1,1)*a(2,4)*a(3,3)
     &                   +a(1,1)*a(2,2)*a(3,4)
     &                   -a(1,1)*a(2,2)*a(3,3)
     &                   -a(1,1)*a(3,2)*a(2,4)
     &                   +a(1,1)*a(3,2)*a(2,3)
     &                   +a(2,1)*a(1,3)*a(3,4)
     &                   -a(2,1)*a(1,4)*a(3,3)
     &                   -a(2,1)*a(1,2)*a(3,4)
     &                   +a(2,1)*a(1,2)*a(3,3)
     &                   +a(2,1)*a(3,2)*a(1,4)
     &                   -a(2,1)*a(3,2)*a(1,3)
     &                   -a(3,1)*a(1,3)*a(2,4)
     &                   +a(3,1)*a(1,4)*a(2,3)
     &                   +a(3,1)*a(1,2)*a(2,4)
     &                   -a(3,1)*a(1,2)*a(2,3)
     &                   -a(3,1)*a(2,2)*a(1,4)
     &                   +a(3,1)*a(1,3)*a(2,2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine StiffnessMatrix(nbnodes,nodes,
     &                           nbelem,elements,
     &                           material,Ik,Jk,Vk)
      implicit none 

      integer nbnodes,nbelem
      double precision nodes(nbnodes,3)
      integer elements(nbelem,4)
      double precision CC(6,6),Vk(12*12*nbelem)
      integer Ik(12*12*nbelem),Jk(12*12*nbelem)
      double precision material(2),young,nu,lambda,mu

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) material
Cf2py intent(out) Ik
Cf2py intent(out) Jk
Cf2py intent(out) Vk

      double precision ke(12,12)
      double precision X(4),Y(4),Z(4)
      integer idnodes(4),dofx(4),dofy(4),dofz(4),dofelem(12)
      integer p,e,i,j
      double precision det_of_sys,Vol
      double precision a34(3,4),a23(2,3)
      double precision beta(4),gamm(4),delt(4)
      double precision det44_ligne_de_un
      double precision det33_ligne_de_un
      double precision BB(6,12)
    
      young = material(1)
      nu    = material(2)
      do i=1,6
        do j=1,6
          CC(i,j)=0.0d0
        enddo
      enddo 
      lambda = nu*young/((1+nu)*(1-2*nu))
      mu     = young/(2*(1+nu))
      do i=1,3
        CC(i,i)=lambda+2*mu
        CC(i+3,i+3)=mu
      enddo
      CC(1,2)=lambda;CC(2,1)=lambda
      CC(1,3)=lambda;CC(3,1)=lambda
      CC(2,3)=lambda;CC(3,2)=lambda

      do i=1,6
        do j=1,12
          BB(i,j)=0.0d0
        enddo
      enddo

      p=1
      do e=1,nbelem

        do i=1,4
          idnodes(i) = elements(e,i)
c                          python indexing
          dofx(i)    = (idnodes(i)-1)*3
          dofy(i)    = (idnodes(i)-1)*3+1
          dofz(i)    = (idnodes(i)-1)*3+2
        enddo
        do i=1,4
          dofelem(1+3*(i-1)) = dofx(i)
          dofelem(2+3*(i-1)) = dofy(i)
          dofelem(3+3*(i-1)) = dofz(i)
        enddo

        do i=1,4
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)
        enddo

        a34(1,1)=X(1);a34(1,2)=X(2);a34(1,3)=X(3);a34(1,4)=X(4)
        a34(2,1)=Y(1);a34(2,2)=Y(2);a34(2,3)=Y(3);a34(2,4)=Y(4)
        a34(3,1)=Z(1);a34(3,2)=Z(2);a34(3,3)=Z(3);a34(3,4)=Z(4)
        det_of_sys = det44_ligne_de_un(a34)
        Vol        = abs(det_of_sys/6)

        a23(1,1) = Y(2);a23(1,2) = Y(3);a23(1,3) = Y(4)
        a23(2,1) = Z(2);a23(2,2) = Z(3);a23(2,3) = Z(4)
        beta(1)=-det33_ligne_de_un(a23)

        a23(1,1) = X(2);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Z(2);a23(2,2) = Z(3);a23(2,3) = Z(4)
        gamm(1)=+det33_ligne_de_un(a23)

        a23(1,1) = X(2);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Y(2);a23(2,2) = Y(3);a23(2,3) = Y(4)
        delt(1)=-det33_ligne_de_un(a23)

        a23(1,1) = Y(1);a23(1,2) = Y(3);a23(1,3) = Y(4)
        a23(2,1) = Z(1);a23(2,2) = Z(3);a23(2,3) = Z(4)
        beta(2)=+det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Z(1);a23(2,2) = Z(3);a23(2,3) = Z(4)
        gamm(2)=-det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Y(1);a23(2,2) = Y(3);a23(2,3) = Y(4)
        delt(2)=+det33_ligne_de_un(a23)

        a23(1,1) = Y(1);a23(1,2) = Y(2);a23(1,3) = Y(4)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(4)
        beta(3)=-det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(4)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(4)
        gamm(3)=+det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(4)
        a23(2,1) = Y(1);a23(2,2) = Y(2);a23(2,3) = Y(4)
        delt(3)=-det33_ligne_de_un(a23)

        a23(1,1) = Y(1);a23(1,2) = Y(2);a23(1,3) = Y(3)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(3)
        beta(4)=+det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(3)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(3)
        gamm(4)=-det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(3)
        a23(2,1) = Y(1);a23(2,2) = Y(2);a23(2,3) = Y(3)
        delt(4)=+det33_ligne_de_un(a23)


        do j=1,4
          BB(1,1+3*(j-1)) = beta(j)/(6.0d0*Vol)
          BB(2,2+3*(j-1)) = gamm(j)/(6.0d0*Vol)
          BB(3,3+3*(j-1)) = delt(j)/(6.0d0*Vol)
          BB(4,2+3*(j-1)) = delt(j)/(6.0d0*Vol)
          BB(4,3+3*(j-1)) = gamm(j)/(6.0d0*Vol)
          BB(5,1+3*(j-1)) = delt(j)/(6.0d0*Vol)
          BB(5,3+3*(j-1)) = beta(j)/(6.0d0*Vol)
          BB(6,1+3*(j-1)) = gamm(j)/(6.0d0*Vol)
          BB(6,2+3*(j-1)) = beta(j)/(6.0d0*Vol)
        enddo

        ke = matmul(transpose(BB),matmul(CC,BB))*Vol

        do i=1,12
          do j=1,12
            Ik(p)=dofelem(i)
            Jk(p)=dofelem(j)
            Vk(p)=ke(i,j)
            p=p+1
          enddo
        enddo

      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE MassMatrix(nbnodes, 
     &                      nodes,
     &                      nelem,
     &                      elements,
     &                      rho,
     &                      Ik,
     &                      Jk,
     &                      Vk)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: nelem,nbnodes
      INTEGER,INTENT(IN) :: elements(nelem,4)
      DOUBLE PRECISION,INTENT(IN) :: nodes(nbnodes,3),rho
      !!
      INTEGER,INTENT(OUT) :: Ik(12*12*nelem),Jk(12*12*nelem)
      DOUBLE PRECISION,INTENT(OUT) :: Vk(12*12*nelem)
      !!
      INTEGER :: npgt,p,e,i,j,g
      DOUBLE PRECISION :: det_of_sys, det44_ligne_de_un
      INTEGER :: dofelem(12),dofx(4),dofy(4),dofz(4),idnodes(4)
      DOUBLE PRECISION :: RG(5),SG(5),TG(5),WG(5),NN(4)
      DOUBLE PRECISION :: A34(3,4),PHI(3,12),me(12,12)
      DOUBLE PRECISION :: Vol
      DOUBLE PRECISION :: X(4),Y(4),Z(4)


    ! Define Gauss points in reference tetrahedral
      npgt=5
      RG=(/1.0/4.0,1.0/2.0,1.0/6.0,1.0/6.0,1.0/6.0/)
      SG=(/1.0/4.0,1.0/6.0,1.0/2.0,1.0/6.0,1.0/6.0/)
      TG=(/1.0/4.0,1.0/6.0,1.0/6.0,1.0/2.0,1.0/6.0/)
      WG=(/-4.0/(5.0*6.0),9.0/(20.0*6.0),9.0/(20.0*6.0)
     & ,9.0/(20.0*6.0),9.0/(20.0*6.0)/)

      p=1

      DO i = 1,3
          DO j =1,12
            PHI(i,j)=0.0d0
          ENDDO
      ENDDO

      DO e =1,nelem
          DO i = 1,4
            idnodes(i) = elements(e,i)
    !                          python indexing
            dofx(i)    = (idnodes(i)-1)*3
            dofy(i)    = (idnodes(i)-1)*3+1
            dofz(i)    = (idnodes(i)-1)*3+2
          ENDDO
          DO i = 1,4
            dofelem(1+3*(i-1)) = dofx(i)
            dofelem(2+3*(i-1)) = dofy(i)
            dofelem(3+3*(i-1)) = dofz(i)
          ENDDO

          DO i = 1,4
            X(i)=nodes(idnodes(i),1)
            Y(i)=nodes(idnodes(i),2)
            Z(i)=nodes(idnodes(i),3)
          ENDDO

          a34(1,1) = X(1);a34(1,2) = X(2)
          a34(1,3) = X(3);a34(1,4) = X(4)
          a34(2,1) = Y(1);a34(2,2) = Y(2)
          a34(2,3) = Y(3);a34(2,4) = Y(4)
          a34(3,1) = Z(1);a34(3,2) = Z(2)
          a34(3,3) = Z(3);a34(3,4) = Z(4)
          det_of_sys = det44_ligne_de_un(a34)
          Vol        = ABS(det_of_sys/6.0)

          DO i = 1,12
            DO j =1,12
                me(i,j)=0.0d0
            ENDDO
          ENDDO

    !       loop over Gauss Points
          DO g=1,npgt

            NN(1)=RG(g);NN(2)=SG(g);NN(3)=TG(g)
            NN(4)=1-RG(g)-SG(g)-TG(g)

            PHI(1,1)=NN(1);PHI(1,4)=NN(2);PHI(1,7)=NN(3);PHI(1,10)=NN(4)
            PHI(2,2)=NN(1);PHI(2,5)=NN(2);PHI(2,8)=NN(3);PHI(2,11)=NN(4)
            PHI(3,3)=NN(1);PHI(3,6)=NN(2);PHI(3,9)=NN(3);PHI(3,12)=NN(4)

            me = me + MATMUL(TRANSPOSE(PHI),PHI)*rho*WG(g)*Vol*6.0

          ENDDO

          DO i = 1,12
            DO j =1,12
                Ik(p)=dofelem(i)
                Jk(p)=dofelem(j)
                Vk(p)=me(i,j)
                p=p+1
            ENDDO
          ENDDO
      ENDDO

      RETURN
      END SUBROUTINE MassMatrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C stress smoothing
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine Compute_stress_strain_error(nbnodes,nodes,
     &                                       nbelem,elements,
     &                                       material,
     &                                       QQ,
     &                                       Sigma,
     &                                       sig_smooth,
     &                                       EpsilonElem,
     &                                       EpsilonNodes,
     &                                       ErrElem,
     &                                       ErrGlob)
      implicit none 

      integer nbnodes,nbelem
      double precision nodes(nbnodes,3)
      integer elements(nbelem,4)
      double precision QQ(3*nbnodes),sigma(nbelem,7)
      double precision sig_smooth(nbnodes,7),ErrElem(nbelem),ErrGlob
      double precision nodalweight(nbnodes),CCinv(6,6),CC(6,6)
      double precision EpsilonElem(nbelem,6),material(2),young,nu
      double precision lambda,mu,EpsilonNodes(nbnodes,6)
      double precision tmpeps(6),tmpsig(6)

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) material
Cf2py intent(in) QQ
Cf2py intent(out) Sigma
Cf2py intent(out) sig_smooth
Cf2py intent(out) EpsilonElem
Cf2py intent(out) EpsilonNodes
Cf2py intent(out) ErrElem
Cf2py intent(out) ErrGlob


      double precision X(4),Y(4),Z(4)
      integer idnodes(4)
      integer p,e,i,j,g,npgt,dofx(4),dofy(4),dofz(4),dofelem(12)
      double precision det_of_sys,Vol,det44_ligne_de_un
      double precision det33_ligne_de_un
      double precision a34(3,4),a23(2,3)
      double precision Cm1S(6),sigdiff(6)
      double precision NormSig,NormSigElt(nbelem)
      double precision RG(5),SG(5),TG(5),WG(5)
      double precision NN(4),sig(7),Q(12),sigvm
      double precision beta(4),gamm(4),delt(4),BB(6,12)

      young = material(1)
      nu    = material(2)
      do i=1,6
        do j=1,6
          CC(i,j)=0.0d0
          CCinv(i,j)=0.0d0
        enddo
      enddo 
      lambda = nu*young/((1+nu)*(1-2*nu))
      mu     = young/(2*(1+nu))
      do i=1,3
        CC(i,i)=lambda+2*mu
        CC(i+3,i+3)=mu
        CCinv(i,i)=1.0/young
        CCinv(i+3,i+3)=2.0*(1+nu)/young
      enddo
      CC(1,2)=lambda;CC(2,1)=lambda
      CC(1,3)=lambda;CC(3,1)=lambda
      CC(2,3)=lambda;CC(3,2)=lambda
      CCinv(1,2)=-nu/young;CCinv(2,1)=-nu/young
      CCinv(1,3)=-nu/young;CCinv(3,1)=-nu/young
      CCinv(2,3)=-nu/young;CCinv(3,2)=-nu/young

      do e=1,nbelem
        do i=1,4
          idnodes(i) = elements(e,i)
c                          fortran indexing
          dofx(i)    = (idnodes(i)-1)*3+1
          dofy(i)    = (idnodes(i)-1)*3+2
          dofz(i)    = (idnodes(i)-1)*3+3
        enddo
        do i=1,4
          dofelem(1+3*(i-1)) = dofx(i)
          dofelem(2+3*(i-1)) = dofy(i)
          dofelem(3+3*(i-1)) = dofz(i)
        enddo

        do i=1,4
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)
        enddo

        do i=1,12
          Q(i)=QQ(dofelem(i))
        enddo

        a34(1,1)=X(1);a34(1,2)=X(2);a34(1,3)=X(3);a34(1,4)=X(4)
        a34(2,1)=Y(1);a34(2,2)=Y(2);a34(2,3)=Y(3);a34(2,4)=Y(4)
        a34(3,1)=Z(1);a34(3,2)=Z(2);a34(3,3)=Z(3);a34(3,4)=Z(4)
        det_of_sys = det44_ligne_de_un(a34)
        Vol        = abs(det_of_sys/6)

        a23(1,1) = Y(2);a23(1,2) = Y(3);a23(1,3) = Y(4)
        a23(2,1) = Z(2);a23(2,2) = Z(3);a23(2,3) = Z(4)
        beta(1)=-det33_ligne_de_un(a23)

        a23(1,1) = X(2);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Z(2);a23(2,2) = Z(3);a23(2,3) = Z(4)
        gamm(1)=+det33_ligne_de_un(a23)

        a23(1,1) = X(2);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Y(2);a23(2,2) = Y(3);a23(2,3) = Y(4)
        delt(1)=-det33_ligne_de_un(a23)
  

        a23(1,1) = Y(1);a23(1,2) = Y(3);a23(1,3) = Y(4)
        a23(2,1) = Z(1);a23(2,2) = Z(3);a23(2,3) = Z(4)
        beta(2)=+det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Z(1);a23(2,2) = Z(3);a23(2,3) = Z(4)
        gamm(2)=-det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(3);a23(1,3) = X(4)
        a23(2,1) = Y(1);a23(2,2) = Y(3);a23(2,3) = Y(4)
        delt(2)=+det33_ligne_de_un(a23)

        a23(1,1) = Y(1);a23(1,2) = Y(2);a23(1,3) = Y(4)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(4)
        beta(3)=-det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(4)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(4)
        gamm(3)=+det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(4)
        a23(2,1) = Y(1);a23(2,2) = Y(2);a23(2,3) = Y(4)
        delt(3)=-det33_ligne_de_un(a23)

        a23(1,1) = Y(1);a23(1,2) = Y(2);a23(1,3) = Y(3)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(3)
        beta(4)=+det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(3)
        a23(2,1) = Z(1);a23(2,2) = Z(2);a23(2,3) = Z(3)
        gamm(4)=-det33_ligne_de_un(a23)

        a23(1,1) = X(1);a23(1,2) = X(2);a23(1,3) = X(3)
        a23(2,1) = Y(1);a23(2,2) = Y(2);a23(2,3) = Y(3)
        delt(4)=+det33_ligne_de_un(a23)

        do i=1,6
          do j=1,12
            BB(i,j)=0.0d0
          enddo
        enddo

        do j=1,4
          BB(1,1+3*(j-1)) = beta(j)/(6.0d0*Vol)
          BB(2,2+3*(j-1)) = gamm(j)/(6.0d0*Vol)
          BB(3,3+3*(j-1)) = delt(j)/(6.0d0*Vol)
          BB(4,2+3*(j-1)) = delt(j)/(6.0d0*Vol)
          BB(4,3+3*(j-1)) = gamm(j)/(6.0d0*Vol)
          BB(5,1+3*(j-1)) = delt(j)/(6.0d0*Vol)
          BB(5,3+3*(j-1)) = beta(j)/(6.0d0*Vol)
          BB(6,1+3*(j-1)) = gamm(j)/(6.0d0*Vol)
          BB(6,2+3*(j-1)) = beta(j)/(6.0d0*Vol)
        enddo

        Sig(1:6)=matmul(CC,matmul(BB,Q))

        sigvm  = 1.5*(Sig(1)**2+Sig(2)**2+Sig(3)**2
     &         +2.0*Sig(4)**2+2.0*Sig(5)**2+2.0*Sig(6)**2)
     &         -0.5*(Sig(1)+Sig(2)+Sig(3))**2
        sigvm  = sqrt(sigvm)
        Sig(7) = sigvm

        do i=1,7
          Sigma(e,i)=Sig(i)
        enddo
        do i=1,6
          tmpsig(i)=Sig(i)
        enddo
        tmpeps=matmul(CCinv,tmpsig)
        do i=1,6
          EpsilonElem(e,i)=tmpeps(i)
        enddo
      enddo

      do i=1,nbnodes
        nodalweight(i)=0.0d0
        do j=1,7
          sig_smooth(i,j)=0.0d0
        enddo
      enddo

      p=1
      do e=1,nbelem

        do i=1,4
          idnodes(i) = elements(e,i)
        enddo
        do i=1,4
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)
        enddo

        a34(1,1)=X(1); a34(1,2)=X(2); a34(1,3)=X(3); a34(1,4)=X(4)
        a34(2,1)=Y(1); a34(2,2)=Y(2); a34(2,3)=Y(3); a34(2,4)=Y(4)
        a34(3,1)=Z(1); a34(3,2)=Z(2); a34(3,3)=Z(3); a34(3,4)=Z(4)
        det_of_sys = det44_ligne_de_un(a34)
        Vol        = abs(det_of_sys/6.0d0)
 
        do i=1,4
          nodalweight(idnodes(i))=nodalweight(idnodes(i))+Vol
          do j=1,7
            sig_smooth(idnodes(i),j)=sig_smooth(idnodes(i),j)
     &           +sigma(e,j)*Vol
          enddo
        enddo

      enddo

c norm the smooth stress
      do i=1,nbnodes
        do j=1,7
          sig_smooth(i,j)=sig_smooth(i,j)/nodalweight(i)
        enddo
        do j=1,6
          tmpsig(j)=sig_smooth(i,j)
        enddo
        tmpeps=matmul(CCinv,tmpsig)
        do j=1,6
          EpsilonNodes(i,j)=tmpeps(j)
        enddo
      enddo

c Compute Error


c     Define Gauss points in reference tetrahedral

c      npgt=4
c      RG(1)=0.58541020;SG(1)=0.13819660;TG(1)=0.13819660;WG(1)=0.25/6.0
c      RG(2)=0.13819660;SG(2)=0.58541020;TG(2)=0.13819660;WG(2)=0.25/6.0
c      RG(3)=0.13819660;SG(3)=0.13819660;TG(3)=0.58541020;WG(3)=0.25/6.0
c      RG(4)=0.13819660;SG(4)=0.13819660;TG(4)=0.13819660;WG(4)=0.25/6.0

      npgt=5
      RG(1)=1.0/4.0; SG(1)=1.0/4.0; TG(1)=1.0/4.0; WG(1)=-4.0/(5.0*6.0)
      RG(2)=1.0/2.0; SG(2)=1.0/6.0; TG(2)=1.0/6.0; WG(2)=9.0/(20.0*6.0)
      RG(3)=1.0/6.0; SG(3)=1.0/2.0; TG(3)=1.0/6.0; WG(3)=9.0/(20.0*6.0)
      RG(4)=1.0/6.0; SG(4)=1.0/6.0; TG(4)=1.0/2.0; WG(4)=9.0/(20.0*6.0)
      RG(5)=1.0/6.0; SG(5)=1.0/6.0; TG(5)=1.0/6.0; WG(5)=9.0/(20.0*6.0)

c     compute energy norm of uh and sigma-smooth
      NormSig=0.0d0
      do e=1,nbelem

        do i=1,4
          idnodes(i) = elements(e,i)
        enddo

        do i=1,4
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)
        enddo

        a34(1,1)=X(1); a34(1,2)=X(2); a34(1,3)=X(3); a34(1,4)=X(4)
        a34(2,1)=Y(1); a34(2,2)=Y(2); a34(2,3)=Y(3); a34(2,4)=Y(4)
        a34(3,1)=Z(1); a34(3,2)=Z(2); a34(3,3)=Z(3); a34(3,4)=Z(4)
        det_of_sys = det44_ligne_de_un(a34)
        Vol        = abs(det_of_sys/6.0d0)

        NormSigElt(e) = 0.0d0
        ErrElem(e)    = 0.0d0

c       loop over Gauss Points
        do g=1,npgt
          NN(1)=RG(g);NN(2)=SG(g);NN(3)=TG(g);NN(4)=1-RG(g)-SG(g)-TG(g)
c         evaluate sigma at gauss point
          do i=1,6
            sig(i)=0.0d0
          enddo
          do i=1,4
            do j=1,6
              sig(j)=sig(j)+sig_smooth(idnodes(i),j)*NN(i)
            enddo
          enddo

          do i=1,6
            sigdiff(i)=sig(i)-sigma(e,i)
          enddo

          do i=1,6
            Cm1S(i)=0.0d0
            do j=1,6
              Cm1S(i)=Cm1S(i)+CCinv(i,j)*sig(j)
            enddo
          enddo

          do i=1,6
            NormSigElt(e)=NormSigElt(e)+sig(i)*Cm1S(i)*WG(g)*Vol*6.0
          enddo

          do i=1,6
            Cm1S(i)=0.0d0
            do j=1,6
              Cm1S(i)=Cm1S(i)+CCinv(i,j)*sigdiff(j)
            enddo
          enddo

          do i=1,6
            ErrElem(e)=ErrElem(e)+sigdiff(i)*Cm1S(i)*WG(g)*Vol*6.0
          enddo

        enddo

        NormSig = NormSig+NormSigElt(e)

      enddo

      ErrGlob=0.0d0
      do e=1,nbelem
        ErrElem(e) = ErrElem(e)/NormSig
        ErrGlob    = ErrGlob + ErrElem(e)
      enddo

      ErrGlob    = sqrt(ErrGlob)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine forceonsurface(nbnodes,nodes,
     &                          nbelem,elements,
     &                          Press,
     &                          direction,
     &                          Fp)
      implicit none 

      integer nbnodes,nbelem
      double precision  nodes(nbnodes,3),Fp(3*nbnodes)
      integer elements(nbelem,3)
      double precision  Press

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) Press
Cf2py intent(in) direction
Cf2py intent(out) Fp

      double precision X(3),Y(3),Z(3),v12(3),v13(3),VecN(3),NormVecN
      double precision Ae,FpPerNode(3),NormVector,surf,direction(3)
      integer i,e
      integer idnodes(3),dofx(3),dofy(3),dofz(3)
      double precision tmp
      integer pressure_flag

      tmp=0.0d0
      do i=1,3
        tmp=tmp+direction(i)*direction(i)
      enddo
      tmp=dsqrt(tmp)

      if (tmp.le.1.0e-5) then 
        pressure_flag=1
      else
        pressure_flag=0
        do i=1,3
          direction(i)=direction(i)/tmp
        enddo
      endif

      do i=1,3*nbnodes
        Fp(i)=0.0d0
      enddo
      surf=0.0d0
      do e=1,nbelem
        do i=1,3
          idnodes(i) = elements(e,i)
c                          fortran indexing
          dofx(i)    = (idnodes(i)-1)*3+1
          dofy(i)    = (idnodes(i)-1)*3+2
          dofz(i)    = (idnodes(i)-1)*3+3
        enddo
        do i=1,3
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
          Z(i)=nodes(idnodes(i),3)
        enddo

        v12(1)=X(2)-X(1)
        v12(2)=Y(2)-Y(1)
        v12(3)=Z(2)-Z(1)

        v13(1)=X(3)-X(1)
        v13(2)=Y(3)-Y(1)
        v13(3)=Z(3)-Z(1)

        call CrossProduct(v12,v13,VecN)
        Ae = 0.5*NormVector(VecN,3) ! area of the triangle
        NormVecN=NormVector(VecN,3)

        do i=1,3
          VecN(i)=VecN(i)/NormVecN
        enddo

        if (pressure_flag.eq.1) then
          do i=1,3
            FpPerNode(i) = VecN(i)*press*Ae/3.0d0
          enddo
        else
          do i=1,3
            FpPerNode(i) = direction(i)*press*Ae/3.0d0
          enddo
        endif

        do i=1,3
          Fp(dofx(i))=Fp(dofx(i))+FpPerNode(1)
          Fp(dofy(i))=Fp(dofy(i))+FpPerNode(2)
          Fp(dofz(i))=Fp(dofz(i))+FpPerNode(3)
        enddo
        surf=surf+Ae
      enddo
      !write(*,*) ' --- >  surf ',surf

      return
 
      end

