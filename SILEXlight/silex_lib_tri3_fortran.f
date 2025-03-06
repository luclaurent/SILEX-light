cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           SILEX CODE
c                   3-node-triangle element
c
c                  Antoine Legay - CNAM - Paris
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c to compile this fortran routines to a python library :
c     f2py3 -c -m silex_lib_tri3_fortran silex_lib_tri3_fortran.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
      subroutine ElementalStiffness(X,Y,thickness,CC,ke)
      implicit none 

      double precision ke(6,6)
      double precision X(3),Y(3)
      double precision det_of_sys,Area
      double precision a23(2,3)
      double precision b(3),c(3)
      double precision det33_ligne_de_un
      double precision BB(3,6),CC(3,3),CB(3,6)
      double precision thickness

      integer i,j,k

      a23(1,1) = X(1)
      a23(1,2) = X(2)
      a23(1,3) = X(3)
      a23(2,1) = Y(1)
      a23(2,2) = Y(2)
      a23(2,3) = Y(3)
      det_of_sys=det33_ligne_de_un(a23)
      Area=abs(0.5*det_of_sys)

      b(1)=-(Y(3)-Y(2))/det_of_sys
      c(1)=+(X(3)-X(2))/det_of_sys

      b(2)=+(Y(3)-Y(1))/det_of_sys
      c(2)=-(X(3)-X(1))/det_of_sys

      b(3)=-(Y(2)-Y(1))/det_of_sys
      c(3)=+(X(2)-X(1))/det_of_sys

      BB(1,1)=b(1)
      BB(1,2)=b(2)
      BB(1,3)=b(3)
      BB(1,4)=0.0d0
      BB(1,5)=0.0d0
      BB(1,6)=0.0d0
      BB(2,1)=0.0d0
      BB(2,2)=0.0d0
      BB(2,3)=0.0d0
      BB(2,4)=c(1)
      BB(2,5)=c(2)
      BB(2,6)=c(3)
      BB(3,1)=c(1)
      BB(3,2)=c(2)
      BB(3,3)=c(3)
      BB(3,4)=b(1)
      BB(3,5)=b(2)
      BB(3,6)=b(3)

      do i=1,3
        do j=1,6
          CB(i,j)=0.0d0
          do k=1,3
            CB(i,j)=CB(i,j)+CC(i,k)*BB(k,j)
          enddo
        enddo
      enddo

      do i=1,6
        do j=1,6
          ke(i,j)=0.0d0
          do k=1,3
            ke(i,j)=ke(i,j)+BB(k,i)*CB(k,j)*Area*thickness
          enddo
        enddo
      enddo

      return
 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine StiffnessMatrix(nbnodes,nodes,
     &                           nbelem,elements,
     &                           material,Ik,Jk,Vk)
      implicit none 

      integer nbnodes,nbelem
      double precision nodes(nbnodes,2)
      integer elements(nbelem,3)
      double precision CC(3,3)
      integer Ik(6*6*nbelem),Jk(6*6*nbelem)
      double precision Vk(6*6*nbelem)
      double precision thickness
      double precision material(3),young,nu

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) material
Cf2py intent(out) Ik
Cf2py intent(out) Jk
Cf2py intent(out) Vk

      double precision ke(6,6)
      double precision X(3),Y(3)
      integer idnodes(3),dofx(3),dofy(3),dofelem(6)
      integer p,e,i,j

      young     = material(1)
      nu        = material(2)
      thickness = material(3)
      CC(1,1)   = young/(1-nu**2)
      CC(1,2)   = nu*young/(1-nu**2)
      CC(1,3)   = 0.0d0
      CC(2,1)   = nu*young/(1-nu**2)
      CC(2,2)   = young/(1-nu**2)
      CC(2,3)   = 0.0d0
      CC(3,1)   = 0.0d0
      CC(3,2)   = 0.0d0
      CC(3,3)   = young/(2*(1+nu))

      p=1
      do e=1,nbelem

        do i=1,3
          idnodes(i) = elements(e,i)
c                          python indexing
          dofx(i)    = (idnodes(i)-1)*2
          dofy(i)    = (idnodes(i)-1)*2+1
        enddo
        
        dofelem(1) = dofx(1)
        dofelem(2) = dofx(2)
        dofelem(3) = dofx(3)
        dofelem(4) = dofy(1)
        dofelem(5) = dofy(2)
        dofelem(6) = dofy(3)

        do i=1,3
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
        enddo

        call ElementalStiffness(X,Y,thickness,CC,ke)

        do i=1,6
          do j=1,6
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
      subroutine ComputeElementalStress(X,Y,CC,Q,Sig,Eps,area)
      implicit none 

      double precision X(3),Y(3)
      double precision det_of_sys
      double precision a23(2,3)
      double precision a(3),b(3),c(3)
      double precision det33_ligne_de_un
      double precision BB(3,6),CC(3,3),CB(3,6)
      double precision Q(6),Sig(4),sigvm,Eps(3)
      double precision area
      
      integer i,j,k

      a23(1,1) = X(1)
      a23(1,2) = X(2)
      a23(1,3) = X(3)
      a23(2,1) = Y(1)
      a23(2,2) = Y(2)
      a23(2,3) = Y(3)
      det_of_sys=det33_ligne_de_un(a23)
      area=det_of_sys/2.0d0

      a(1)=+(X(2)*Y(3)-X(3)*Y(2))/det_of_sys
      b(1)=-(Y(3)-Y(2))/det_of_sys
      c(1)=+(X(3)-X(2))/det_of_sys

      a(2)=-(X(1)*Y(3)-X(3)*Y(1))/det_of_sys
      b(2)=+(Y(3)-Y(1))/det_of_sys
      c(2)=-(X(3)-X(1))/det_of_sys

      a(3)=+(X(1)*Y(2)-X(2)*Y(1))/det_of_sys
      b(3)=-(Y(2)-Y(1))/det_of_sys
      c(3)=+(X(2)-X(1))/det_of_sys

      BB(1,1)=b(1)
      BB(1,2)=b(2)
      BB(1,3)=b(3)
      BB(1,4)=0.0d0
      BB(1,5)=0.0d0
      BB(1,6)=0.0d0
      BB(2,1)=0.0d0
      BB(2,2)=0.0d0
      BB(2,3)=0.0d0
      BB(2,4)=c(1)
      BB(2,5)=c(2)
      BB(2,6)=c(3)
      BB(3,1)=c(1)
      BB(3,2)=c(2)
      BB(3,3)=c(3)
      BB(3,4)=b(1)
      BB(3,5)=b(2)
      BB(3,6)=b(3)

      do i=1,3
        do j=1,6
          CB(i,j)=0.0d0
          do k=1,3
            CB(i,j)=CB(i,j)+CC(i,k)*BB(k,j)
          enddo
        enddo
      enddo

      do i=1,3
        Sig(i)=0.0d0
        do k=1,6
          Sig(i)=Sig(i)+CB(i,k)*Q(k)
        enddo
      enddo

      do i=1,3
        Eps(i)=0.0d0
        do k=1,6
          Eps(i)=Eps(i)+BB(i,k)*Q(k)
        enddo
      enddo

      sigvm  = 1.5*(Sig(1)**2+Sig(2)**2
     &         +2.0*Sig(3)**2)
     &         -0.5*(Sig(1)+Sig(2))**2
      sigvm  = sqrt(sigvm)
      Sig(4) = sigvm


      return

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Compute_stress_strain_error(nbnodes,nodes,
     &                                       nbelem,elements,
     &                                       material,QQ,
     &                                       Sigma,
     &                                       sig_smooth,
     &                                       EpsilonElem,
     &                                       EpsilonNodes,
     &                                       ErrElem,
     &                                       ErrGlob)
      implicit none 

      integer nbnodes,nbelem
      double precision nodes(nbnodes,2)
      integer elements(nbelem,3)
      double precision CC(3,3),CCinv(3,3),QQ(2*nbnodes),Sigma(nbelem,4)
      double precision EpsilonElem(nbelem,3),sig_smooth(nbnodes,4)
      double precision ErrElem(nbelem),ErrGlob,EpsilonNodes(nbnodes,3)
      double precision material(3),young,nu,thickness

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

      double precision X(3),Y(3)
      double precision Q(6),Sig(4),Eps(3)
      integer e,i,j,g
      integer idnodes(3),dofx(3),dofy(3),dofelem(6)
      double precision rg(3),sg(3),tg(3),wg(3)
      double precision area,nodalweight(nbnodes)
      double precision NN(3),sig_smooth_Gauss(3),sigdiff(3)
      double precision NormSigElt(nbelem)
      double precision a23(2,3)
      double precision NormSig,det33_ligne_de_un
      double precision tmpeps(3),tmpsig(3)

      young     = material(1)
      nu        = material(2)
      thickness = material(3)
      CC(1,1)   = young/(1-nu**2)
      CC(1,2)   = nu*young/(1-nu**2)
      CC(1,3)   = 0.0d0
      CC(2,1)   = nu*young/(1-nu**2)
      CC(2,2)   = young/(1-nu**2)
      CC(2,3)   = 0.0d0
      CC(3,1)   = 0.0d0
      CC(3,2)   = 0.0d0
      CC(3,3)   = young/(2*(1+nu))

      CCinv(1,1)   = 1.0/young
      CCinv(1,2)   = -nu/young
      CCinv(1,3)   = 0.0d0
      CCinv(2,1)   = -nu/young
      CCinv(2,2)   = 1.0/young
      CCinv(2,3)   = 0.0d0
      CCinv(3,1)   = 0.0d0
      CCinv(3,2)   = 0.0d0
      CCinv(3,3)   = (2*(1+nu))/young
     

      rg(1)=2.0/3.0
      sg(1)=1.0/6.0
      tg(1)=1.0/6.0
      wg(1)=0.5*1.0/3.0

      rg(2)=1.0/6.0
      sg(2)=2.0/3.0
      tg(2)=1.0/6.0
      wg(2)=0.5*1.0/3.0

      rg(3)=1.0/6.0
      sg(3)=1.0/6.0
      tg(3)=2.0/3.0
      wg(3)=0.5*1.0/3.0

      do i=1,nbnodes
        nodalweight(i)=0.0d0
        do j=1,4
          sig_smooth(i,j)=0.0d0
        enddo
      enddo


      do e=1,nbelem
        do i=1,3
          idnodes(i) = elements(e,i)
c                          fortran indexing
          dofx(i)    = (idnodes(i)-1)*2+1
          dofy(i)    = (idnodes(i)-1)*2+2
        enddo
        
        dofelem(1) = dofx(1)
        dofelem(2) = dofx(2)
        dofelem(3) = dofx(3)
        dofelem(4) = dofy(1)
        dofelem(5) = dofy(2)
        dofelem(6) = dofy(3)

        do i=1,3
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
        enddo

        do i=1,6
          Q(i)=QQ(dofelem(i))
        enddo

        call ComputeElementalStress(X,Y,CC,Q,Sig,Eps,area)

        do i=1,4
          Sigma(e,i)=Sig(i)
        enddo

        do i=1,3
          nodalweight(idnodes(i))=nodalweight(idnodes(i))+area
          do j=1,4
            sig_smooth(idnodes(i),j)=sig_smooth(idnodes(i),j)
     &           +Sigma(e,j)*area
          enddo
        enddo

        do i=1,3
          EpsilonElem(e,i)=Eps(i)
        enddo

      enddo

c norm the smooth stress
      do i=1,nbnodes
        do j=1,4
          sig_smooth(i,j)=sig_smooth(i,j)/nodalweight(i)
        enddo
        do j=1,3
          tmpsig(j)=sig_smooth(i,j)
        enddo
        tmpeps=matmul(CCinv,tmpsig)
        do j=1,3
          EpsilonNodes(i,j)=tmpeps(j)
        enddo
      enddo

c compute error
      do e=1,nbelem
        ErrElem(e)=0.0d0
        NormSigElt(e)=0.0d0
      enddo
      NormSig=0.0d0

      do e=1,nbelem
        do i=1,3
          idnodes(i) = elements(e,i)
c                          fortran indexing
          dofx(i)    = (idnodes(i)-1)*2+1
          dofy(i)    = (idnodes(i)-1)*2+2
        enddo
        
        dofelem(1) = dofx(1)
        dofelem(2) = dofx(2)
        dofelem(3) = dofx(3)
        dofelem(4) = dofy(1)
        dofelem(5) = dofy(2)
        dofelem(6) = dofy(3)

        do i=1,3
          X(i)=nodes(idnodes(i),1)
          Y(i)=nodes(idnodes(i),2)
        enddo


        a23(1,1) = X(1)
        a23(1,2) = X(2)
        a23(1,3) = X(3)
        a23(2,1) = Y(1)
        a23(2,2) = Y(2)
        a23(2,3) = Y(3)
        area=det33_ligne_de_un(a23)/2.0d0

        do g=1,3

          NN(1)=rg(g)
          NN(2)=sg(g)
          NN(3)=tg(g)

          do j=1,3
            sig_smooth_Gauss(j)=0.0d0
            do i=1,3
              sig_smooth_Gauss(j)=sig_smooth_Gauss(j)
     &           +sig_smooth(idnodes(i),j)*NN(i)
            enddo
          enddo

          do j=1,3
            sigdiff(j)=sig_smooth_Gauss(j)-Sigma(e,j)
          enddo
          ErrElem(e)=ErrElem(e)
     &              +dot_product(sigdiff,matmul(CCinv,sigdiff))
     &                 *area/0.5d0*wg(g)

          ! compute norm of smooth solution in the element
          NormSigElt(e)=NormSigElt(e)
     &    +dot_product(sig_smooth_Gauss,matmul(CCinv,sig_smooth_Gauss))
     &                 *area/0.5d0*wg(g)

        enddo ! end Gauss point loop 
        NormSig = NormSig+NormSigElt(e)

      enddo ! end element loop

      ErrGlob=0.0d0
      do e=1,nbelem
        ErrElem(e) = ErrElem(e)/NormSig
        ErrGlob    = ErrGlob + ErrElem(e)
c        ErrElem(i) = sqrt(ErrElem(i))
      enddo

      ErrGlob    = sqrt(ErrGlob)

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine forceonline(nbnodes,
     &    nodes,
     &    nbelem,
     &    elements,
     &    fs,
     &    pts,
     &    F)


      implicit none

      integer nbnodes,nbelem
      double precision nodes(nbnodes,2)
      double precision fs(4)
      double precision pts(4)
      integer elements(nbelem,2)
      double precision F(2*nbnodes)

Cf2py intent(in) nbnodes
Cf2py intent(in) nodes
Cf2py intent(in) nbelem
Cf2py intent(in) elements
Cf2py intent(in) fs
Cf2py intent(in) pts
Cf2py intent(out) F

      double precision rg(2)
      double precision wg(2)
      double precision xnodes(2)
      double precision ynodes(2)
      double precision lelem
      double precision N(2)
      double precision forceelem(4)
      double precision x
      double precision y
      double precision l0
      double precision l1
      double precision Phi(2,4)
      double precision Phit(4,2)
      double precision forcegausspt(2)
      integer i
      integer e
      integer g
      integer idnodes(2)
      integer dofx(2)
      integer dofy(2)
      double precision forcex
      double precision forcey

      rg(1)=-1.0/dsqrt(3.0d0)
      rg(2)=1.0/dsqrt(3.0d0)
      wg(1)=1.0
      wg(2)=1.0

      do i=1,2*nbnodes
        F(i)=0.0d0
      enddo

      do e=1,nbelem

        idnodes(1)=elements(e,1)
        idnodes(2)=elements(e,2)
        dofx(1)=(idnodes(1)-1)*2+1
        dofy(1)=(idnodes(1)-1)*2+2
        dofx(2)=(idnodes(2)-1)*2+1
        dofy(2)=(idnodes(2)-1)*2+2

        xnodes(1)=nodes(idnodes(1),1)
        xnodes(2)=nodes(idnodes(2),1)
        ynodes(1)=nodes(idnodes(1),2)
        ynodes(2)=nodes(idnodes(2),2)

        lelem=dsqrt((xnodes(2)-xnodes(1))**2+(ynodes(2)-ynodes(1))**2)

        do i=1,4
          forceelem(i)=0.0d0
        enddo

        do g=1,2
          N(1) = (1-rg(g))*0.5d0
          N(2) = (1+rg(g))*0.5d0
          x = N(1)*xnodes(1)+N(2)*xnodes(2)
          y = N(1)*ynodes(1)+N(2)*ynodes(2)
          !pts=[ pt 1 x, pt 1 y , pt 2 x , pt 2 y]
          l0=dsqrt((pts(1)-x)**2+(pts(2)-y)**2)
          l1=dsqrt((pts(3)-x)**2+(pts(4)-y)**2)

          ! fs=[  fs-1-x , fs-1-y , fs-2-x  , fs-2-y   ]
          forcex=(l1*fs(1)+l0*fs(3))/(l0+l1)
          forcey=(l1*fs(2)+l0*fs(4))/(l0+l1)
          
          Phi(1,1)=N(1)
          Phi(1,2)=0.0d0
          Phi(1,3)=N(2)
          Phi(1,4)=0.0d0
          Phi(2,1)=0.0d0
          Phi(2,2)=N(1)
          Phi(2,3)=0.0d0
          Phi(2,4)=N(2)
          Phit = transpose(Phi)

          forcegausspt(1)=forcex
          forcegausspt(2)=forcey
          
          forceelem=forceelem
     &       +matmul(Phit,forcegausspt)*wg(g)*0.5d0*lelem

        enddo
        F(dofx(1))=F(dofx(1))+forceelem(1)
        F(dofy(1))=F(dofy(1))+forceelem(2)
        F(dofx(2))=F(dofx(2))+forceelem(3)
        F(dofy(2))=F(dofy(2))+forceelem(4)
      enddo
      

      return
      end