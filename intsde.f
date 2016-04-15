* demo-fortran.F.
*********************************************************************
* This program has adapted the Cuba1.5 Integrator  to Calculate     *
* Sigmai. Any cuestion eduardor@fisica.unam.mx                      *
********************************************************************* 
*To change some                                                     *
*parameters see below. 31 october 2010                              *
*********************************************************************
* test program for the Cuba library
* last modified 2 Feb 05 th


	subroutine intsde(sigmai)
        implicit none
        include 'commonsde.f'
        integer ndim, ncomp,i,j
	double precision xx(2), ff(2)
        double precision z1,jac,cf,q2,pr2,mathcalg,kr,g2,k2r,
     .  p2r,pkr
        double complex A2,B2,A,B,k2,p2,tr
        double precision sigmai(2)

        external vertex

        pi = 3.14159265358979323846D0
        
*=============================================================================
*     ESTE PROGRAMA DEFINE EL INTEGRANDO PARA CUBA
*=============================================================================

      ff(1)=0d0
      ff(2)=0d0

      do 100 i=1,np
      do 101 j=1,nzq
    
      kr= (lambda-dsqrt(pmin2))*xp(i)+dsqrt(pmin2)
       
c       k2r = (lambda**2-pmin2)*xx(1)+pmin2
       k2r=kr*kr
       z1 = 2d0*xzq(j)-1d0

       jac =2*pi*2d0*(lambda-dsqrt(pmin2))*kr**3*2*dsqrt(1d0-z1**2)*
     .  wp(i)*wzq(j) 
       cf = 4d0/3d0

c       kr = abs(dsqrt(k2r))
       p2r = pr*pr
       pkr = dsqrt(p2r*k2r)*z1
       q2 = p2r+k2r-2*pkr


      
        mathcalg = g2(q2,z2g)

        
        call vertex(k2r,q2,z1,mathcalg,tr)
        ff(1)  = dble( tr/(2*pi)**4*jac)+ff(1)
        ff(2)  = dimag(tr/(2*pi)**4*jac)+ff(2)

101     continue
100     continue


       sigmai(1) =  ff(1)
       sigmai(2) =  ff(2)
	end



