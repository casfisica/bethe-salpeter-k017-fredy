
        double precision function g2(q2,z2)
        implicit none 
        double precision  mt2,D,Nf,gammam,
     .  lambdaqcd,omega,tau,mathcalf,q2,LPV
        double precision pi,zz1,m1,m2,m3
        double precision  CA,gone2,rho1,mu2,rho2,z2
        parameter (pi = 3.14159265358979323846D0)
   
 

c----------------------------mathcalg calculation ------------------------------------------
       mt2 = (0.5d0)**2
c       D = 0.93d0 
       Nf = 4d0
       gammam = 12d0/(33d0-2*Nf)
       lambdaqcd = 0.234d0    ! in GeV
       omega =0.4d0  ! 0.6d0 ! 0.5d0
       
   
      
       tau = (2.718281828d0)**2-1d0 
       mathcalf = (1d0-dexp(-q2/(4*mt2)))/q2

c       if(dsqrt(q2).lt.4d0)then
c       call alphastrong(4d0,alphas)
c       elseif(dsqrt(q2).ge.4d0)then
c       call alphastrong(dsqrt(q2),alphas)
c       endif

       D= 0.372d0/omega

       g2 = q2*(4*pi**2*D/(omega)**6*q2*dexp(-q2/omega**2)    ! *q2  MT
     . +4*pi*gammam*pi
     . /(  1d0/2d0*dlog(tau+(1d0+q2/lambdaqcd**2)**2)  )      ! q2( 1/(log())
     . *mathcalF)
      
c         D= (0.87d0)**3/omega
c        D= (0.8d0)**3/omega
c        D= (1.1d0)**3/omega
c         D= (0.55d0)**3/omega 


c       g2 =z2**2*q2*(8*pi**2*D/(omega)**4*dexp(-q2/omega**2)    ! *q2 qin chang liu roberts wilson 2011
c     . +8*pi**2*gammam
c     . /( dlog(tau+(1d0+q2/lambdaqcd**2)**2)  )      ! q2( 1/(log())
c     . *mathcalF)
         


c---------------------------------------------------------------------------------------------- 
c       ZZ1 =0.8290  !0.83333d0       !0.821d0  ! 0.8196d0 
c       M1 = 3.94d0   !4.473d0       !4.09d0  ! 4.22d0
c       M2 = 0.583d0  !0.704d0       !0.558d0  ! 0.631d0
c       M3 = 0.3224d0  !0.3959d0       !0.380d0  ! 0.3177d0
c        g2 =4d0*pi*q2*(q2+M1)
c     . *0.295d0*ZZ1/(q2**2+M2*q2+M3)		    


        
c       CA = 3d0
c       mu2 = (4.3d0)**2
c       rho1 = 8.55d0
c       rho2 = 1.91d0
c       m2 = (0.52d0)**4/(q2+rho2*(0.52d0)**2)
c       gone2 = 5.68d0

c       g2 = 1d0/ (
c     . m2+q2*(1d0+  13d0/(96d0*pi**2)*CA*gone2
c     .  *dlog(   (q2+rho1*m2)/mu2) 
c     .      )    )
c     . *4d0*pi*q2*0.295d0

         lpv = 2000d0
         g2 = g2*1d0/(1d0+Q2/LPV**2)

        return
        end
