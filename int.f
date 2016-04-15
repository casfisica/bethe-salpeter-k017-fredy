* demo-fortran.F.
*********************************************************************
* This program has adapted the Cuba1.5 Integrator.                  * 
* any question rojas4000@yahoo.com                                  *
*********************************************************************
* test program for the Cuba library



	subroutine int(print2,mn2,ppf2,FR2,FI2,sigmai)
	implicit none         
        integer ndim, ncomp, mineval, maxeval, verbose, last
        integer introut 
        double precision epsrel, epsabs, sigmai(2)
        double precision pi
        logical print,print2
       	parameter (pi = 3.14159265358979323846D0)
    
        double precision ppf2(64),FR2(4,0:5,64),ppf(64),FR(4,0:5,64),
     .  mn2
        double precision FI2(4,0:5,64),FI(4,0:5,64),mn
        integer i,j,k
        common /int1/ppf(64),FR(4,0:5,64),FI(4,0:5,64),mn
        common /int3/print

             
	parameter (ndim = 2)
	parameter (ncomp = 2)
	parameter (epsrel = 1D-8)
	parameter (epsabs = 1D-12)
	parameter (verbose = 2)
	parameter (last = 4)
	parameter (mineval = 0)
	parameter (maxeval = 50000)

	integer nstart, nincrease
	parameter (nstart = 10000)
	parameter (nincrease = 5000)

	integer nnew
	double precision flatness
	parameter (nnew = 100000)
	parameter (flatness = 1D0)

	integer key1, key2, key3, maxpass
	double precision border, maxchisq, mindeviation
	integer ngiven, ldxgiven, nextra
	parameter (key1 = 47)
	parameter (key2 = 1)
	parameter (key3 = 1)
	parameter (maxpass = 5)
	parameter (border = 0D0)
	parameter (maxchisq = 10D0)
	parameter (mindeviation = .25D0)
	parameter (ngiven = 0)
	parameter (ldxgiven = ndim)
	parameter (nextra = 0)

	integer key
	parameter (key = 0)

	external integrand3

	double precision integral(ncomp), error(ncomp), prob(ncomp)
	integer nregions, neval, fail

	integer c               
**************************************************
*       Introut = integration soubroutine        *
*       1 vegas, 2 suave, 3 divonne, 4 Cuhre     *
**************************************************

        introut = 1


C	 "------Integration subroutines--------------------"

           
      do 4 i = 1,8
      do 5 j = 0,5        !nt
      do 6 k=1,64
      FR(i,j,k)= FR2(i,j,k) 
      FI(i,j,k)= FI2(i,j,k)
6     continue 
5     continue
4     continue

 
      do 7 k=1,64
      ppf(k)=ppf2(k)
7     continue
    
      mn = mn2
      print = print2


        if(introut.eq.1) then
           call vegas(ndim, ncomp, integrand3,
     &    epsrel, epsabs, verbose, mineval, maxeval,
     &    nstart, nincrease,
     &    neval, fail, integral, error, prob)

        elseif(introut.eq.2) then
          call suave(ndim, ncomp, integrand3,
     &    epsrel, epsabs, verbose + last, mineval, maxeval,
     &    nnew, flatness,
     &    nregions, neval, fail, integral, error, prob)
        
        elseif(introut.eq.3) then
          call divonne(ndim, ncomp, integrand3,
     &    epsrel, epsabs, verbose, mineval, maxeval,
     &    key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    nregions, neval, fail, integral, error, prob)
        
         elseif(introut.eq.4) then
        call cuhre(ndim, ncomp, integrand3,
     &    epsrel, epsabs, verbose + last, mineval, maxeval,
     &    key,
     &    nregions, neval, fail, integral, error, prob)
        endif
          
        sigmai(1) =  integral(1) 
        sigmai(2) =  integral(2)

	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &  (integral(c), error(c), prob(c), c = 1, ncomp)

              


        

 
    
       return

	end


************************************************************************

	subroutine integrand3(ndim, xx, ncomp, ff)
	implicit none
        include 'commonf.f'
        integer ndim, ncomp 
	double precision xx(*), ff(*),qq,zq,jac,qv(4),q
    
        double precision pi,x,ppf(64),FR(4,0:5,64),FI(4,0:5,64),
     .  mn,w1,w2,lamb_bs,pmin2_
        integeri,j,k
        double complex pc4,pcv(4),traza1,dAqp,dBqp,dAqm,dBqm,
     .  svp,ssp,svm,ssm,dsvp,dssp,dsvm,dssm,F(4,0:5)

        logical print
        common /int1/ppf(64),FR(4,0:5,64),FI(4,0:5,64),mn 
        common /int3/print
        pi = 3.14159265358979323846D0
                      

*=============================================================================
*     ESTE PROGRAMA DEFINE EL INTEGRANDO PARA CUBA
*=============================================================================

             
       pmin2_ =pmin2

       qq = (lambda**2-pmin2_)*xx(1)+pmin2_
       zq = 2*xx(2)-1d0

       q = abs(dsqrt(qq))
       


       do 102 i=1,4
       PCV(i) =0d0
102    continue
       PCV(4)= (0d0,1d0)*mn

     


       qv(1)=0d0 ! q*dsqrt(1d0-zq*zq)*dsqrt(1d0-zqp*zqp)/dsqrt(2d0)
       qv(2)=0d0 ! q*dsqrt(1d0-zq*zq)*dsqrt(1d0-zqp*zqp)/dsqrt(2d0)
       qv(3)= q*dsqrt(1d0-zq*zq)   !*zqp
       qv(4)= q*zq

       ! d^4q = q^2dq^2/2d0*dsqrt(1d0-zq^2)*dzqp*dphi
       !dqq = (lambda**2-pmin2)
       !dzq = 2d0*dsqrt(1d0-zq**2)
       !dzqp = 2d0
       !dphi = 2d0*pi
       !jac = qq/2d0*dqq*dzq*dzqp*dphi
        jac = qq/2d0*(lambda**2-pmin2_)*2d0*dsqrt(1d0-zq**2)*
     .        2d0*2d0*pi  

      lamb_bs = lambda

      call interT(64,ppf,FR,FI,q,F)

      if(ffit.eqv..true.)then                 ! nakanishi
      lambda =  dsqrt(lamb_bs**2+mh**2/4d0)
      w1 = mh/(mh+mh2)
      pc4 =  (0d0,1d0)*mn*w1
c     CALL    dscomplex(print,mn,pc4,qq,zq,svp,ssp,dsvp,dssp)
      call dscomplexfit(print,mn,pc4,qq,zq,xval,npol,svp,ssp,
     . dsvp,dssp)
      w2 = mh2/(mh+mh2)
      pc4 = -(0d0,1d0)*mn*w2
      lambda =  dsqrt(lamb_bs**2+mh2**2/4d0)
c     call   dscomplex2(print,mn,pc4,qq,zq,svm,ssm,dsvm,dssm)
      call dscomplexfit(print,mn,pc4,qq,zq,xval2,npol2,svm,ssm,
     . dsvm,dssm)
      else                                      !  SDE in complex plane
      lambda =  dsqrt(lamb_bs**2+mh**2/4d0)
      w1 = mh/(mh+mh2)
      pc4 =  (0d0,1d0)*mn*w1
      CALL dscomplex(print,mn,pc4,qq,zq,svp,ssp,dsvp,dssp)
      w2 = mh2/(mh+mh2)
      pc4 = -(0d0,1d0)*mn*w2
      lambda =  dsqrt(lamb_bs**2+mh2**2/4d0) 
      call dscomplex2(print,mn,pc4,qq,zq,svm,ssm,dsvm,dssm)
      endif
      lambda = lamb_bs


       if(print.eqv..true.)then
       write(2,*) "int variables ____________________________________"
       write(2,*) "xx(1),xx(2),qq,zq,jac",xx(1),xx(2),qq,zq,jac
       write(2,*)" qv(1),qv(2),qv(3),qv(4)",qv(1),qv(2),qv(3),qv(4)
       write(2,*) "q,lambda",q,lambda
       write(2,*)"F"
       do  103 i =1,4
       write(2,*)(F(i,j),j=0,4)
103    continue
       write(2,*) "pi,x,mn",pi,x,mn
       write(2,*)"ppf"
       write(2,*)(ppf(i),i=1,64)
       write(2,*)"FR"
       do  104 i=1,4
       write(2,*)(FR(i,j,1),j=0,4)
104    continue
       write(2,*)"FI"
       do  105 i=1,4
       write(2,*)(FI(i,j,1),j=0,4)
105    continue
       write(2,*)" pc4,pcv(4),traza1",pc4,pcv(4),traza1
       write(2,*)" dAqp,dBqp,dAqm,dBqm",dAqp,dBqp,dAqm,dBqm
       write(2,*)" svp,ssp,svm,ssm",svp,ssp,svm,ssm
       write(2,*)" dsvp,dssp,dsvm,dssm",dsvp,dssp,dsvm,dssm
       write(2,*)"pcv(1),pcv(2),pcv(3),pcv(4)",pcv(1),pcv(2),
     . pcv(3),pcv(4)
       else
       endif

     
 
      if(ntr.eq.1)then
      call traza(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,dssp,dssm,
     . dsvp,dsvm,F,traza1)
       ff(1)  =  jac*dble(traza1)/(2d0*pi)**4 
       ff(2)  =  jac*dimag(traza1)/(2d0*pi)**4
      elseif(ntr.eq.2)then
      call trazad(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,F,traza1)
      ff(1)  =  jac*dimag(traza1)/(2d0*pi)**4
      ff(2)  =  jac*dble(traza1)/(2d0*pi)**4
      elseif(ntr.eq.3)then
      call trazarh(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,F,traza1)
      ff(1)  =  jac*dimag(traza1)/(2d0*pi)**4
      ff(2)  =  jac*dble(traza1)/(2d0*pi)**4
      elseif(ntr.eq.4)then
      call traza2(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,dssp,dssm,
     . dsvp,dsvm,F,traza1)
       ff(1)  =  jac*dble(traza1)/(2d0*pi)**4
       ff(2)  =  jac*dimag(traza1)/(2d0*pi)**4
      else
      endif




       

200    format(D25.16,D25.16,D25.16)

	end




