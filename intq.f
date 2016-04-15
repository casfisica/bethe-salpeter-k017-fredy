*********************************************************************
* any question rojas4000@yahoo.com                                  *
*********************************************************************
* test program for the Cuba library



	subroutine intq(print2,mn2,ppf2,FR2,FI2,sigmai)
	implicit none         
        double precision  sigmai(2)
        double precision pi
        logical print,print2
       	parameter (pi = 3.14159265358979323846D0)
        double precision ppf2(64),FR2(8,0:5,64),ppf(64),FR(8,0:5,64),
     .  mn2,lamb_bs
        double precision FI2(8,0:5,64),FI(8,0:5,64),mn
        integer i,j,k
        include 'commonf.f' 
	double precision qq,zq,jac,qv(4),q
        double precision x,w1,w2,dq,dzq,ff(2)
        integer is,it
        double complex pc4,pcv(4),traza1,dAqp,dBqp,dAqm,dBqm,
     .  svp,ssp,svm,ssm,dsvp,dssp,dsvm,dssm,F(8,0:5)

      ff(1)=0d0
      ff(2)=0d0
                                            
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
      do 102 i=1,4
       PCV(i) =0d0
102    continue
       PCV(4)= (0d0,1d0)*mn

      print = print2

       do 106 is =1,ns
       do 107 it =1,nt 
        
       q = (1d0+xs(is))/(1d0-xs(is))    !(lambda-dsqrt(pmin2))*xs(is)+dsqrt(pmin2)
       zq = xt(it)
       qq = q*q
       dq = 2d0/(1d0-xs(is))**2*ws(is)   !(lambda-dsqrt(pmin2))*ws(is)
       dzq= wt(it)
       
       
       qv(1)=0d0 ! q*dsqrt(1d0-zq*zq)*dsqrt(1d0-zqp*zqp)/dsqrt(2d0)
       qv(2)=0d0 ! q*dsqrt(1d0-zq*zq)*dsqrt(1d0-zqp*zqp)/dsqrt(2d0)
       qv(3)= q*dsqrt(1d0-zq*zq)   !*zqp
       qv(4)= q*zq

       !d^4q = q^3dsqrt(1d0-zq^2)*dq*dzq*dzqp*dphi
       jac   = q**3*dsqrt(1d0-zq**2)*dq*dzq*2d0*2d0*pi  

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
       write(2,*) "xs(1),xt(2),qq,zq,jac",xs(1),xt(2),qq,zq,jac
       write(2,*)" qv(1),qv(2),qv(3),qv(4)",qv(1),qv(2),qv(3),qv(4)
       write(2,*) "q",q
       write(2,*)"F"
       do  103 i =1,8
       write(2,*)(F(i,j),j=0,4)
103    continue
       write(2,*) "pi,x,mn",pi,x,mn
       write(2,*)"ppf"
       write(2,*)(ppf(i),i=1,64)
       write(2,*)"FR"
       do  104 i=1,8
       write(2,*)(FR(i,j,1),j=0,4)
104    continue
       write(2,*)"FI"
       do  105 i=1,8
       write(2,*)(FI(i,j,1),j=0,4)
105    continue
       write(2,*)" pc4,pcv(4),traza1",pc4,pcv(4),traza1
       write(2,*)" dAqp,dBqp,dAqm,dBqm",dAqp,dBqp,dAqm,dBqm
       write(2,*)" svp,ssp,svm,ssm",svp,ssp,svm,ssm
c       write(2,*)" dsvp,dssp,dsvm,dssm",dsvp,dssp,dsvm,dssm
c       write(2,*)"pcv(1),pcv(2),pcv(3),pcv(4)",pcv(1),pcv(2),
c     . pcv(3),pcv(4)
       else
       endif

     
 
      if(ntr.eq.1)then
      call traza(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,dssp,dssm,
     . dsvp,dsvm,F,traza1)
       ff(1)  =  jac*dble(traza1)/(2d0*pi)**4+ff(1) 
       ff(2)  =  jac*dimag(traza1)/(2d0*pi)**4+ff(2)
      elseif(ntr.eq.2)then
      call trazad(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,F,traza1)
      ff(1)  =  jac*dimag(traza1)/(2d0*pi)**4+ff(1)
      ff(2)  =  jac*dble(traza1)/(2d0*pi)**4+ff(2)
      elseif(ntr.eq.3)then
      call trazarh(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,F,traza1)
      ff(1)  =  jac*dimag(traza1)/(2d0*pi)**4+ff(1)
      ff(2)  =  jac*dble(traza1)/(2d0*pi)**4+ff(2)
      elseif(ntr.eq.4)then
      call traza2(print,w1,zq,qv,pcv,ssp,ssm,svp,svm,dssp,dssm,
     . dsvp,dsvm,F,traza1)
       ff(1)  =  jac*dble(traza1)/(2d0*pi)**4+ff(1)
       ff(2)  =  jac*dimag(traza1)/(2d0*pi)**4+ff(2)
      else
      endif

107   continue
106   continue
      
      sigmai(1)=ff(1)
      sigmai(2)=ff(2)      


       

200    format(D25.16,D25.16,D25.16)

	end




