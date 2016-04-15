********************************************************************* 
* any question eduardor@fisica.unam.mx                              *
*********************************************************************
	subroutine kernint(mn,sigmai)
	implicit none
        include 'common.f'
        integer i,j,ip,iq,iqp,i1,i2
	double precision ff(3),w1
        double precision pi,zq,jac,pq,k2,q2,p2,
     .  q,zp,zqp,qv(4),pv(4),u,g2,sigmai(2),mn,kv(4)

        double complex suma,sigmasp,
     .  sigmasm,sigmavp,sigmavm,
     .  Pcv(4),tr,pc4

     
        external u
        external g2
        external kabt

        pi = 3.14159265358979323846D0
        

        w1= mh/(mh+mh2)
          

*==============================
*     kernel 
*==============================


        
     
        q2 = q2E
                    
       if(init.eqv..true.)then
 
        suma = 0d0

       
      
       do 405 i1=1,ns
       do 406 i2=1,ns

       do 402 ip =1,nt
       do 403 iq =1,nt
       do 404 iqp=1,nzq

       if(ffit.eqv..true.)then 
       p =     (1d0+xs(i1))/(1d0-xs(i1))
       q2 = (  (1d0+xs(i2))/(1d0-xs(i2)) )**2
       else
       p=   (lambda-dsqrt(pmin2))*xs(i1)+dsqrt(pmin2)           
       q2= ((lambda-dsqrt(pmin2))*xs(i2)+dsqrt(pmin2))**2
       endif
      

      
       zp = xt(ip)
       zq = xt(iq)
       zqp= xq(iqp)

       
       jac = 2d0*dsqrt(1d0-zp**2)*dsqrt(1d0-zq**2) !*u(m,zp)*u(n,zq)
     .      *wt(ip)*wt(iq)*wq(iqp)
            

 

 
       
       q = abs(dsqrt(q2))
       p2 = p*p
            

       do 102 i=1,4
       PCV(i) =0d0
102    continue
       PCV(4)= (0d0,1d0)*mn

       pv(1)= 0d0
       pv(2)= 0d0
       pv(3)= dsqrt(p2)*dsqrt(1d0-zp*zp)
       pv(4)= dsqrt(p2)*zp


       qv(1)= q*dsqrt(1d0-zq*zq)*dsqrt(1d0-zqp*zqp)/dsqrt(2d0)
       qv(2)= q*dsqrt(1d0-zq*zq)*dsqrt(1d0-zqp*zqp)/dsqrt(2d0)
       qv(3)= q*dsqrt(1d0-zq*zq)*zqp
       qv(4)= q*zq

       k2=0d0
       do 101 i =1,4
        kv(i)= pv(i)-qv(i)
        k2= kv(i)**2+k2
101    continue
  
          
       



       sigmavp = ssigmavp(i2,iq)
       sigmasp = ssigmasp(i2,iq)
       sigmavm = ssigmavm(i2,iq) 
       sigmasm = ssigmasm(i2,iq)
      
     

       call  kabT (printkabt,w1,pcv,pv,qv,zp,zq,kv,sigmavp,sigmasp,
     . sigmavm,sigmasm,g2(k2,z2k),alpha,beta,tr)

        kabR(i1,i2,ip,iq,iqp)=tr*jac  
        
 404    continue
 403    continue
 402    continue

 406    continue
 405    continue
         
        init = .false.
        write(*,*)"init ready"
        else
        endif    

        suma=0d0

        do 501 ip =1,nt
        do 502 iq =1,nt
        do 503 iqp=1,nzq
        zp = xt(ip)
        zq = xt(iq)


        suma=u(m,zp)*kabR(jp,jq,ip,iq,iqp)*u(n,zq)+suma
      
       
503     continue
502     continue
501     continue

c        endif


        sigmai(1) =   dble(suma)/(2d0*pi)**4 
        sigmai(2) =   dimag(suma)/(2d0*pi)**4

        

        return  
	end

