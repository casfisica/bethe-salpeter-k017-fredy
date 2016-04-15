       subroutine  vertexr(flag1,nvertex,ikernel,np,pp,aain,
     .            bbin,ccin,p2,lambda,pmin2,z2,k2,z1,tr)
       implicit none
c       include 'common.f'
       double precision La(4),Lb(4),tr,k2,kp,Ak2,Bk2,Ak,Bk,mathcalg,
     . Ta(8),Tb(8),d,suma,llambda(4),ttau(8),delta2,cf,Bp,Ap,Bp2,Ap2,
     . Au,Au2,Bu,Bu2,u,v,u2,v2,CA,Fgh,g2,Cq,Ck,chi0
     
       double precision q2,z1,h1,p,k,lpv
       integer  i,j
       

        double precision p2,lambda,pmin2,z2
        integer np,nvertex,ikernel
        double precision pp(128),aain(128),bbin(128),ccin(128)
        logical flag1

        



       external massf
c       external chi0
       external h1
c       external inter

       cf = 4d0/3d0
 
       p=dsqrt(p2)
       kp = dsqrt(k2*p2)*z1
       delta2 = kp**2-k2*p2
       q2= p2+k2-2*kp       
       u= dsqrt(p2+k2+2*kp)
       v= dsqrt(p2/4d0+k2+2*kp/2d0)
       u2= u*u
       v2= v*v
       k = abs(dsqrt(k2))

       call massf(np,pp,BBin,p,Bp)
c       call inter(np,pp,BBin,p,FAK1p,FAK2p,FAk3p,Bp)
        
       call massf(np,pp,AAin,p,Ap)
       call massf(np,pp,BBin,k,Bk)
       call massf(np,pp,AAin,k,Ak)
 
c       call massf(np,pp,AAin,u,Au)
c       call massf(np,pp,BBin,u,Bu) 
c       call massf(np2,ppmil,CC,dsqrt(q2),Cq)
c       call massf(npc,ppc,CCc,dsqrt(q2),Cq)   ! for a precualculated x0
c       call massf(np,pp,CCin,dsqrt(q2),Cq)
c        cq = chi0(q2)
       if(ikernel.eq.3)then
       call massf(np,pp,CCin,dsqrt(k2),Ck)
       else
       endif


       Bk2 = Bk*Bk
       Ak2 = Ak*Ak
       Bp2 = Bp*Bp
       Ap2 = Ap*Ap
       Au2= Au*Au
       Bu2= Bu*Bu
      
       d = 4d0




       do 102 i = 1,3
       llambda(i) = 0d0
102    continue

       do 103 i = 2,8
       ttau(i) = 0d0
103    continue

c-------------------------- vertex----------------

        if(nvertex.eq.1)then
  
        llambda(1) = 1d0
        elseif(nvertex.eq.2)then

        llambda(1) = (Ak+Ap)/2d0
        llambda(2) = (Ak-Ap)/(k2-p2+1d-10)/2d0  ! lambdaE
        llambda(3) =  -(Bk-Bp)/(k2-p2+1d-10)   
        elseif(nvertex.eq.3)then
        cq = cq*Fgh(q2)
        ck = ck*Fgh(k2)
        llambda(1) = Cq*Fgh(q2)*(Ak+Ap)/2d0
        llambda(2) = Cq*Fgh(q2)*(Ak-Ap)/(k2-p2)/2d0  ! lambdaE
        llambda(3) =  -Cq*Fgh(q2)*(Bk-Bp)/(k2-p2)
        else
        endif


       Lb(1) = (d-1)*Bk
       Lb(2) = -4d0*Bk*delta2/q2
       Lb(3) = -2d0*Ak*delta2/q2
       Tb(2) = -2d0*Bk*delta2
       Tb(3) = -(d-1)*Bk*q2
       Tb(6) =  (d-1)*Bk*(k2-p2)


        
       La(1) = (2d0-d)*q2/2d0+(d-3d0)*(k2+p2)/2+(k2-p2)**2/(2*q2)      
c       La(1) = 3*kp+2d0*(delta2)/q2  !pk
       La(2) = 2*(k2+p2)*Delta2/q2 !  RW I 
c       La(2) = -2*(k2+p2)*Delta2/q2  ! pk
       La(3) = -2*Bk*delta2/(Ak*q2)  ! RW
c       Ta(2) = (k2+p2)*delta2
        Ta(2) = (k2+p2)*delta2  
c       Ta(3) =-(-(d-2)*(q2)**2/2d0+(d-3)*(k2+p2)*q2/2d0
c     .        +(k2-p2)**2/2d0)
        Ta(3) = -(2*delta2+3*q2*kp)
c       Ta(6) = (d-1)*kp*(k2-p2)
        Ta(6) = -(-3d0*(k2-p2)*kp)
c       Ta(8) = -(d-2)*delta2 
       Ta(8) = -(2*delta2)
                  
             
        mathcalg = g2(q2,z2)

c---------------------------------------------------------------------------------------------------(ikernel 1)
        LPV = 2000d0

       if(ikernel.eq.1)then 
       
       suma = 0d0

       do 100 i = 1,3
       suma = Lb(i)*llambda(i)+suma
 100   continue

        
       tr = cf*mathcalg*1d0/q2*                        !1d0/(1d0+q2/LPV**2)*
     .      1d0/(Ak2*k2+Bk2)*suma  ! mathcalg = 4*pi*alpha = g^2

c---------------------------------------------------------------------------------------- (ikernel 2)
       elseif(ikernel.eq.2)then

     
       suma = 0d0

       do 104 i = 1,3
       suma = La(i)*llambda(i)+suma
 104   continue

       
       tr = (cf/p2)*mathcalg*1d0/q2*          !1d0/(1d0+q2/LPV**2)*
     .      Ak/(Ak2*k2+Bk2)*suma  ! mathcalg = 4*pi*alpha
c---------------------------------------------------------------------------------------------(ikernel 3)
       elseif(ikernel.eq.3)then
        
       CA = 3d0 

       tr = CA*g2(k2,z2)/(8d0*k2)*(-delta2)/k2*Fgh(v2)/v2*Fgh(k2)    !**2
     . *Au*(Au+Ap)/(Au2*u2+Bu2)*Ck                   !*h1(u2)                                              

c      tr = CA*g2(k2,z2)/(8d0*k2)*(-delta2)/k2*Fgh(v2)/v2*Fgh(k2)
c     . *2d0/(u2)          

c      tr = CA*g2(k2,z2)/(4d0*k2)*(-delta2)/k2*Fgh(u2)*Fgh(k2)**2   !3.17  
c     . *1d0/(u2)**2



       else
       endif

       

    
         


       
      
       if((p2.eq.pmin2).and.(flag1.eqv..true.))then  
       write(2,*)'tr,suma,mathcalg,k2,Bk,Ak',tr,suma, 
     .  mathcalg,k2,Bk,Ak
       write(2,*)'1/k2,,p2,1d0/p2,q2,z1',1d0/k2,p2,1d0/p2,q2,z1
       write(2,*)'La(1),llambda(1)',La(1),llambda(1)


       write(2,*) 'suma', suma
       flag1 = .false.
       endif

    

       return
       end




       double precision function chi0(x)
       ! x = q^2
       double precision x,xpar1(9)
       data  xpar1/0.5391848796769366d0,
     .  1.0733422786909068d0, 7.698064712961662d0,
     .  1.3501603716181783d0, 2.327608552776829d0,
     .  5.77041495082685d0, 2.681835086130385d0,
     .  0.005945186487675137d0, 0.012694192214748047d0/


       chi0 =   xpar1(1) + (xpar1(2) + x*xpar1(3) + x**2*xpar1(4)
     .  + x**3*xpar1(8))/(xpar1(5) + x*xpar1(6) + x**2*xpar1(7)
     .  + x**3*xpar1(9))
       return
       end



  

