       subroutine vertex(k2r,q2,z1,mathcalg,tr)
       implicit none
       include 'commonsde.f'
       include 'commoncmaker.f'
       double precision mathcalg,kpr,
     . d,cf, Au,Au2,Bu,Bu2,u,v,u2,v2,CA,Fgh,g2,Cq,Ck,chi0,z(3),dz(3)
       double precision q2,z1,h1,q,kr,k2r,p2r,krv(4),prv(4)
       double complex p2,pIv(4),k2,kp,ak2,bk2,ak,bk,delta2,bp,ap,bp2,
     . ap2,sumav(6),llambda(4),La(4),Lb(4),ttau(8),Ta(8),Tb(8),tr,
     .  kx2(3),dkx2(3),suma
       integer  i,j,li

       logical flag1       

       external massf
c       external chi0
       external h1
c       external inter

       cf = 4d0/3d0
       p2 = p*p

       kr = abs(dsqrt(k2r))
       p2r = pr*pr
       q = dsqrt(q2)
       kpr = dsqrt(p2r*k2r)*z1

       krv(3) =kr*dsqrt(1d0-z1**2)
       krv(4) =kr*z1
       prv(3) =0d0
       prv(4) =pr  
       pIv(3) =0d0
       pIv(4) =pc4

c       k2 =0d0
c       do 101 i=3,4
c       k2=(krv(i)+pIv(i))*(krv(i)+pIv(i))+k2
c101    continue
      
       k2 = (krv(3)+pIv(3))*(krv(3)+pIv(3))
     .     +(krv(4)+pIv(4))*(krv(4)+pIv(4))

       kp = (krv(4)+pIv(4))*(prv(4)+pIv(4)) 

c       kp = zsqrt(k2*p2)*z1
       delta2 = kp**2-k2*p2
c        delta2 =(kpr)**2-k2r*p2r     

c       call interij(np,contour,pp,BBin,kr,Bk)
c       call interij(np,contour,pp,AAin,kr,Ak)
c       call interij(np,contour,pp,BBin,pr,Bp)
c       call interij(np,contour,pp,AAin,pr,Ap)
c       call massf(np,pp,CCin,dsqrt(q2),Cq)

       sumav(1) =0d0
       sumav(2) =0d0
       sumav(3) =0d0
       sumav(4) =0d0
       sumav(5) =0d0
       sumav(6) =0d0
 
       i = 1
       flag1 = .false. 
   
30         if( (( abs(kx2(1)-k2).ge.pmin2).or.
     .          ( abs(kx2(2)-k2).ge.pmin2).or.
     .          ( abs(kx2(3)-k2).ge.pmin2)).and.
     .          (flag1.eqv..false.)   )then
 

       if(linear.eqv..false.)then
       z(1) = xp(i)*lambda**2+(1d0-xp(i))*pmin2 !+(0,1d0)*1d-10
       z(2) = xv(i)*eplus*mh-(1d0-xv(i))*eminus*mh ! +(0,1d0)*1d-10 
       z(3) = xp(i)*lambda**2+(1d0-xp(i))*pmin2 !+(0,1d0)*1d-10          
    
       dz(1)  =  (lambda**2-pmin2)
       dz(2)  =  (eplus+eminus)*mh
       dz(3)  = -(lambda**2-pmin2)

       dkx2(1)=(1d0-(0,1d0)*eminus*mh/dsqrt(z(1)))*dz(1)
       dkx2(2)=2d0*(-z(2)+(0,1d0)*lambda)*dz(2)
       dkx2(3)=(1d0+ (0,1d0)*eplus*mh/dsqrt(z(3)))*dz(3)

       kx2(1) = z(1)-(eminus*mh)**2-2d0*(0,1d0)*dsqrt(z(1))*eminus*mh
       kx2(2) = lambda**2-z(2)**2  +2d0*(0,1d0)*z(2)*lambda
       kx2(3) = z(3)-(eplus*mh)**2 +2d0*(0,1d0)*dsqrt(z(3))*eplus*mh

       else

       z(1) = xp(i)*lambda+(1d0-xp(i))*dsqrt(pmin2) !+(0,1d0)*1d-10
       z(2) = xv(i)*eplus*mh-(1d0-xv(i))*eminus*mh  ! +(0,1d0)*1d-10
       z(3) = xp(i)*lambda+(1d0-xp(i))*dsqrt(pmin2) !+(0,1d0)*1d-10

       dz(1)  =  (lambda-dsqrt(pmin2))
       dz(2)  =  (eplus+eminus)*mh
       dz(3)  = -(lambda-dsqrt(pmin2))

       dkx2(1)=2d0*(z(1)-(0,1d0)*eminus*mh)*dz(1)
       dkx2(2)=2d0*(-z(2)+(0,1d0)*lambda)*dz(2)
       dkx2(3)=2d0*(z(3)+ (0,1d0)*eplus*mh)*dz(3)

             
       kx2(1) = z(1)**2-(eminus*mh)**2-2d0*(0,1d0)*z(1)*eminus*mh
       kx2(2) = lambda**2-z(2)**2     +2d0*(0,1d0)*z(2)*lambda
       kx2(3) = z(3)**2-(eplus*mh)**2 +2d0*(0,1d0)*z(3)*eplus*mh

       endif

  
        sumav(1) =1d0*(
     .  wp(i)*BBin(i,1)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*BBin(i,2)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*BBin(i,3)*dkx2(3)/(kx2(3)-k2))+sumav(1)


       
        sumav(2) =1d0*(
     .  wp(i)*AAin(i,1)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*AAin(i,2)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*AAin(i,3)*dkx2(3)/(kx2(3)-k2))+sumav(2)

    
         sumav(3) =1d0*(
     .  wp(i)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*dkx2(3)/(kx2(3)-k2))+sumav(3)


        

         
       
        i = i+1
        li =i
        if(i.eq.np)flag1 = .true.
        goto 30
        elseif(flag1.eqv..true.)then
        Bk = sumav(1)/sumav(3)
        Ak = sumav(2)/sumav(3)
        else
        Bk = 0d0
        Ak = 0d0
        endif


      if((abs(Ak).eq.0).and.(abs(Bk).eq.0))goto 40

      if(flag2.eqv..true.)then

        


      write(3,*)dble(k2),dimag(k2),dble(Bk)
      write(4,*)dble(k2),dimag(k2),dimag(Bk)
      write(8,*)dble(k2),dimag(k2),dble(Ak)
      write(9,*)dble(k2),dimag(k2),dimag(Ak)
      
      pr2v(ia,ib,1)=dble(k2)
      pr2v(ia,ib,2)=dimag(k2)

      br(ia,ib)=dble(Bk)
      bi(ia,ib)=dimag(Bk)
      ar(ia,ib)=dble(Ak)
      ai(ia,ib)=dimag(Ak)
      write(2,*)pr2v(ia,ib,2) 

       else
       endif

       Bk2 = Bk*Bk
       Ak2 = Ak*Ak
c       Bp2 = Bp*Bp
c       Ap2 = Ap*Ap
      
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
        llambda(1) = Cq*Fgh(q2)*(Ak+Ap)/2d0
        llambda(2) = Cq*Fgh(q2)*(Ak-Ap)/(k2-p2+1d-10)/2d0  ! lambdaE
        llambda(3) =  -Cq*Fgh(q2)*(Bk-Bp)/(k2-p2+1d-10)
        else
        endif


       Lb(1) = (d-1)*Bk
       Lb(2) = -4d0*Bk*delta2/q2
       Lb(3) = -2d0*Ak*delta2/q2
       Tb(2) = -2d0*Bk*delta2
       Tb(3) = -(d-1)*Bk*q2
       Tb(6) =  (d-1)*Bk*(k2-p2)


        
c        La(1) =(2d0-d)*q2/2d0+(d-3d0)*(k2+p2)/2+(k2-p2)*(k2-p2)/(2*q2)   !RW   
        La(1) = 3*kp+2d0*(delta2)/q2  ! pk
        La(2) = 2*(k2+p2)*Delta2/q2   ! RW I 
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
                  



c---------------------------------------------------------------------------------------------------(ikernel 1)

       if(ikernel.eq.1)then 
       
       suma = 0d0

       do 100 i = 1,3
       suma = Lb(i)*llambda(i)+suma
 100   continue


       tr = cf*mathcalg/q2*1d0/(Ak2*k2+Bk2)*suma  ! mathcalg = 4*pi*alpha = g^2

c---------------------------------------------------------------------------------------- (ikernel 2)
       elseif(ikernel.eq.2)then

     
       suma = 0d0

       do 104 i = 1,3
       suma = La(i)*llambda(i)+suma
 104   continue

       
       tr = (cf/p2)*mathcalg/q2*Ak/(Ak2*k2+Bk2)*suma  ! mathcalg = 4*pi*alpha
c---------------------------------------------------------------------------------------------(ikernel 3)
       elseif(ikernel.eq.3)then
        
       CA = 3d0 

c       tr = CA*g2(k2,z2g)/(8d0*k2)*(-delta2)/k2*Fgh(v2)/v2*Fgh(k2)    !**2
c     . *Au*(Au+Ap)/(Au2*u2+Bu2)*Ck  !*h1(u2)                                              
c      tr = CA*g2(k2,z2g)/(8d0*k2)*(-delta2)/k2*Fgh(v2)/v2*Fgh(k2)
c     . *2d0/(u2)          
c      tr = CA*g2(k2,z2g)/(4d0*k2)*(-delta2)/k2*Fgh(u2)*Fgh(k2) !**2
c     . *1d0/(u2)**2

        tr =1d0

       else
       endif

c       write(*,*)"vertex inside"
       
  
      
       if(print.eqv..true.)then  
       write(2,*)'tr,suma,mathcalg,k2,Bk,Ak',tr,suma, 
     . mathcalg,k2,Bk,Ak
       write(2,*)'mathcalg',mathcalg
       write(2,*)'q2',q2
       write(2,*)'Ak',Ak
       write(2,*)'Bk',Bk
       write(2,*)'k2',k2
       write(2,*)'z1',z1
       write(2,*)'p2',p2
       write(2,*)'pmin2',pmin2
       write(2,*)'AAinm'
       write(2,*)(AAin(i,1),i=1,np)
       write(2,*)'AAinv'
       write(2,*)(AAin(i,2),i=1,np)
       write(2,*)'AAinp'
       write(2,*)(AAin(i,3),i=1,np)
       write(2,*)'bbinm'
       write(2,*)(bbin(i,1),i=1,np)
       write(2,*)'bbinv'
       write(2,*)(bbin(i,2),i=1,np)
       write(2,*)'bbinp'
       write(2,*)(bbin(i,3),i=1,np)
       write(2,*)'ccin'
       write(2,*)(ccin(i),i=1,np)
       write(2,*)'li',li
       write(2,*)'z(i)'
       write(2,*)(z(i),i=1,3)
       write(2,*)'kx2(1),z(1)',kx2(1),z(1)
       write(2,*)'mathcalg,q2,(Ak2*k2+Bk2)',
     . mathcalg,q2,(Ak2*k2+Bk2)
       write(2,*)'La(1),llambda(1)',La(1),llambda(1)
       write(2,*)(llambda(i),i=1,3)
       write(2,*)'mh',mh
       write(2,*)'pc4',pc4   
       write(2,*)'kx2'
       write(2,*)(kx2(i),i=1,3)
       write(2,*) 'suma', suma
       write(2,*)'sumav(1)',sumav(1)
       write(2,*)'sumav(2)',sumav(2)
       print = .false.
       endif

40     continue
       if((abs(Ak).eq.0).and.(abs(Bk).eq.0))tr=0d0
 

       return
       end




       


      

