       subroutine dscomplex(print,mn,pc4,k2r,z1,sigmav,sigmas,
     .                       dsigmav,dsigmas)
       implicit none
       include 'commonf.f'
       double precision mathcalg,kpr,p2,
     . g2,z(3),dz(3),mn
       double precision q2,z1,kr,k2r,p2r,krv(4)
       double complex pIv(4),k2,ak2,bk2,ak,bk,
     . sumav(6),tr,kv(4),dk2,dsigmav,dsigmas,D,
     .  kx2(3),dkx2(3),suma,sigmav,sigmas,pc4,dbk,dak
       integer  i,li

       logical flag1,print       

      
      

      
c       p2 = p*p

       kr = abs(dsqrt(k2r))
       kpr = dsqrt(p2r*k2r)*z1

       krv(3) =kr*dsqrt(1d0-z1**2)
       krv(4) =kr*z1

       pIv(3) =0d0
       pIv(4) =pc4

       kv(1)= 0d0 
       kv(2)= 0d0
       kv(3)= krv(3)+pIv(3)
       kv(4)= krv(4)+pIv(4) 
      
       k2 = (krv(3)+pIv(3))*(krv(3)+pIv(3))
     .     +(krv(4)+pIv(4))*(krv(4)+pIv(4))


c      dk2= (0d0,1d0)*kv(4)*pIv(4)/((0,1d0)*abs(pIv(4))) ! i_ k.uh
       dk2=       2d0*kv(4)*pIv(4)/(mn)

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


  
        sumav(1) =1d0*(
     .  wp(i)*BB(i,1)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*BB(i,2)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*BB(i,3)*dkx2(3)/(kx2(3)-k2))+sumav(1)


       
        sumav(2) =1d0*(
     .  wp(i)*AA(i,1)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*AA(i,2)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*AA(i,3)*dkx2(3)/(kx2(3)-k2))+sumav(2)

    
         sumav(3) =1d0*(
     .  wp(i)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*dkx2(3)/(kx2(3)-k2))+sumav(3)


        sumav(4) =1d0*(
     .  wp(i)*BB(i,1)*dkx2(1)/(kx2(1)-k2)**2+
     .  wv(i)*BB(i,2)*dkx2(2)/(kx2(2)-k2)**2+
     .  wp(i)*BB(i,3)*dkx2(3)/(kx2(3)-k2)**2)+sumav(4)
 

        sumav(5) =1d0*(
     .  wp(i)*AA(i,1)*dkx2(1)/(kx2(1)-k2)**2+
     .  wv(i)*AA(i,2)*dkx2(2)/(kx2(2)-k2)**2+
     .  wp(i)*AA(i,3)*dkx2(3)/(kx2(3)-k2)**2)+sumav(5)



         sumav(6) =1d0*(
     .  wp(i)*dkx2(1)/(kx2(1)-k2)**2+
     .  wv(i)*dkx2(2)/(kx2(2)-k2)**2+
     .  wp(i)*dkx2(3)/(kx2(3)-k2)**2)+sumav(6)

         

       
        i = i+1
        li =i
        if(i.eq.np)flag1 = .true.
        goto 30
        elseif(flag1.eqv..true.)then
        Bk = sumav(1)/sumav(3)
        Ak = sumav(2)/sumav(3)
        dBk =(sumav(4)-Bk*sumav(6))/sumav(3)
        dAk =(sumav(5)-Ak*sumav(6))/sumav(3)
        else
        Bk = 0d0
        Ak = 0d0
        dbk =0d0
        dak = 0d0
        endif


      if((abs(Ak).eq.0).and.(abs(Bk).eq.0))goto 40

c      if(flag2.eqv..true.)then

        
c      write(14,*)dble(k2),dimag(k2),dble(Bk/Ak)
c      write(15,*)dble(k2),dimag(k2),dble(1d0/Ak)
c      write(16,*)dble(k2),dimag(k2),dble(Bk)
c      write(17,*)dble(k2),dimag(k2),dble(Ak)
c      write(18,*)dble(k2),dimag(k2),dimag(Bk/Ak)
c      write(19,*)dble(k2),dimag(k2),dimag(1d0/Ak)
c      write(20,*)dble(k2),dimag(k2),dimag(Bk)
c      write(21,*)dble(k2),dimag(k2),dimag(Ak)

c       else
c       endif

       Bk2 = Bk*Bk
       Ak2 = Ak*Ak

       D= (k2*Ak2+Bk2+pmin2) 
       sigmas= Bk/D
       sigmav= Ak/D
      

       dsigmas= (dBk*D-Bk*(2d0*k2*Ak*dAk+2d0*Bk*dBK+Ak2))*dk2/D**2
       dsigmav= (dak*D-Ak*(2d0*k2*Ak*dAk+2d0*Bk*dBk+Ak2))*dk2/D**2


c        dBk =(sumav(4)-Bk*sumav(6))/sumav(3) 
c        dAk =(sumav(5)-Ak*sumav(6))/sumav(3) 

c       Bp2 = Bp*Bp
c       Ap2 = Ap*Ap
      
      ! I remove this part of vertex


  
 
             
       if( !c     . (p2.eq.pmin2).and.
     . print.eqv..true.)then  
       write(2,*)"scomplex"
       write(2,*)'tr,suma,k2,Bk,Ak',tr,suma, 
     . k2,Bk,Ak
        write(2,*)'g2',g2
       write(2,*)'q2',q2
       write(2,*)'Ak',Ak
       write(2,*)'Bk',Bk
       write(2,*)'k2',k2
       write(2,*)'z1',z1
       write(2,*)'p2',p2
       write(2,*)'pmin2',pmin2
       write(2,*)'AAinm'
       write(2,*)(AA(i,1),i=1,np)
       write(2,*)'AAinv'
       write(2,*)(AA(i,2),i=1,np)
       write(2,*)'AAinp'
       write(2,*)(AA(i,3),i=1,np)
       write(2,*)'bbinm'
       write(2,*)(bb(i,1),i=1,np)
       write(2,*)'bbinv'
       write(2,*)(bb(i,2),i=1,np)
       write(2,*)'bbinp'
       write(2,*)(bb(i,3),i=1,np)
       write(2,*)'li',li
       write(2,*)'z(i)'
       write(2,*)(z(i),i=1,3)
       write(2,*)'kx2(1),z(1)',kx2(1),z(1)
       write(2,*)'mathcalg,q2,(Ak2*k2+Bk2)',
     . mathcalg,q2,(Ak2*k2+Bk2)
c       write(2,*)'La(1),llambda(1)',La(1),llambda(1)
c       write(2,*)(llambda(i),i=1,3)
       write(2,*)'mh',mh
       write(2,*)'pc4',pc4   
       write(2,*)'kx2'
       write(2,*)(kx2(i),i=1,3)
       write(2,*) 'suma', suma
       write(2,*)'sumav(1)',sumav(1)
       write(2,*)'sumav(2)',sumav(2)
       write(2,*)'sumav(3)',sumav(3)
c       print = .false.
       endif

40     continue
       if((abs(Ak).eq.0).and.(abs(Bk).eq.0))then
       sigmas= 0d0
       sigmav= 0d0
       dAk=0d0
       dBk=0d0    
       else
       endif
 

       return
       end




       


      

