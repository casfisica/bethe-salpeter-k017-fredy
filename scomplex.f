       subroutine scomplex(pc4,k2r,z1,mhl,AAl,BBl,sigmav,sigmas)
       implicit none
       include 'common.f'
       double precision mhl,z(3),dz(3)
       double precision z1,kr,k2r,p2r,krv(4),prv(4)
       double complex pIv(4),k2,ak2,bk2,ak,bk,sumav(6),tr,
     .  kx2(3),dkx2(3),suma,sigmav,sigmas,pc4,BBl(128,3),AAl(128,3)
       integer  i,j,li

       logical flag1       

       external massf
       external inter

       lambda=dsqrt(lamb_bs**2+mhl**2/4d0)
      
       

       kr = abs(dsqrt(k2r))
       

       krv(3) =kr*dsqrt(1d0-z1**2)
       krv(4) =kr*z1

       pIv(3) =0d0
       pIv(4) =pc4

        
      
       k2 = (krv(3)+pIv(3))*(krv(3)+pIv(3))
     .     +(krv(4)+pIv(4))*(krv(4)+pIv(4))


       sumav(1) =0d0
       sumav(2) =0d0
       sumav(3) =0d0
       sumav(4) =0d0
       sumav(5) =0d0
       sumav(6) =0d0
 
       i = 1
       flag1 = .false. 
   
30     continue

c         if( (( abs(kx2(1)-k2).ge.pmin2).or.
c     .          ( abs(kx2(2)-k2).ge.pmin2).or.
c     .          ( abs(kx2(3)-k2).ge.pmin2)).and.
c     .          (flag1.eqv..false.)   )then
 

       

       z(1) = xp(i)*lambda+(1d0-xp(i))*dsqrt(pmin2) !+(0,1d0)*1d-10
       z(2) = xv(i)*eplus*mhl-(1d0-xv(i))*eminus*mhl  ! +(0,1d0)*1d-10
       z(3) = xp(i)*lambda+(1d0-xp(i))*dsqrt(pmin2) !+(0,1d0)*1d-10

       dz(1)  =  (lambda-dsqrt(pmin2))
       dz(2)  =  (eplus+eminus)*mhl
       dz(3)  = -(lambda-dsqrt(pmin2))

       dkx2(1)=2d0*(z(1)-(0,1d0)*eminus*mhl)*dz(1)
       dkx2(2)=2d0*(-z(2)+(0,1d0)*lambda)*dz(2)
       dkx2(3)=2d0*(z(3)+ (0,1d0)*eplus*mhl)*dz(3)

             
       kx2(1) = z(1)**2-(eminus*mhl)**2-2d0*(0,1d0)*z(1)*eminus*mhl
       kx2(2) = lambda**2-z(2)**2     +2d0*(0,1d0)*z(2)*lambda
       kx2(3) = z(3)**2-(eplus*mhl)**2 +2d0*(0,1d0)*z(3)*eplus*mhl

             if( (( abs(kx2(1)-k2).ge.pmin2).or.
     .          ( abs(kx2(2)-k2).ge.pmin2).or.
     .          ( abs(kx2(3)-k2).ge.pmin2)).and.
     .          (flag1.eqv..false.)   )then


  
        sumav(1) =1d0*(
     .  wp(i)*BBl(i,1)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*BBl(i,2)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*BBl(i,3)*dkx2(3)/(kx2(3)-k2))+sumav(1)


       
        sumav(2) =1d0*(
     .  wp(i)*AAl(i,1)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*AAl(i,2)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*AAl(i,3)*dkx2(3)/(kx2(3)-k2))+sumav(2)

    
         sumav(3) =1d0*(
     .  wp(i)*dkx2(1)/(kx2(1)-k2)+
     .  wv(i)*dkx2(2)/(kx2(2)-k2)+
     .  wp(i)*dkx2(3)/(kx2(3)-k2))+sumav(3)

         
             sumav(4) =1d0*(
     .  wp(i)*BBl(i,1)*dkx2(1)+
     .  wv(i)*BBl(i,2)*dkx2(2)+
     .  wp(i)*BBl(i,3)*dkx2(3))+sumav(4)


        sumav(5) =1d0*(
     .  wp(i)*AAl(i,1)*dkx2(1)+
     .  wv(i)*AAl(i,2)*dkx2(2)+
     .  wp(i)*AAl(i,3)*dkx2(3))+sumav(5)

        sumav(6) =1d0*(
     .  wp(i)*dkx2(1)+
     .  wv(i)*dkx2(2)+
     .  wp(i)*dkx2(3))+sumav(6)
 


         
       
        i = i+1
        li =i
        if(i.eq.np)flag1 = .true.
        goto 30
        elseif(flag1.eqv..true.)then
c        Kx2(1)=pmin2-(eminus*mhl)**2 
c        sumav(3)=sumav(3)-sumav(6)/(kx2(1)-k2)
c        Bk = (sumav(1)-sumav(6)*(BBL(1,1)+BBL(1,3))/2)/sumav(3)
c        Ak = (sumav(2)-sumav(6)*(AAL(1,1)+AAL(1,3))/2)/sumav(3)
        Bk = sumav(1)/sumav(3)
        Ak = sumav(2)/sumav(3)
 

        else
        Bk = 0d0
        Ak = 0d0
        endif


      if((abs(Ak).eq.0).and.(abs(Bk).eq.0))goto 40


       Bk2 = Bk*Bk
       Ak2 = Ak*Ak

       sigmas= Bk/(k2*Ak2+Bk2)
       sigmav= Ak/(k2*Ak2+Bk2)

     
      ! I remove this part of vertex


  
 
             
       if(print.eqv..true.)then  
       write(2,*)"scomplex"
       write(2,*)"linear", linear
       write(2,*)'tr,suma,k2,Bk,Ak',tr,suma, 
     . k2,Bk,Ak
       write(2,*)'Ak',Ak
       write(2,*)'Bk',Bk
       write(2,*)'k2',k2
       write(2,*)'z1',z1
       write(2,*)'pmin2',pmin2
       write(2,*)'AAinm'
       write(2,*)(AAl(i,1),i=1,np)
       write(2,*)'AAinv'
       write(2,*)(AAl(i,2),i=1,np)
       write(2,*)'AAinp'
       write(2,*)(AAl(i,3),i=1,np)
       write(2,*)'bbinm'
       write(2,*)(bbl(i,1),i=1,np)
       write(2,*)'bbinv'
       write(2,*)(bb(i,2),i=1,np)
       write(2,*)'bbinp'
       write(2,*)(bb(i,3),i=1,np)
       write(2,*)'li',li
       write(2,*)'z(i)'
       write(2,*)(z(i),i=1,3)
       write(2,*)'kx2(1),z(1)',kx2(1),z(1)
       write(2,*)'(Ak2*k2+Bk2)',
     . (Ak2*k2+Bk2)
       write(2,*)'mh',mhl
       write(2,*)'pc4',pc4   
       write(2,*)'kx2'
       write(2,*)(kx2(i),i=1,3)
       write(2,*) 'suma', suma
       write(2,*)'sumav(1)',sumav(1)
       write(2,*)'sumav(2)',sumav(2)
       write(2,*)'sumav(3)',sumav(3)
       write(2,*)"sumav(1)/sumav(3)",sumav(1)/sumav(3)
       write(2,*)"sumav(2)/sumav(3)",sumav(2)/sumav(3)
       write(2,*)'sumav(4)',sumav(4)
       write(2,*)'sumav(5)',sumav(5)
       write(2,*)'sumav(6)',sumav(6)
       print = .false.
       endif

40     continue
       if((abs(Ak).eq.0).and.(abs(Bk).eq.0))then
       sigmas= 0d0
       sigmav= 0d0
       else
       endif
  
       lambda=lamb_bs

       return
       end




       


      

