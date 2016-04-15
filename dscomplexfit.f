       subroutine dscomplexfit(print,mn,pc4,q2E,z1,xval,npol,sigv,sigs,
     . dsigv,dsigs)
      implicit none
      double precision zr(3),zi(3),mr(3),mi(3),xval(26),q2E,z1 
      double precision krv(4),kr,pmin2,mn,ai(0:6),bi(0:6),lamb,lamb2
      double complex pc4,pIv(4),k2,dk2,Bk,Ak,dBk,dAk,bk2,ak2
      double complex sigs,sigv,zk(3),zkc(3),mk(3),mkc(3),dsigv,dsigs,D
      integer k,npol,i
      logical print  
 
      pmin2 = 1d-14


             do 102 k = 1,3

            zr(k) = xval(k)
            zi(k) = xval(k+3)
            mr(k) = xval(k+6)
            mi(k) = xval(k+9)
102         continue
 
             do 103  k = 0,6
            bi(k)  = xval(k+13)
            ai(k)  = xval(k+20)
103         continue
            
        


       kr = abs(dsqrt(q2E))

       krv(3) =kr*dsqrt(1d0-z1**2)
       krv(4) =kr*z1

       pIv(3) =0d0
       pIv(4) =pc4
       
       k2 = (krv(3)+pIv(3))*(krv(3)+pIv(3))
     .     +(krv(4)+pIv(4))*(krv(4)+pIv(4))

           




            zk(1)=  zr(1)+zi(1)*(0d0,1d0)
            zk(2)=  zr(2)+zi(2)*(0d0,1d0)
            zk(3)=  zr(3)+zi(3)*(0d0,1d0)

            zkc(1)= zr(1)-zi(1)*(0d0,1d0)
            zkc(2)= zr(2)-zi(2)*(0d0,1d0)
            zkc(3)= zr(3)-zi(3)*(0d0,1d0)




            mk(1) = mr(1)+mi(1)*(0d0,1d0)
            mk(2) = mr(2)+mi(2)*(0d0,1d0)
            mk(3) = mr(3)+mi(3)*(0d0,1d0)


            mkc(1)= mr(1)-mi(1)*(0d0,1d0)
            mkc(2)= mr(2)-mi(2)*(0d0,1d0)
            mkc(3)= mr(3)-mi(3)*(0d0,1d0)





        sigs = 0d0
        sigv = 0d0
        dsigs = 0d0
        dsigv = 0d0
        Ak = 0d0
        Bk = 0d0
        dAk = 0d0
        dBk = 0d0

        lamb = 10d0
        lamb2 = lamb*lamb

        do 101 k = 1,npol
       if(dble(k2).le.lamb2)then
        sigs   =  zk(k)*mk(k)/(k2+mk(k)**2)
     .           +zkc(k)*mkc(k)/(k2+mkc(k)**2)+sigs
        sigv  =  zk(k)/(k2+mk(k)**2)
     .         + zkc(k)/(k2+mkc(k)**2)+sigv

        dsigs   = -zk(k)*mk(k)/(k2+mk(k)**2)**2
     .            -zkc(k)*mkc(k)/(k2+mkc(k)**2)**2+dsigs
        dsigv  =  -zk(k)/(k2+mk(k)**2)**2
     .            -zkc(k)/(k2+mkc(k)**2)**2+dsigv
        else
        endif
101      continue
       
        do 104 k = 1,6
        if(dble(k2).gt.lamb2)then
         Bk = (abs(Bi(k))/k2)**k*Bi(k)/abs(bi(k)+1d-10)+Bk
        Ak = (abs(ai(k))/k2)**k*ai(k)/abs(ai(k)+1d-10)+ak
        dBk = -k*(abs(Bi(k))/k2)**(k+1)*Bi(k)/abs(bi(k)+1d-10)**2+dBk  ! dBk/dk2
        dAk = -k*(abs(ai(k))/k2)**(k+1)*ai(k)/abs(ai(k)+1d-10)**2+dak
        else
        endif
104     continue

        if(dble(k2).gt.lamb2)then
        bk = bi(0)+bk
        ak = ai(0)+ak
        Bk2 = Bk*BK
        Ak2 = Ak*Ak
        sigs = Bk/(Ak*Ak*k2+Bk*BK) 
        sigv = Ak/(Ak*Ak*k2+Bk*BK)
        else
        endif

         dk2=    2d0*(krv(4)+pIv(4))*pIv(4)/(mn)  ! dk^2/dk^\mu*k^\mu

          
 

         if(dble(k2).le.lamb2)then
         dsigs  = dsigs*dk2
         dsigv  = dsigv*dk2
         Bk =   sigs/(sigs**2 + k2*sigv**2)
         Ak =   sigv/(sigs**2 + k2*sigv**2)
         else
       D= (k2*Ak2+Bk2+pmin2)
       dsigs= (dBk*D-Bk*(2d0*k2*Ak*dAk+2d0*Bk*dBK+Ak2))*dk2/D**2
       dsigv= (dak*D-Ak*(2d0*k2*Ak*dAk+2d0*Bk*dBk+Ak2))*dk2/D**2
       endif


         if( !c     . (p2.eq.pmin2).and.
     . print.eqv..true.)then
       write(2,*)"dscomplexfit"
       write(2,*)'q2E',q2E
       write(2,*)'Ak',Ak
       write(2,*)'Bk',Bk
       write(2,*)'k2',k2
       write(2,*)'z1',z1
       write(2,*)'pmin2',pmin2
       write(2,*)"D",D
       write(2,*)"zr(3),zi(3),mr(3),mi(3)"
       write(2,*)(zr(i),i=1,3)
       write(2,*)(zi(i),i=1,3)
       write(2,*)(mr(i),i=1,3)
       write(2,*)(mi(i),i=1,3)
       write(2,*)"xval"
       write(2,*)(xval(i),i=1,12)
       write(2,*)"sigs,sigv,dsigv,dsigs"
       write(2,*) sigs,sigv,dsigv,dsigs
       write(2,*)"zk(3),zkc(3),mk(3),mkc(3)"
       write(2,*)(zk(i),i=1,3)
       write(2,*)(zkc(i),i=1,3)
       write(2,*)(mk(i),i=1,3)
       write(2,*)(mkc(i),i=1,3)
       write(2,*)'(Ak2*k2+Bk2)',(Ak2*k2+Bk2)
       write(2,*)'pc4',pc4
       write(2,*)"k,npol",k,npol
c       print = .false.
       endif


      
      end
