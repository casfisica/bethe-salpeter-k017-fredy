      subroutine fcn(npar,grad,fval,xval,iflag,chi2)

      implicit none
       include 'commonfit.f'
      integer npar,iflag,i,j
c     double precision grad(npar),fval,xval(npar),chi2
      double precision grad(26),fval,xval(26),chi2,pull(250),
     . pull2(250)
      double complex smval(250),smval2(250) 
      external chi2
      npar =26
c      fval = chi2(xval,npar,smval,pull)
       fval = chi2(xval,npar,smval,smval2,pull,pull2)

       if(iflag.eq.3)then
       write(7,*)'Bpulls'
       write(7,*)"smval(i),pull(i)"
       write(7,104)(dble(smval(i)), dble(pp128(i)),pull(i),i=1,128)
       write(7,*)'Apulls'
       write(7,104)(dble(smval2(i)),dble(pp128(i)),pull2(i),i=1,128)
       write(7,*)'Bpulls'
       write(7,*)"smval(i),pull(i)"
       write(7,104)(dble(smval(i)),dble(B(i,2)),
     . abs((smval(i)-B(i,2))/(B(i,2)*0.01d0)),i=1,128)
       write(7,*)'Apulls'
       write(7,104)(dble(smval2(i)),dble(A(i,2)),
     . abs((smval2(i)-A(i,2))/(A(i,2)*0.01d0)),i=1,128)

       do 100 i = 1,26
       xval2(i)= xval(i)
100    continue

       else
       endif
c104    format(D18.9,D18.9) 
104       format(F18.9,F18.9,F18.9)  
      return
      end

      double precision function chi2(xval,npar,smval,smval2,pull,pull2)
      
      implicit none
      include 'commonfit.f'
      integer i,j,k,npar
      double precision xval(npar),ai(0:6),bi(0:6),lamb
      double precision p2,zr(6),zi(6),mr(6),mi(6),pull(250),pull2(250)
      double complex zk(6),zkc(6),mk(6),mkc(6),sigv,sigs,smval(250),
     . value(250),error(250),smval2(250),Ak,Bk
c      data value/  
c     .     0.3072d0, 0.3093d0, 0.2917d0, 0.0298d0, 0.2820d0, 0.0440d0/
c      
c      data error/
c     .     0.0033d0, 0.0031d0, 0.0093d0, 0.0087d0, 0.0142d0, 0.0136d0/

            do 102 k = 1,npol
            zr(k) = xval(k)
            zi(k) = xval(k+npol)
            mr(k) = xval(k+2*npol)
            mi(k) = xval(k+3*npol)
102         continue

 
             do 103 k = 0,6 
            bi(k)  = xval(k+13)
            ai(k)  = xval(k+20)
103         continue
            

            do 105 k = 1,npol
            zk(k)=  zr(k)+zi(k)*(0d0,1d0)
            zkc(k)= zr(k)-zi(k)*(0d0,1d0)
            mk(k) = mr(k)+mi(k)*(0d0,1d0)
            mkc(k)= mr(k)-mi(k)*(0d0,1d0)
105         continue
      
 
           

            

        chi2 = 0d0
 
        do 100 i = 1,128


        p2 = pp128(i)**2
        
        sigs = 0d0
        sigv = 0d0
        Ak = 0d0
        Bk = 0d0
        lamb = 10d0

        if(pp128(i).le.lamb)then
        do 101 k = 1,npol
        sigs   =  zk(k)*mk(k)/(p2+mk(k)**2)
     .           +zkc(k)*mkc(k)/(p2+mkc(k)**2)+sigs
        sigv  =  zk(k)/(p2+mk(k)**2)
     .         + zkc(k)/(p2+mkc(k)**2)+sigv
101     continue
        else
        endif
       

        if(pp128(i).gt.lamb)then
        bk=0d0
        ak=0d0
        do 104 k = 1,6
        Bk = (abs(Bi(k))/p2)**k*Bi(k)/abs(bi(k)+1d-10)+Bk
        Ak = (abs(ai(k))/p2)**k*ai(k)/abs(ai(k)+1d-10)+ak
104     continue
        else
        endif

        if(pp128(i).le.lamb)then
        smval(i)= sigs/(sigs**2 + pp128(i)**2*sigv**2)
        smval2(i)=sigv/(sigs**2 + pp128(i)**2*sigv**2)
        else
        endif

        if(pp128(i).gt.lamb)then
        smval(i) = Bk+bi(0)             !dlog(pp128(i)**2)
        smval2(i)= Ak+ai(0)
        else
        endif

                  
        pull(i) = abs( (smval(i)-B(i,2))/(0.01d0*B(i,2)))
        pull2(i)= abs((smval2(i)-A(i,2))/(0.01d0*A(i,2)))
c        pull(i) =  abs((smval(i)-B(i,2))/0.01d0)
c        pull2(i)= abs((smval2(i)-A(i,2))/0.01d0)

          
c        if(pp128(i).lt.pp128(128))then  ! 104
        chi2 = chi2+pull(i)**2+pull2(i)**2
c        else
c       endif
100     continue

        
      return
      end
