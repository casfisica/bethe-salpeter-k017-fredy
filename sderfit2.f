      subroutine  sderfit(ppr,BB,AA,xval2_,npol2)

      implicit none
      include 'commonfit.f'
      double precision chi2,AA(128),bb(128),
     . ppr(128),xval2_(26)
      integer i,npol2
      external fcn
      external chi2
c      open (35,file="B-Dw0p87-w0p5-mu2-mhat6p4.dat",status="old")
c      open (36,file="A-Dw0p87-w0p5-mu2-mhat6p4.dat",status="old")
     
      
      do 1 i = 1,128
c      read(1,*)B(i,1),B(i,2)
c      read(2,*)A(i,1),A(i,2)
       B(i,2)= BB(i)
       A(i,2)= AA(i)
      pp128(i)= ppr(i)

      sigmas(i)=B(i,2)/(pp128(i)**2*A(i,2)**2+B(i,2)**2) 
      sigmav(i)=A(i,2)/(pp128(i)**2*A(i,2)**2+B(i,2)**2)
                
      
1     continue

       
         
        npol= npol2


      call mintio(5,7,9)

      open (5,file='smfit_2.dat',status='old')
      open (6,file='/dev/null',status='unknown')
      open (7,file='smfit_2.out',status='unknown')

      call minuit(fcn,chi2)

       do 100 i = 1,26
       xval2_(i)= xval2(i)
100    continue


c      stop
      return
      end
