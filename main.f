      program main
      implicit none
      include 'common.f'
       double precision tol1,Bmu,Amu,mu,ppr(128),bbr(128),AAr(128),
     . aaf1r(128,3),aaf1i(128,3),bbf1r(128,3),bbf1i(128,3),
     . aaf2r(128,3),aaf2i(128,3),bbf2r(128,3),bbf2i(128,3),
     . ppf1r(128,3),ppf1i(128,3),ppf2r(128,3),ppf2i(128,3),
     . mcurr(0:10),ppc(250),cc(250),ccin(128),alatt,blatt,
     . mhv(0:10),rtsec,z2r,z2rh,norm,func,fpi,rh,xpa(128),wpa(128),
     . xpb(128),wpb(128),xpc(128),wpc(128),xsa(128),wsa(128),
     . xsb(128),wsb(128),err,dlamb,z4r,z4rh,gammam,mcurr19(2),
     . mhat1,mhat2,minfty,pinfty,ppr2(128),bbrh(128),aarh(128)


       double complex z2,z4,z2h,z4h,norm22
                
      integer i,j,im,im0,nres,counter,nm_,units(19)
      logical finfty
 
      external kernint
      external massf
      external ktr
      external eigenv
      external sde
      external sder

      open (1,file="salida1.out",status="new")
      open (2,file="salida2.out",status="new")
      ! real axis sde files
      open (3,file="salida3.out",status="new")
      open (4,file="salida4.out",status="new")
      open (8,file="salida8.out",status="new")
      open (9,file="salida9.out",status="new")
      ! end  real axis sde

      open (10,file=
     ."datafiles/bbr-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
     . ,status="old")
      open (11,file=
     ."datafiles/bbi-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
     . ,status="old")
      open (12,file=
     ."datafiles/aar-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
     . ,status="old")
      open (13,file=
     ."datafiles/aai-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
     . ,status="old")
 


c      open (31,file=
c     ."datafiles/bbr-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
c     ., status="unknown")
c      open (32,file=
c     ."datafiles/bbi-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
c     . ,status="unknown")
c      open (33,file=
c     ."datafiles/aar-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
c     . ,status="unknown")
c      open (34,file=
c     ."datafiles/aai-u3p4-mu2-mh0p3-cd2011-L10-ddd.dat"
c     . ,status="unknown")

      open (35,file="datafiles/F.dat",
     . status="old")
      open (36,file="datafiles/B-Dw0p8-w0p6-mu2-mq4p0-pinfty10to3.dat",
     . status="old")
      open (37,file="datafiles/A-Dw0p8-w0p6-mu2-mq4p0-pinfty10to3.dat",
     . status="old")
       open (38,file="FT.out",status="new")
      open (39,file="datafiles/B-Dw0p8-w0p6-mu2-mq82-pinfty10to3.dat",
     . status="old")
      open (40,file="datafiles/A-Dw0p8-w0p6-mu2-mq82-pinfty10to3.dat",
     . status="old")

      ! complex sde files
      open (14,file="salida14.out",status="new")
      open (15,file="salida15.out",status="new")
      open (16,file="salida16.out",status="new")
      open (17,file="salida17.out",status="new")
      open (18,file="salida18.out",status="new")
      open (19,file="salida19.out",status="new")
      open (20,file="datafiles/x0_lr.dat",status="old")
      !  end sde 
      open (22,file="kbsr.out",status="new")
      open (23,file="kbsi.out",status="new")
      open (24,file="datafiles/kbsr.dat",status="old")
      open (25,file="datafiles/kbsi.dat",status="old")
      open (26,file="kijr.out",status="unknown")
      open (27,file="kiji.out",status="unknown")
      open (28,file="kijr2.out",status="unknown")
      open (29,file="kiji2.out",status="unknown")
      open (30,file="F.out",status="unknown")

      np = 128 ! for fermions internal calculation
      np2 =128 ! "
      ns= 64   ! for kernel      
      nt=  20    !16 ! change in common.f ! kernel
      nzq= 20  ! "
      nres = 1 !  resonancia 1,2,3
      nm_ = 0  ! number of chevychev moments  for kernel go to 20 to undertand
      nm  = nm_
 
      lamb_bs=2900d0                !dsqrt(40000d0)  ! UV cutuff for the bethe salpeter
      z2r =0.78201468395236029d0 !                !0.78201426196155177d0     !.782481338d0      !.971340499d0
      z4r =0.38571706705012776d0 !                      !0.37335477729364813d0     !.384259014        !.638852374d0
      z2rh =  0.83495261586028746d0
      z4rh =  0.41182927243446449d0
      do 114 i=1,128
      read(36,*)ppr(i),bbr(i)
      read(37,*)ppr(i),aar(i)
      read(39,*)ppr2(i),bbrh(i)
      read(40,*)ppr2(i),aarh(i)
114   continue
   
      do 118 j = 1,3
      do 119 i = 1,np

      read(10,200)ppf1r(i,j),ppf1i(i,j),BBf1r(i,j)
      read(11,200)ppf1r(i,j),ppf1i(i,j),BBf1i(i,j)
      read(12,200)ppf1r(i,j),ppf1i(i,j),AAf1r(i,j)
      read(13,200)ppf1r(i,j),ppf1i(i,j),AAf1i(i,j)



c      read(31,200)ppf2r(i,j),ppf2i(i,j),BBf2r(i,j)
c      read(32,200)ppf2r(i,j),ppf2i(i,j),BBf2i(i,j)
c      read(33,200)ppf2r(i,j),ppf2i(i,j),AAf2r(i,j)
c      read(34,200)ppf2r(i,j),ppf2i(i,j),AAf2i(i,j)



119   continue
118   continue

      print     = .true.
      flag2     = .false.
      init      = .true.
      printkabt = .true.
      linear    = .true.
      ffit      = .true.

       
       mu = 2d0  ! pseudo renormalization scale
       eplus = 1d0/2d0
       eminus= 1d0/2d0
      call  quadrature(lamb_bs,xp,wp,np)
      call  gauleg(0d0,1d0,xv,wv,np2) 
      if(ffit.eqv..true.)then     
      call  gauleg(-1d0,1d0,xs,ws,ns)
      else
       call  quadrature(lamb_bs,xs,ws,ns)
      endif
      call  gauleg(-1d0,1d0,xq,wq,nzq)
      call  gauleg(-1d0,1d0,xt,wt,nt)
      
      tol1 = 0.001d0

      pmin2= 1d-14


      im0 = 0
c      do 116 im = im0,10

      
      mhv(1) = 1.4d0  !+im*0.1d0
      mhv(2) = 1.4d0  !2.6d0  !10d0   !4.2d0 
      

      mh = mhv(1)
      mh2 =mhv(2)

  

      mcurr19(1)= 3.4d0/1000d0 !    3.4d0/1000d0  !3.145d0/1000d0      !5d0/1000d0 !  0.1d0/1000d0
      mcurr19(2)= 3.4d0/1000d0  !    3.4d0/1000d0  !3.145d0/1000d0      !3.71d0           !           3.4d0/1000d0     !0.2d0
     

      data units/1,2,3,4,8,9,14,15,16,17,18,19,22,23,26,
     .           27,28,29,30/         

      do  100 i= 1,19
      write(units(i),*)"# mcurr(1),mh(1)"   ,mcurr(1),mh
100   continue

          

         gammam = 12d0/(33d0-2d0*4d0)
         mhat1 = mcurr19(1)*(1d0/2d0*Log((19d0/0.234d0)**2))**gammam 
         mcurr(1)= mhat1/(1d0/2d0*Log((2d0/0.234d0)**2))**gammam
         mhat2 = mcurr19(2)*(1d0/2d0*Log((19d0/0.234d0)**2))**gammam
         mcurr(2)= mhat2/(1d0/2d0*Log((2d0/0.234d0)**2))**gammam
         pinfty = 1000d0
         
 

        ! fermion 1
         finfty = .true.
         minfty  = mhat1/(1d0/2d0*Log((pinfty/0.234d0)**2))**gammam
         write(2,*)"minfty,mu,mhat1",minfty,mu
         Bmu = mcurr(1)
cc        sder(finfty,pinfty,minfty,mu,mmu,B(mu),A(mu),lambda,128,ppr,bbr,aar,z2r,z4r)
cc        if(im0.eq.im)then
c        call sder(finfty,pinfty,minfty,mu,mcurr(1),Bmu,1d0,3000d0,128,
c     .  ppr,bbr,aar,z2r,z4r)
c        call massf(128,ppr,bbr,mu,Bmu)
c        call massf(128,ppr,aar,mu,Amu)
        finfty  = .false. 
         mu =19d0
         mcurr(1)= 3.4d0/1000d0
         mcurr(2)= 3.4d0/1000d0
         Bmu =     3.4d0/1000d0
         Amu = 1d0                                        !lamb_bs
        call sder(finfty,pinfty,minfty,mu,mcurr(1),Bmu,1d0,lamb_bs,128,
     .  ppr,bbr,aar,z2r,z4r)

            
        npol = 3
        call sderfit(ppr,bbr,aar,xval,npol)
 
        z2 = z2r
        z4 = z4r

c        else  
c        endif

       
       write(2,*)'mu,Bmu'
       write(2,10)mu,Bmu
       write(2,*)'mu,Amu'
       write(2,10)mu,Amu
c       call sde(mh,mu,Bmu,Amu,mcurr(1),dsqrt(lamb_bs**2+mh**2/4d0),
c     . np2,ccin,z2,z4,bbf1r,bbf1i,aaf1r,aaf1i)
      write(2,*)"z2",z2
      write(2,*)"z4",z4 
      write(*,*)"complex sde ready"
      write(2,*)"complex sde ready"
      do 120 j = 1,3
      do 121 i = 1,np
      BB(i,j)= BBf1r(i,j)+(0d0,1d0)*BBf1i(i,j)
      aa(i,j)= aaf1r(i,j)+(0d0,1d0)*aaf1i(i,j)
c       BB(i,j)= BBf2r(i,j)+(0d0,1d0)*BBf2i(i,j)
c       aa(i,j)= aaf2r(i,j)+(0d0,1d0)*aaf2i(i,j)

121   continue
120   continue


       ! fermion 2
          finfty  = .true.
          minfty  = mhat2/(1d0/2d0*Log((pinfty/0.234d0)**2))**gammam
         Bmu = mcurr(2)
c        sder(mu,B(mu),A(mu),lambda,128,ppr,bbr,aar,z2rh)
c         call sder(finfty,pinfty,minfty,mu,mcurr(2),Bmu,1d0,3000d0,128,
c     .   ppr2,bbrh,aarh,z2rh,z4rh)

       call massf(128,ppr2,bbrh,mu,Bmu)
       call massf(128,ppr2,aarh,mu,Amu)
       write(2,*)'mu,Bmu',mu,Bmu
       write(2,*)'mu,Amu',mu,Amu
c          finfty  = .false. 
c       call sder(finfty,pinfty,minfty,mu,mcurr(2),Bmu,1d0,lamb_bs,128,
c     .  ppr,bbr,aar,z2rh,z4rh)
   
           npol2 = 3
c          call sderfit(ppr2,bbrh,aarh,xval2,npol2)
         write(2,*)"xval2"
         do 129 i = 1,26
         xval2(i)= xval(i)
         write(2,*)i,xval2(i)
129      continue

        z2h = z2rh
        z4h = z4rh



c      call sde(mh2,mu,Bmu,Amu,mcurr(2),dsqrt(lamb_bs**2+mh2**2/4d0),
c     . np2,ccin,z2h,z4h,bbf2r,bbf2i,aaf2r,aaf2i) 
      write(*,*)"complex sde ready"
      write(2,*)"complex sde ready"
      do 122 j = 1,3
      do 123 i = 1,np
c      BB2(i,j)= BBf2r(i,j)+(0d0,1d0)*BBf2i(i,j)  
c      aa2(i,j)= aaf2r(i,j)+(0d0,1d0)*aaf2i(i,j)
      BB2(i,j)= BBf1r(i,j)+(0d0,1d0)*BBf1i(i,j)  ! comment if f1 diff of f2
      aa2(i,j)= aaf1r(i,j)+(0d0,1d0)*aaf1i(i,j)
123   continue
122   continue




       lambda = lamb_bs
       write(2,*)"np,np2,nzq,ns,nt,nm",np,np2,nzq,ns,nt,nm
       write(2,*)"mh,mu,Bmu,Amu,lambda,np2",mh,mu,Bmu,Amu,
     . lambda,np2   

116    continue


          z2rh = z2r
          z2k = dsqrt(z2r*z2rh)
          write(2,*)"z2k",z2k
         
      
         rtsec =0.770d0      !0.13178664167051821d0              !.133157632d0           !0.128404692d0           !.133157632d0       
        call sec(rtsec*(1d0-0.0001d0),rtsec*(1d0+0.0005),
     .  tol1,nres,rtsec)
        write(2,*)"improved pion mass",rtsec 
c        rtsec = 0.12618063d0      
c        call  eigenv(rtsec,nres,func)
c        write(2,*)"eigenv ready"  
cc        call derivative(dlamb,rtsec,0.005d0,err,nres)
cc        write(2,*)"dlamb,rtsec,err",dlamb,rtsec,err

          
         
        call   decayc(ffit,rtsec,lambda,mh,mh2,aa,bb,
     . aa2,bb2,xval,npol,xval2,npol2,z2r,z2rh,norm,fpi,rh,norm22)
        write(2,*)"fpi",fpi
       write(2,*)"fpi-from-rh/z4",rh*(mcurr(1)+mcurr(2))/rtsec**2
       write(2,*)"norm^2 dlamb/ds",norm22*dlamb
       write(2,*)"real part of the norm dlamb/ds",
     .  dsqrt(abs(dble(norm22*dlamb)))   
        write(2,*)"real part of the norm ",
     .  dsqrt(abs(dble(norm22)))

       write(2,*)"z2,z4,z2h,z4h",z2,z4,z2h,z4h


c        call derivative(dlamb,rtsec,0.005d0,err,nres)
c        write(2,*)"dlamb,rtsec,err",dlamb,rtsec,err
        


           

c116    continue


 
      
10     format(D25.16,D25.16)
200    format(D25.16,D25.16,D25.16)


      end

