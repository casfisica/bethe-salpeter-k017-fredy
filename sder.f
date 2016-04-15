
       subroutine sder(finfty,pinfty,minfty,mu,mmu,Bmu,Amu,llambda,np3,
     . ppr,bb,aa,z2r,z4r)
       implicit none
      double precision sigmai,tol1,L,xmin,u,umin,umax,errormax,
     . error,A19,B19,z2r,z4r,mmu,mu,alatt,blatt,Bmu,Amu,pi,
     . minfty,pinfty

       double precision aa(np3),bb(np3),cc(128),p,ainfty

        double precision p2,lambda,pmin2,llambda,xp(128),wp(128)
        integer np,nvertex,ikernel,np3
        double precision ppr(np3),aain(128),bbin(128),ccin(128),
     .  pp(128)
        logical flag1,finfty



       

              
      integer i,j,maxit
      logical flagsalir
      external  intsder
      external massf
      external alatt
      external blatt
c      external inter
c      external SPLGR1
      external Fgh
      external g2
     

c      open (1,file="salida1.out",status="new")
c      open (2,file="salida2.out",status="new")
c      open (3,file="salida3.out",status="unknown")
c      open (4,file="salida4.out",status="unknown")
c      open (5,file="salida5.out",status="unknown")
c      open (7,file="salida6.out",status="new") ! avoid open (6,
c      open (8,file="X0_costura_1.0.dat",status="old")
     
  
       write(2,*)"sde for the real x"

       pi = 3.14159265358979323846D0

       flag1 = .true.
       nvertex = 1

       lambda =llambda   !5d0

      np =np3 ! 128
      L = dlog(10d0)/20d0   !  log(10d0)/m
      xmin = 1d-10
c       xmin =1d-6
      umin = dlog(1d0/(1d0-exp(-L)))
      umax = dlog(1d0/(1d0-exp(-L))/xmin)
      tol1 = 0.0001d0

       call gauleg(-1d0,1d0,xp,wp,np)
     
     
      do 100 i = 1, np
c      u = umin*(np-i*1d0)/(np-1d0)+umax*(i-1d0)/(np-1d0)
c      pp(i) =lambda*xmin*dexp(u)*(1d0-dexp(-L))
c      pp(i) = lambda*(i/(np*1d0))**(1.5d0)   
      pp(i) = (1d0+xp(i))/(1d0-xp(i))
      ppr(i)= pp(i)   
c      AA(i)= alatt(pp(i)**2) !  0.5d0       
c      BB(i)= blatt(pp(i)**2) !  1d0
      AAin(i)= AA(i)
      BBin(i)= BB(i)
c      write(3,*)ppr(i),BBin(i)
c      write(4,*)ppr(i),AAin(i)
      CC(i) = 1d0
      CCin(i)= CC(i)
100   continue
      pmin2 = 1d-14      !pp(1)**2
     

      


      flagsalir = .false.
      maxit =  100
      j = 1 
c      z2r = 1d0
c      z4r = 1d0

 30   IF ((.FALSE..EQV.flagsalir).and.(j.le.maxit)) THEN

c       mu = 4.3d0  ! imput
          
        ikernel = 2
         p = mu
         p2= p*p
       call  intsder(flag1,nvertex,ikernel,np,pp,aain,
     .   bbin,ccin,p2,lambda,pmin2,z2r,A19)
       z2r =(Amu-A19)




       if(finfty.eqv..true.)then
       p = pinfty
       p2 = p*p

       ikernel = 1

       call  intsder(flag1,nvertex,ikernel,np,pp,aain,
     .   bbin,ccin,p2,lambda,pmin2,z2r,B19)
        

        ikernel = 2
        call  intsder(flag1,nvertex,ikernel,np,pp,aain,
     .  bbin,ccin,p2,lambda,pmin2,z2r,A19)


       write(2,*)"minfty,z2r,A19,B19,mmu",minfty,z2r,A19,B19,mmu  
       z4r = (minfty*(z2r+A19)-B19)/mmu   !(mmu-B19)/mmu



       elseif(finfty.eqv..false.)then
       p = mu
       p2 = p*p
        ikernel = 1
        call  intsder(flag1,nvertex,ikernel,np,pp,aain,
     .   bbin,ccin,p2,lambda,pmin2,z2r,B19)
        z4r = (Bmu-B19)/mmu   !(mmu-B19)/mmu
       else
       endif

c       ikernel = 2
c         p = mu
c         p2= p*p
c       call  intsder(flag1,nvertex,ikernel,np,pp,aain,
c     .   bbin,ccin,p2,lambda,pmin2,z2r,A19)
c       z2r =(Amu-A19)

      do 105 i = 1,np
      p = pp(i)
      p2 = p*p
      ikernel = 1

      call  intsder(flag1,nvertex,ikernel,np,pp,aain,
     .   bbin,ccin,p2,lambda,pmin2,z2r,sigmai)

      BB(i) =z4r*mmu+sigmai
105   continue     
      
   
      do 106 i = 1,np
      p = pp(i)
      p2 = p*p
      ikernel = 2

      call  intsder(flag1,nvertex,ikernel,np,pp,aain,
     .   bbin,ccin,p2,lambda,pmin2,z2r,sigmai)

      AA(i) = z2r+sigmai
106   continue


c      do 110 i = 1,np
c      p = pp(i)
c      p2 = p*p
c      ikernel = 3
c
c       call  intsder(flag1,nvertex,ikernel,np,pp,aain,
c     .   bbin,ccin,p2,lambda,pmin2,z2r,sigmai)
c
c      CC(i) = 1d0+sigmai
c110   continue




      errormax = abs((bb(1)-bbin(1))/bbin(1))

      do 107 i = 1,np
      error = abs((bb(i)-bbin(i))/bbin(i))
      if(error.gt.errormax) errormax = error
      error= abs((aa(i)-aain(i))/aain(i))
      if(error.gt.errormax) errormax = error
c      error = abs((cc(i)-ccin(i))/ccin(i))
c      if(error.gt.errormax) errormax = error
107   continue
      
      write(2,*)'errormax', errormax

     

      write(2,*)'B19,A19,before-loop',b19,a19
      call massf(np,pp,AAin,mu,sigmai)
      write(2,*)'Amu-in-after-loop',sigmai
      call massf(np,pp,BBin,mu,sigmai)
      write(2,*)'Bmu-in-alfer-lopp',sigmai
      write(2,*)'z2r,z4r',z2r,z4r
      call massf(np,pp,BB,mu,sigmai)
      write(2,*)'Bmu',sigmai 
      call massf(np,pp,AA,mu,sigmai)
      write(2,*)'Amu',sigmai


   

        
      do 108 i = 1,np
      BBin(i) = BB(i)
      AAin(i) = AA(i)
      CCin(i) = CC(i)
108   continue
     
      write(2,*)'errormax', errormax
      do 111 i = 1,np
c      write(5,*)pp(i)**2,bb(i),BB(i)/AA(i),1d0/AA(i),cc(i)
c     . g2(pp(i)**2)/(pp(i)**2*4*3.1415d0*0.295),Fgh(pp(i)**2)
111   continue
       

       write(2,*)'tol12',tol1
      IF(errormax.LT.tol1)flagsalir = .true.
       j=j+1
       write(2,*)'tol1',tol1
       write(2,*)'flagsalir',flagsalir
       GOTO 30
      ElSEIF(.TRUE..EQV.flagsalir) THEN
      ENDIF
      write(2,*)'j',j
      write(3,*)'sde new'
      write(4,*)'sde new'
      do 109 i = 1,np
      write(3,*)pp(i),BB(i)
      write(4,*)pp(i),AA(i)
c      write(5,*)pp(i)**2,CC(i) 
c      write(7,*)pp(i)**2,CC(i)*Fgh(pp(i)**2)**2*g2(pp(i)**2)/(4d0*pi)     
109   continue

c      do 113 i = 1,np
c      write(5,*)pp(i)**2,CC(i)
c113   continue

      end

