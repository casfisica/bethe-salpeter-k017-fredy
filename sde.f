      subroutine sde(mn,mu,Bmu,Amu,mmu,llambda,np3,cc0,z2r,z4r,
     .                  bbijr,bbiji,aaijr,aaiji)
      implicit none
      include 'commonsde.f'
      double precision sigmai(2),tol1,xmin,u,umin,umax,errormax,
     . error,mmu,mu,g2,Fgh,alatt,blatt,ppr(65),
     . aar(65),bbr(65),a,b,Bmu,Amu,kr,z1,k2r,p2r,pkr,q2,mathcalg,
     . aaijr(np3,3),aaiji(np3,3),bbijr(np3,3),bbiji(np3,3),
     . ppij2r(128,3),ppij2i(128,3),mn,llambda,cc0(128),
     . xpa(128),wpa(128),xpb(128),wpb(128),xpc(128),wpc(128),ab(2),
     . lamb_bs

      
      double complex A19,B19,z2r,z4r,ppij2(128,3),tr
        
      integer i,j,maxit,l,np3,nab(3)
      logical flagsalir
      external intsde
      external massf
      external alatt
      external blatt
      external Fgh
    

c      open (1,file="salida1.out",status="new")
      open (2,file="salida2.out",status="unknown")
c      open (3,file="salida3.out",status="new")
c      open (4,file="salida4.out",status="new")
c      open (5,file="salida5.out",status="new")
c      open (7,file="salida7.out",status="new") ! avoid open (6,
c      open (8,file="aa.dat",status="old")
c      open (9,file="bb.dat",status="old")
c      open (10,file="salida10.out",status="new") ! avoid open (6,
c      open (11,file="salida11.out",status="new")
c      open (12,file="salida12.out",status="new")
c      open (13,file="salida13.out",status="new")
c      open (14,file="salida14.out",status="new")
c      open (15,file="salida15.out",status="new")
c      open (16,file="salida16.out",status="new
   
  



       pi = 3.14159265358979323846D0   


       print = .true.
       linear= .true.
       flag2 = .false. 

       nvertex = 1

       lambda =llambda    ! dsqrt(10d0)
c       mu = 0.99d0*lambda
       eplus = 1d0/2d0
       eminus= 1d0/2d0
       mh = mn       !0.4977d0    ! 0.1385d0
       lamb_bs = dsqrt(lambda**2-mh**2/4d0)
 
      np = np3   !64
      np2 =np3      !64
      nzq = 128
      z2g = 1d0

      call  quadrature(lamb_bs,xp,wp,np)
      call  gauleg(0d0,1d0,xv,wv,np2)     
      call  gauleg(0d0,1d0,xzq,wzq,nzq) 


      tol1 = 0.001d0

   





      pmin2=1d-14

     
c      do 101 j= 1,3
      do 100 i = 1, np
      if(linear.eqv..false.)then
      pp(i) = dsqrt((lambda**2-pmin2)*xp(i)+pmin2)
      ppv(i) = (eplus+eminus)*mh*xv(i)-eminus*mh
      else
      pp(i) = (lambda-dsqrt(pmin2))*xp(i)+dsqrt(pmin2)
      ppv(i) = (eplus+eminus)*mh*xv(i)-eminus*mh
      endif
100   continue
c101   continue

       do 114 i =1,np
      ccin(i)=cc0(i)
115   continue

114    continue       

          
       do 102 j=1,3
       do 103 i=1,np
       BB(i,j)= BBijr(i,j)+(0,1d0)*BBiji(i,j)
       BBin(i,j)= BB(i,j)
       AA(i,j)= aaijr(i,j)+(0,1d0)*aaiji(i,j)
       AAin(i,j)= AA(i,j)
103    continue
102    continue

       write(2,*)"lambda,np,mu,bmu",lambda,np,mu,bmu

       


      
     
c      if(flag2.eqv..false.)goto 40
 
      flagsalir = .false.

      maxit =  20
      l = 1
      write(*,*)'sde.f'
      write(2,*)'sde.f'
 
 30   IF ((.FALSE..EQV.flagsalir).and.(l.le.maxit)) THEN

      write(2,*)"iteration",l
      write(*,*)"iteration, errormax",l,errormax
     
       write(2,*)'Bmu',Bmu
      
       ikernel = 1
       p= mu            !0.99d0*lambda
       pr = mu          !0.99d0*lambda
c       mmu = Bmu
       pc4 =0d0
       call intsde(sigmai)
       B19 = sigmai(1)+(0,1d0)*sigmai(2)
c       z4r = (Bmu-B19)/mmu
       write(2,*)"B19,z4r",B19,z4r
       ikernel = 2
       call intsde(sigmai)
       A19 = sigmai(1)+(0,1d0)*sigmai(2)
c       z2r =Amu-A19
c       z2g = z2r
        z2g = z2r
       write(2,*)'z2r,A19',z2r,A19

       
      do 110 j =1,3
      contour = j
      do 105 i = 1,np
      
      if(j.eq.1)then
      pc4 = -(0d0,1d0)*eminus*mh
      p = pp(i)+pc4
      ppij2(i,j)=p*p
      ppij2r(i,j)=dble(p*p)
      ppij2i(i,j)=dimag(p*p) 
      pr = pp(i)
      elseif(j.eq.2)then
      pc4 = (0d0,1d0)*ppv(i)
      p =lambda+pc4
      ppij2(i,j)=p*p
      ppij2r(i,j)=dble(p*p)
      ppij2i(i,j)=dimag(p*p)
      pr =lambda
      elseif(j.eq.3)then
      pc4 = (0d0,1d0)*eplus*mh
      p = pp(i)+pc4
      pr = pp(i)
      ppij2(i,j)=p*p
      ppij2r(i,j)=dble(p*p)
      ppij2i(i,j)=dimag(p*p)
      else
      endif

      if((j.eq.1).or.
     . (  (j.eq.2).and.(i.le.64)  ))then 

      ikernel = 1
      call intsde(sigmai)
      BB(i,j) =z4r*mmu+sigmai(1)+(0,1d0)*sigmai(2)
      
      ikernel = 2
      call intsde(sigmai)
      AA(i,j) = z2r+sigmai(1)+(0,1d0)*sigmai(2)
   
      elseif(j.eq.3)then
      BB(i,j) = dble(BB(i,1))-(0d0,1d0)*dimag(BB(i,1))
      AA(i,j) = dble(AA(i,1))-(0d0,1d0)*dimag(AA(i,1))
      elseif((j.eq.2).and.(i.gt.64))then
      BB(i,j) = dble(BB(128-i+1,2))-(0d0,1d0)*dimag(BB(128-i+1,2))
      AA(i,j) = dble(AA(128-i+1,2))-(0d0,1d0)*dimag(AA(128-i+1,2))
      else
      endif

      
       
   

105    continue
       write(2,1000)(dble(BB(i,j)),dimag(BB(i,j)),
     . pp(i),j,i=1,np)
      write(2,1001)(dble(AA(i,j)),dimag(AA(i,j)),
     . pp(i),j,i=1,np)
110    continue


1000  format('B,pr,j','(',D18.9,D18.9,')',D18.9,i2)
1001  format('A,pr,j','(',D18.9,D18.9,')',D18.9,i2)



      errormax = abs((bb(1,1)-bbin(1,1))/bbin(1,1))

      do 106 j = 1,3
      do 107 i = 1,np
      error = abs((bb(i,j)-bbin(i,j))/bbin(i,j))
      if(error.gt.errormax) errormax = error
      error= abs((aa(i,j)-aain(i,j))/aain(i,j))
      if(error.gt.errormax) errormax = error
107   continue
106   continue
      
      write(2,*)'errormax', errormax


      do 112 j = 1,3  
      do 108 i = 1,np
      BBin(i,j) = BB(i,j)
      AAin(i,j) = AA(i,j)
108   continue
112   continue
     
c      write(3,*)'errormax', errormax
c      do 111 i = 1,np
c      write(3,*)pp(i),bb(i,j),BB(i,j)/AA(i,j),1d0/AA(i,j)
c111   continue
       

       IF(errormax.LT.tol1)flagsalir = .true.
       l=l+1
       GOTO 30
      ElSEIF(.TRUE..EQV.flagsalir) THEN
      ENDIF
      write(2,*)'j',j
      
     
      write(14,*)"Bmu",bmu
      write(15,*)"Bmu",bmu
      write(16,*)"Bmu",bmu
      write(17,*)"Bmu",bmu
300    format(D25.16,D25.16,D25.16)            

      do 113 j = 1,3
      do 109 i = 1,np  
      write(14,300)dble(ppij2(i,j)),dimag(ppij2(i,j)),
     . dble(BB(i,j))
       bbijr(i,j)=dble(BB(i,j))
       write(15,300)dble(ppij2(i,j)),dimag(ppij2(i,j)),
     . dimag(BB(i,j))
       bbiji(i,j)=dimag(BB(i,j))
       write(16,300)dble(ppij2(i,j)),dimag(ppij2(i,j)),
     . dble(AA(i,j)) 
       aaijr(i,j)=dble(aa(i,j))   
       write(17,300)dble(ppij2(i,j)),dimag(ppij2(i,j)),
     . dimag(AA(i,j))
       aaiji(i,j)=dimag(aa(i,j))
109    continue
113    continue

        
c      call cplotter      
   
      return
      end

