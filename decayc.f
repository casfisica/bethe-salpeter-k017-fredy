       subroutine decayc(ffit_,mn,lambda_,mh_,mh2_,aa_,bb_,
     . aa2_,bb2_,xval_,npol_,xval2_,npol2_,z2,z2h,norm,fpi,rh,norm22)
                   
       implicit none
       include 'commonf.f' 
       include 'common2.f'
       double precision sigmai(2),FR(8,0:5,64),FI(8,0:5,64)
        double precision mn,norm,fpi,z2,z2h,mh_,mh2_,lambda_,rh,ab(2),
     .  xpa(128),wpa(128),xpb(128),wpb(128),xpc(128),wpc(128),lamb_bs,
     .  norm2,xval_(26),xval2_(26)


       double complex  AA_(128,3),BB_(128,3),AA2_(128,3),BB2_(128,3),
     . normc2,fpi2,rh2,norm22

       integer i,j,k,np2,nab(3),npol_,npol2_
       logical print,ffit_
       external int
c       open (38,file="FT.out",status="old") 
c       open (2,file="salidas2.out",status="unknown")

      do 101 i = 1,8
      do 102 j = 0,5        !nm chebyshev moments 
      do 103 k=1,64
      FR(i,j,k)=0d0
      FI(i,j,k)=0d0
      FR(i,j,k)=  dble(F(i,j,k)/F(1,0,1))
      FI(i,j,k)= dimag(F(i,j,k)/F(1,0,1))
103   continue
102   continue
101   continue



       print = .true. 
        

c      do 1 i = 1,4
c      do 2 j = 0,0        !nm
c      do 3 k=1,64
c      read(35,200)pph(k),FR(i,j,k),FI(i,j,k)
c      FI(i,j,k)= 0d0
c      if((j.eq.1).or.(j.eq.3))then
c      FR(i,j,k)=0d0 
c      FI(i,j,k)=0d0
c      else 
c      endif
c3     continue
c2     continue
c1     continue
 
      write(38,*)"decay.f"
      do 4 i = 1,8
      do 5 j = 0,4        !nmp
      do 6 k=1,64
      write(38,*)i,j,k,pph(k),FR(i,j,k),FI(i,j,k)
6     continue
5     continue
4     continue



      np = 128 ! to change it, modify bbrij & bb dimensions
      np2 =128
      nt = 20
      ns = 64
      lambda = lambda_
      lamb_bs = lambda
      write(2,*)"lambda",lambda
      pmin2 = 1d-14
      eplus = 1d0/2d0
      eminus= 1d0/2d0
      mh = mh_
      mh2= mh2_
      npol = npol_
      npol2= npol2_
      ffit = ffit_
      do 127 i = 1,26
      xval(i)  = xval_(i)
      xval2(i) = xval2_(i)
127   continue
      
      call quadrature(lamb_bs,xp,wp,np)
      call gauleg(0d0,1d0,xv,wv,np2)
      call gauleg(-1d0,1d0,xs,ws,ns)
      call gauleg(-1d0,1d0,xt,wt,nt)
 
      do 118 j = 1,3
      do 119 i = 1,np


       BB(i,j)= BB_(i,j)                  !BBf1r(i,j)+(0d0,1d0)*BBf1i(i,j)
       aa(i,j)= aa_(i,j)                  !aaf1r(i,j)+(0d0,1d0)*aaf1i(i,j)
 
       BB2(i,j)= BB2_(i,j)                 
       aa2(i,j)= aa2_(i,j)    

  

119   continue
118   continue


     
      
       ntr = 1
       print = .true.
       call intq(print,mn,pph,FR,FI,sigmai)
       
       norm = dsqrt(abs(3d0*sigmai(1))/(2d0*mn))
       write(2,*)"integral", sigmai(1),sigmai(2)       
       normc2 = 3d0*(sigmai(1)+(0,1d0)*sigmai(2))/(2d0*mn)
       write(2,*)"norm",norm
       write(2,*)"normc2",normc2
      
       ntr = 2 
       print = .true.
       call intq(print,mn,pph,FR,FI,sigmai)
       write(2,*)"integral 2", sigmai(1),sigmai(2)
       fpi = 3d0*sigmai(1)/(dsqrt(1d0)*norm*mn)        ! the 2 is for isospin
       fpi2= 9d0*(sigmai(1)+(0d0,1d0)*sigmai(2))**2/(1d0*normc2*mn**2)
       write(2,*)"fpi2",fpi2 
       write(2,*)"sqrt(abs(Re(fpi2)))",dsqrt(abs(dble(fpi2)))
       write(2,*)"decay constant", fpi
       write(2,*)"z2,z2h",z2,z2h

       ntr = 3
       
       print = .true.
       call intq(print,mn,pph,FR,FI,sigmai)
       write(2,*)"integral 3", sigmai(1),sigmai(2)
       rh = 3d0*sigmai(1)/(dsqrt(1d0)*norm)        ! the 2 is for isospin
       rh2= 9d0*(sigmai(1)+(0d0,1d0)*sigmai(2))**2/(1d0*normc2)
       write(2,*)"rh2",rh2
       write(2,*)"sqrt(abs(Re(rh2)))",dsqrt(abs(dble(rh2)))
       write(2,*)"rh", rh
       write(2,*)"z2,z2h",z2,z2h


       ntr = 4
       print = .true.
       call intq(print,mn,pph,FR,FI,sigmai)

       norm2 = dsqrt(abs(3d0*sigmai(1))/(2d0*mn))
       write(2,*)"integral", sigmai(1),sigmai(2)
       norm22 = -3d0*(sigmai(1)+(0,1d0)*sigmai(2))/(2d0*mn)
       write(2,*)"norm2",norm2
       write(2,*)"norm22",norm22



200    format(D25.16,D25.16,D25.16)
       return
       end
