       subroutine kernel(Mn,nx)
       implicit none
       include 'common.f'
       include 'common2.f'
       double precision mn,zq,sigmai(2),jaq
       double complex  pc4,sigmasp,
     .  sigmasm,sigmavp,sigmavm
       integer i,j,jn,jm,l,k,ia,ib,nx
       external scomplex
       external kernint
 
       write(*,*)Mn
       write(2,*)"kernel"
       do 121 i =1,ns
       do 122 j= 1,nt
       if(ffit.eqv..true.)then
       q2E =  (  (1d0+xs(i))/(1d0-xs(i)))**2
       else
       q2E=  ((lambda-dsqrt(pmin2))*xs(i)+dsqrt(pmin2))**2
       endif
       pph(i)= dsqrt(q2E)
       zq = xt(j)
       if(ffit.eqv..true.)then
       pc4 =  (0d0,1d0)*mn*mh/(mh+mh2)
       call  scomplexfit(print,pc4,q2E,zq,xval,npol,sigmavp,sigmasp)
       pc4 = -(0d0,1d0)*mn*mh2/(mh+mh2)
c       print = .true.
       call scomplexfit(print,pc4,q2E,zq,xval2,npol2,sigmavm,sigmasm)
       else
       pc4 =  (0d0,1d0)*mn*mh/(mh+mh2)
       call  scomplex(pc4,q2E,zq,mh,AA,BB,sigmavp,sigmasp)
       pc4 = -(0d0,1d0)*mn*mh2/(mh+mh2)
       call  scomplex(pc4,q2E,zq,mh2,AA2,BB2,sigmavm,sigmasm)
       endif

       ssigmavp(i,j)=sigmavp
       ssigmasp(i,j)=sigmasp
       ssigmavm(i,j)=sigmavm
       ssigmasm(i,j)=sigmasm

       write(2,*)"B-plus"
       write(2,*)i,j,
     . sigmasp/(sigmasp**2 + q2E*sigmavp**2)
c       write(2,*)"A-plus"
c       write(2,*)i,j,              
c     . sigmavp/(sigmasp**2 + p2*sigmavp**2)

 
 122   continue
 121   continue
      
      
      do 114 i=1,8
      do 115 j=1,8
      init = .true.
      write(*,*)"Mn,i,j",mn,i,j
      write(2,*)"init",init
      write(2,*)"Mn,i,j",mn,i,j
       print *,"Mn,i,j",mn,i,j
       print *, "init",init 
c      write(22,*)"init",init
c      write(22,*)"Mn,i,j",mn,i,j
      do 116 jm=0,nm
      do 117 jn=0,nm
      m=jm
      n=jn

      do 118 k=1,ns
      do 119 l=1,ns
      alpha =i
      beta = j 
      jp=k
      jq=l
c      p=    (lambda-dsqrt(pmin2))*xs(k)+dsqrt(pmin2)
c      q2E=  ((lambda-dsqrt(pmin2))*xs(l)+dsqrt(pmin2))**2
       p = pph(k)
      q2E = (pph(l))**2
      call  kernint(mn,sigmai)
     
      kbsr(i,j,jm,jn,k,l)=sigmai(1)
      kbsi(i,j,jm,jn,k,l)=sigmai(2)

c       read(24,*)kbsr(i,j,jm,jn,k,l)
c       read(25,*)kbsi(i,j,jm,jn,k,l)




c      if(nt.eq.30)then
      write(22,*)kbsr(i,j,jm,jn,k,l)
      write(23,*)kbsi(i,j,jm,jn,k,l)
c      else
c      endif

c      write(*,*)kbsr(i,j,jn,jm,k,l)

119   continue
118   continue

117   continue
116   continue

115   continue
114   continue










      Ia=0
      Ib=0
      do 135 i=1,8
      do 137 jm=0,nm 
      do 139 k=1, ns
      ia=ia+1
 
 
      Ib=0      

   
      do 136 j=1,8
      do 138 jn=0,nm  
      do 140 l=1, ns
      Ib=ib+1
      
 


      m= jm
      n= jn
      alpha =i
      beta = j
      jp=k
      jq=l
c      p=     (lambda-dsqrt(pmin2))*xs(k)+dsqrt(pmin2)
c      q2E=  ((lambda-dsqrt(pmin2))*xs(l)+dsqrt(pmin2))**2

       p = pph(k)
       q2E = pph(l)**2
 
       if(ffit.eqv..true.)then
       jaq = 2d0/(1d0-xs(l))**2 
       else
       jaq = (lambda-dsqrt(pmin2))
       endif

      kijr(Ia,Ib)=kbsr(i,j,jm,jn,k,l)
c     . *ws(l)*2d0*((lambda-dsqrt(pmin2))*xs(l)+dsqrt(pmin2))**3
     . *ws(l)*2d0*(pph(l))**3
     . *jaq


      kiji(Ia,Ib)=kbsi(i,j,jm,jn,k,l)
c     . *ws(l)*2d0*((lambda-dsqrt(pmin2))*xs(l)+dsqrt(pmin2))**3
     . *ws(l)*2d0*(pph(l))**3
     . *jaq


      write(26,*) kbsr(i,j,jm,jn,k,l)
      write(27,*) kbsi(i,j,jm,jn,k,l)
140   continue
138   continue
136   continue

139   continue
137   continue
135   continue

       write(2,*)'kernel'
       write(2,*)'mn',mn  
       write(2,*)"ia,ib",ia,ib
      do 141 i =1,Ia
      do 142 j =1,Ib     
      write(28,*) kijr(i,j)
      write(29,*) kiji(i,j)
142   continue
141   continue


      nmp = nm
      nx = Ia







       return
       end 
