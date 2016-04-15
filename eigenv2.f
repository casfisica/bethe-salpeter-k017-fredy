      subroutine eigenv2(Mn,nres,func)
      implicit none
      include 'common.f'
      include 'common2.f'
      double precision Mn,errormax, error,tol,func
      double complex sum,lp,lpmax,lpmax0,F0(4,0:4,128)
      integer  i,j,maxit,k,l,jn,jm,it,nx,nres
      logical flagsalir

      external kernel
      
      call kernel(Mn,nx)
      tol = 0.0000001d0


      do 157 i = 1,8
      do 158 jm= 0,nm
      do 159 k =1,ns

      p=  (lambda-dsqrt(pmin2))*xs(k)+dsqrt(pmin2)
      F0(i,jm,k) =p*exp(-p)
      F(i,jm,k)=F0(i,jm,k)
c      write(2,*)F0(i,jn,k)
159   continue
158   continue
157   continue

      lp =1d0
      lpmax =1d0
      lpmax0=0.2d0
      
      flagsalir = .false.

      maxit = 100
      it = 1

 30   IF ((.FALSE..EQV.flagsalir).and.(it.le.maxit)) THEN




      do 143 i =1,8 !4
      do 144 jm=0, nm
      do 145 k=1,ns

      sum =0d0

      do 146 j =1,8
      do 147 jn=0,nm
      do 148 l=1,ns

       

      sum =(kbsr(i,j,jm,jn,k,l)+(0d0,1d0)*kbsi(i,j,jm,jn,k,l))
     . *ws(l)*2d0*((lambda-dsqrt(pmin2))*xs(l)+dsqrt(pmin2))**3
     . *(lambda-dsqrt(pmin2))*F0(j,jn,l)+sum           !/1.51985d0+sum  !  ! 3.64811364d0
      
c       write(2,*)sum

148   continue
147   continue
146   continue
     
c      write(2,*)sum   

      F(i,jm,k)= sum
 
145   continue
144   continue
143   continue
    

       

      lpmax =  abs(F(1,0,1)/F0(1,0,1))
       
      errormax = abs((F(1,0,1)-F0(1,0,1))/F0(1,0,1))

 







      do 151 i = 1, 8  ! dont't change
      do 152 jm= 0, nm  ! "

      do 153 k =1,ns

      lp = abs((F(i,jm,k))/F0(i,jm,k))
      if(abs(lp).gt.abs(lpmax)) lpmax =lp



      error = abs((F(i,jm,k)-F0(i,jm,k))/F0(i,jm,k))
      if(error.gt.errormax) errormax =error 
153   continue

152   continue
151   continue

      

       
       
      do 154 i = 1,8
      do 155 jm= 0,nm 

      do 156 k =1,ns

      F0(i,jm,k) = F(i,jm,k) 

      
156   continue

155   continue
154   continue 
      
        write(30,*)'errormax',errormax
       write(30,*)'lpmax,lp',lpmax,lp
       write(30,*)'it',it
     
      

       
c       IF(error.LT.tol)flagsalir = .true.
       IF(abs((lpmax-lpmax0)/lpmax0).LT.tol)flagsalir = .true.
        lpmax0 = lpmax
       it=it+1
       GOTO 30
      ElSEIF(.TRUE..EQV.flagsalir) THEN
      ENDIF



      do 149 i =1,8
      do 150 jm=0,nm
      write(30,*)'F alpha,n',i,jn
      write(30,*)(F(i,jm,k)/F(1,0,1),k=1,ns)
150   continue
149   continue
       write(30,*)'eigenv values'   
       write(30,*)'lpmax,lpmax0',lpmax,lpmax0
       write(*,*)'lpmax,lpmax0',lpmax,lpmax0
       write(30,*)'it',it
       write(30,*)"F(1,0,1)",F(1,0,1)
       write(30,*)'Mn',Mn

       func =dble(lpmax)-1d0   

       return
       end
