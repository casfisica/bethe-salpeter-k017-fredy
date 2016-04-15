      subroutine eigenv(Mn,func)
      implicit none
      include 'common.f'
      double precision Mn,errormax, error,tol,func
      double complex sum,lp,lpmax,lpmax0,F(4,0:4,128),F0(4,0:4,128)
      integer  i,j,maxit,k,l,jn,jm,it
      logical flagsalir

      external kernel
      
      call kernel(Mn)
      tol = 0.00001d0


      do 157 i = 1,4
      do 158 jn= 0,nm
      do 159 k =1,ns

      p=  (lambda-dsqrt(pmin2))*xs(k)+dsqrt(pmin2)
      F0(i,jn,k) =p*exp(-p)
      F(i,jn,k)=F0(i,jn,k)
c      write(2,*)F0(i,jn,k)
159   continue
158   continue
157   continue

      lp =1d0
      lpmax =1d0
      lpmax0=0.2d0
      
      flagsalir = .false.

      maxit = 30
      it = 1

 30   IF ((.FALSE..EQV.flagsalir).and.(it.le.maxit)) THEN




      do 143 i =1,4 !4
      do 144 jn=0, nm
      do 145 k=1,ns

      sum =0d0

      do 146 j =1,4
      do 147 jm=0,nm
      do 148 l=1,ns

       

      sum =(kbsr(i,j,jn,jm,k,l)+(0d0,1d0)*kbsi(i,j,jn,jm,k,l))
     . *ws(l)*2d0*((lambda-dsqrt(pmin2))*xs(l)+dsqrt(pmin2))**3
     . *(lambda-dsqrt(pmin2))*F0(j,jm,l)+sum           !/1.51985d0+sum  !  ! 3.64811364d0
      
c       write(2,*)sum

148   continue
147   continue
146   continue
     
c      write(2,*)sum   

      F(i,jn,k)= sum
 
145   continue
144   continue
143   continue
    

       

      lpmax =  abs(F(1,0,1)/F0(1,0,1))
       
      errormax = abs((F(1,0,1)-F0(1,0,1))/F0(1,0,1))

 







      do 151 i = 1, 1 !4  ! dont't change
      do 152 jn= 0, 0 !nm  ! "

      do 153 k =1,ns

      lp = abs((F(i,jn,k))/F0(i,jn,k))
      if(abs(lp).gt.abs(lpmax)) lpmax =lp



      error = abs((F(i,jn,k)-F0(i,jn,k))/F0(i,jn,k))
      if(error.gt.errormax) errormax =error 
153   continue

152   continue
151   continue

      

       
       
      do 154 i = 1,4
      do 155 jn= 0,nm 

      do 156 k =1,ns

      F0(i,jn,k) = F(i,jn,k) 

      
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



c      do 149 i =1,4
c      do 150 jn=0,nm
c      write(30,*)'F alpha,n',i,jn
c      write(30,*)(F(i,jn,k)/F(1,0,1),k=1,ns)
c150   continue
c149   continue
       write(30,*)'eigenv values'   
       write(30,*)'lpmax,lpmax0',lpmax,lpmax0
       write(*,*)'lpmax,lpmax0',lpmax,lpmax0
       write(30,*)'it',it
       write(30,*)F(1,0,1)
       write(30,*)'Mn',Mn

       func =dble(lpmax)-1d0   

       return
       end
