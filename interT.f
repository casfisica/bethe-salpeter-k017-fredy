      subroutine interT(n,p,FT,FI,q,F)
      implicit none
      integer n,i,ia,j
      double precision p(64),FT(8,0:5,64),q,FI(8,0:5,64),lcut
      double complex F(8,0:5)

       ia = 1
c       n = 64
       do 1 i = 1, n-1
       if((p(i).le.q).and.(q.lt.p(i+1)))ia = i
1      continue
        
       lcut = 200d0

       do 2 i= 1,8
       do 3 j= 0,5
       F(i,j)= 0d0
       F(i,j) = 
     . FT(i,j,ia)+(FT(i,j,ia+1)-FT(i,j,ia))/(p(ia+1)-p(ia))*(q-p(ia))
       F(i,j) = (FI(i,j,ia)+(FI(i,j,ia+1)-FI(i,j,ia))/
     . (p(ia+1)-p(ia))*(q-p(ia)))*(0d0,1d0)+F(i,j)
 
       if(q.gt.lcut)F(i,j)=0d0

3      continue
2      continue


         if(p(n).le.q) then

        do 4 i= 1,8
        do 5 j= 0,5
        F(i,j) = FT(i,j,n-1) !+(FT(i,j,n)-FT(i,j,n-1))/
c     .   (p(n)-p(n-1))*(q-p(n-1))
c        F(i,j) = (FI(i,j,n-1)+(FI(i,j,n)-FI(i,j,n-1))/
c     .   (p(n)-p(n-1))*(q-p(n-1)))*(0d0,1d0)+F(i,j) 
        if(q.gt.lcut)F(i,j)=0d0
5      continue
4      continue
       

        elseif(p(1).ge.q)then
         
        do 6 i= 1,8
        do 7 j= 0,5
        F(i,j) = FT(i,j,1)+(FT(i,j,2)-FT(i,j,1))/
     .   (p(2)-p(1))*(q-p(1))
        F(i,j) = (FI(i,j,1)+(FI(i,j,2)-FI(i,j,1))/
     .   (p(2)-p(1))*(q-p(1)))*(0d0,1d0)+F(i,j)
          if(q.gt.lcut)F(i,j)=0d0
7      continue
6      continue
 
      

         endif

   
      return
      end
