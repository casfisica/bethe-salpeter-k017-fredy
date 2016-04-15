      subroutine massf(n,p,M,q,B)
      implicit none
      integer n,i,jj
      double precision p(n),M(n),q,B

       do 1 i = 1, n-1
       if((p(i).le.q).and.(q.lt.p(i+1)))
     .  B = m(i)+(m(i+1)-m(i))/
     . (p(i+1)-p(i))*(q-p(i))
1      continue

         if(p(n).le.q) then
        B = m(n-1) !+(m(n)-m(n-1))/
c     .   (p(n)-p(n-1))*(q-p(n-1))
         elseif(p(1).ge.q)then
        B = M(1)+(M(2)-M(1))/
     .   (p(2)-p(1))*(q-p(1))
         endif


  
     


   
      return
      end
