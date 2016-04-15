


      double precision FUNCTION ALatt(x)
      implicit none
     
      double precision  x,aux,z,m2,c
      


      z = 0.8443d0
      m2 = 0.395d0
      c = 0.1235d0

         aux = (1d0 +dLOG(x + m2))**c
         aux = z*aux

         ALatt = 1d0/aux

      return
      end




      double  precision FUNCTION BLatt(x)
      implicit none
      double precision  x,aux,bux,cux,z,m2,c,mq0,mq2,mq3,A,rho,gammam
      
      
       z = 0.8443d0
       m2 = 0.395d0
       c = 0.1235d0
       mq0 = 0.0343d0
       mq2 = 0.311d0
       mq3 = 0.122d0 
       A = 0.294d0 
       rho = 38.9d0
       gammam = 0.413793d0

         aux = (1d0 +dLOG(x + m2))**c
         aux = z*aux
         aux = 1d0/aux

         bux = (x + mq2)
         bux = 1d0/bux
         bux = mq3*bux + mq0

         cux = A + dLOG( x + rho*bux*bux )
         cux = cux**gammam
         cux = 1d0/cux

         bux = bux*cux

         BLatt = bux*aux 

       return
       end




       double precision  FUNCTION mqLatt(x)

       double precision  x, bux, cux,m0,m2,m3,A,rho,gammam

        m0 = 0.0343d0
        m2 = 0.311d0
        m3 = 0.122d0
        A = 0.294d0 
        rho = 38.9d0 
        gammam = 0.413793d0

         bux = (x + m2)
         bux = 1d0/bux
         bux = m3*bux + m0

         cux = A + dLOG( x + rho*bux*bux )
         cux = cux**gammam
         cux = 1d0/cux

         mqLatt = bux*cux
       return 
       end




       double precision FUNCTION deltaglue(q)

       IMPLICIT NONE

       double precision q,  m2glue, aux, deltainverse,m,m2,m4,mu,mu2,
     . rho1,rho2,g21

        m = 0.520d0
        m2 = m*m 
        m4 = m2*m2
        mu = 4.3d0
        mu2 = mu*mu
        rho1 = 8.55d0
        rho2 = 1.91d0
        g21 = 0.2337986312547d0 

      m2glue = m4/(q + rho2*m2)

      aux = (q + rho1*m2glue)/mu2
      aux = DLOG(aux)
      aux = 1d0 + g21*aux

      deltainverse = m2glue + q*aux

      deltaglue = 1d0/deltainverse

      return
      end





      doubleprecision FUNCTION Fgh(q) 

      IMPLICIT NONE

      doubleprecision  q,m2ghost,aux, Finverse,m,m2,m4,mu,mu2,
     . rho3,rho4,g22

       m = 0.520d0
       m2 = m*m 
       m4 = m2*m2 
       mu = 4.3d0 
       mu2 = mu*mu 
       rho3 = 0.25d0 
       rho4 = 0.68d0
       g22 = 0.12210785772396d0

      m2ghost = m4/(q + rho4*m2)

      aux = (q + rho3*m2ghost)/mu2
      Finverse = 1d0 + g22*DLOG(aux)

      Fgh = 1d0/Finverse
c      fgh = fgh**2
      return
      end
 


      double precision function H1(k2)
      implicit none
      double precision k2, c,a,k,b,w
      k = dsqrt(k2)
      c = 1.26d0
      a = 0.80d0 ! GeV
      b = 1.3d0  ! GeV
      w = 0.65d0 ! GeV

      H1 = c*(1d0+a**2*k**2/(k**4+b**4))
     .   + (1d0-c)*w**4/(w**4+k**4)


      return
      end

