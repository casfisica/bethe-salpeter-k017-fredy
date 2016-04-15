       double precision function u(n,x)
       implicit none
       double precision x
       integer n
       ! chebyshev second kind
       !u= dsin((1d0+n)*dacos(x))/dsqrt(1d0-x*x)
       if(n.eq.0)u=1d0
       if(n.eq.1)u=2d0*x
       if(n.eq.2)u=4d0*x**2-1d0
       if(n.eq.3)u=8d0*x**3-4d0*x
       if(n.eq.4)u=16d0*x**4-12d0*x**2+1d0
       if(n.eq.5)u=32d0*x**5-32d0*x**3+6d0*x
       if(n.eq.6)u=64d0*x**6-80d0*x**4+24d0*x**2-1d0
       return
       end

