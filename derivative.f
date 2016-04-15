c      FUNCTION derivative(func,x,h,err,nres)
c      subroutine derivative(print,mn,ppf,FR,FI,h,err,dfridr)
      subroutine derivative(dfridr,x,h,err,nres)
      implicit none
      INTEGER NTAB,nres
      double precision dfridr,err,h,x,func,CON,CON2,BIG,SAFE
      double precision mn,xph,xmh
      PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=12,SAFE=2.)
      logical print
c      EXTERNAL func
CU    USES func
      INTEGER i,j
      double precision errt,fac,hh,a(NTAB,NTAB)


c      x =mn

      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
      call eigenv(x+hh,nres,xph)
      call eigenv(x-hh,nres,xmh)
c     a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      a(1,1)=(xph-xmh)/(2.0*hh)
      write(2,*)"2h,a(1,1)",2d0*hh,a(1,1)
      write(*,*)"2h,a(1,1)",2d0*hh,a(1,1)  
      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
      call eigenv(x+hh,nres,xph)
      call eigenv(x-hh,nres,xmh)
c     a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
      a(1,i)=(xph-xmh)/(2.0*hh)
      write(2,*)"2h,a(1,1)",2d0*hh,a(1,i)
      write(*,*)"2h,a(1,1)",2d0*hh,a(1,i)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfridr=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .
