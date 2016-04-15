c      FUNCTION rtsec(func,x1,x2,nres,xacc)
       subroutine sec(x1,x2,xacc,nres,rtsec)
      INTEGER MAXIT,nres
      double precision rtsec,x1,x2,xacc
C,func
c      EXTERNAL func
      EXTERNAL eigenv 
      PARAMETER (MAXIT=30)
      INTEGER j
      double precision dx,f,fl,swap,xl
c      fl=func(x1)
c      f=func(x2)
       call eigenv(x2,nres,f)
       call eigenv(x1,nres,fl)
    

      if(abs(fl).lt.abs(f))then
        rtsec=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
      else
        xl=x1
        rtsec=x2
      endif
      do 11 j=1,MAXIT
        dx=(xl-rtsec)*f/(f-fl)
        xl=rtsec
        fl=f
        rtsec=rtsec+dx
c        f=func(rtsec)
         call eigenv(rtsec,nres,f)
        if(abs(dx).lt.xacc.or.f.eq.0.)return
11    continue
      pause 'rtsec exceed maximum iterations'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .


c       subroutine eigenv(x,nres,f)
c       double precision x,f
c       f= (x-5)**3
c       return
c       end
 
c       program main
c       external sec
c       double precision rtsec
c       call sec(4d0,8d0,0.0001d0,rtsec)
c       write(*,*)rtsec
c       end

