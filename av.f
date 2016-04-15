       subroutine av(nx,v,w)
       implicit none
       include 'common2.f'
       integer nx,i,j
       double complex  v(nx), w(nx)
c       double precision kijr(3000,3000),kiji(3000,3000)


c       open (24,file='kbsr.dat',status='unknown')
c       open (25,file='kbsi.dat',status='onknown')




       do 101 i=1,nx
       w(i)  =0d0
       do 100 j=1,nx
c       w(j)=v(j)
       w(i)= (kijr(i,j)+(0d0,1d0)*kiji(i,j)+
     .      0d0*(kijr(j,i)-0d0*(0d0,1d0)*kiji(j,i))
     .        )*v(j)+w(i)
100    continue
101    continue
       
       return
       end
