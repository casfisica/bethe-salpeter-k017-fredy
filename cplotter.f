      subroutine cplotter
      implicit none 
      include 'commonsde.f'
      include 'commoncmaker.f'  
      double precision kr,z1,k2r,q2,g2,mathcalg,level(4,25)
      double complex tr
      integer npt,nh,i,j,id,iu,l,nfun,ncont
      external gauleg,g2,vertex,cmaker

      npt =100 !300
      nh = 64 !2000
   

      call  gauleg(0d0,1d0,xpt,wpt,npt) 
      call  gauleg(0d0,1d0,xh,wh,nh)


       write(2,*)"vertex"      
       do 114 i=1,nh
       ia=i
       do 116 l=1,npt
       ib=l
       kr  = (lambda-dsqrt(pmin2))*xh(i)+dsqrt(pmin2)
       pc4 = (0d0,1d0)*((eplus+eminus)*mh*xpt(l)-eminus*mh)
       p =kr+pc4


       z1  =1d0 ! 2*xv(l)-1d0 
       k2r = kr*kr
c       pr =0d0 ! it doesn't matter
c       p2r = pr*pr
c       pkr = dsqrt(p2r*k2r)*z1
       q2  = k2r ! p2r+k2r-2*pkr



        mathcalg = g2(q2)
        
        flag2 = .true.    
        call vertex(k2r,q2,z1,mathcalg,tr)
116     continue
114     continue

       do 119 ncont=1,10          
c      call  cmaker(np,np,pr2v,ai,0.06d0,id,iu,cd,cu)
       nfun=4 
       level(nfun,ncont)=0.025d0*ncont 
       call  cmaker(nh,npt,level(nfun,ncont),nfun,id,iu)
       write(2,*)"cmaker ready"   
       write(*,*)"cmaker ready"
       write(18,*)"# vspace"
       Write(19,*)"# vspace"
119    continue
  
c       do 117 i=1,4  !  id
c       write(18,*)cd(i,1),cd(i,2)
c117    continue

          
 
  
c      do 118 i=1,iu
c      write(19,*)cu(i,1),cu(i,2)
c118   continue



      return
      end
