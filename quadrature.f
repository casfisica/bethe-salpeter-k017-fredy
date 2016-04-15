       subroutine  quadrature(lamb_bs,xp,wp,np)
       implicit none 
       integer np,nab(3),i
c       double precision lamb_bs,xpa(np),wpa(np),
c     . xpb(np),wpb(np),xpc(np),wpc(np),ab(2),xp(np),wp(np)
       double precision lamb_bs, xpa(128),wpa(128),
     . xpb(128),wpb(128),xpc(128),wpc(128),ab(2),xp(128),wp(128)
      

c      call  gauleg(0d0,1d0,xp,wp,np)
      ab(1)=  0.4d0/lamb_bs
      ab(2)=  1.5d0/lamb_bs
      if(np.le.64)then
      nab(1) =10 !20          ! any modification here also change in sde.f
      nab(2) =35     !128-nab(1)  ! 108
      nab(3) =np-nab(1)-nab(2)
      else
      nab(1) =20 !20          ! any modification here also change in sde.f
      nab(2) =70     !128-nab(1)  ! 108
      nab(3) =np-nab(1)-nab(2)
      endif


      write(2,*)"ab(1),ab(2)",ab(1),ab(2)
      write(2,*)"nab(1),nab(2),nab(3)",nab(1),nab(2),nab(3)
      call  gauleg(0d0  ,ab(1),xpa,wpa,nab(1))
      call  gauleg(ab(1),ab(2),xpb,wpb,nab(2))
      call  gauleg(ab(2),1d0  ,xpc,wpc,nab(3))
      do 124 i =1,nab(1)
      xp(i)= xpa(i)
      wp(i)= wpa(i)
124   continue
      if(nab(2).gt.0)then
      do 125 i =1,nab(2)
      xp(nab(1)+i)= xpb(i)
      wp(nab(1)+i)= wpb(i)
125   continue
      else
      endif
      if(nab(3).gt.0)then
      do 126 i =1,nab(3)
      xp(i+nab(1)+nab(2))= xpc(i)
      wp(i+nab(1)+nab(2))= wpc(i)
126   continue
      else
      endif
 
       return
       end
 
