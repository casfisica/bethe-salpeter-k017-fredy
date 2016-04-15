       subroutine cmaker(np,nz,hc,nfun,id,iu)
       implicit none
       include 'commoncmaker.f'
       integer i,j,id,iu,np,nz,nfun
       double precision b(2000,300),pv(2000,300,2)

       double precision hc,ee
       
       ee = abs(0.1d0*hc)


       do 102 i=1,np
       do 103 j=1,nz
       pv(i,j,1)=pr2v(i,j,1)
       pv(i,j,2)=pr2v(i,j,2)
       if(nfun.eq.4)b(i,j)=ai(i,j)
103    continue
102    continue 

 
       id=0
       iu=0

c       write(2,*)"cmaker"
       write(2,*)"pv"
       do  100 i =1,np 
 
       do  101 j =1,nz-1
        
       write(2,*)i,j,pv(i,j,2)

       if(   (b(i,j)  .lt.hc).and.
     .       (hc.le.b(i,j+1)).and.
     .       (abs(b(i,j)  -hc).le.ee).and.      
     .       (abs(b(i,j+1)-hc).le.ee)     )then
        id=id+1

        cd(id,1)=pv(i,j,1)+(pv(i,j+1,1)-pv(i,j,1))
     .  *(hc-b(i,j))/(b(i,j+1)-b(i,j))

      

        cd(id,2)=pv(i,j,2)+(pv(i,j+1,2)-pv(i,j,2))
     .  *(hc-b(i,j))/(b(i,j+1)-b(i,j))

        write(18,*)cd(id,1),cd(id,2),pv(i,j,1),pv(i,j+1,1),pv(i,j,2),
     .  pv(i,j+1,2)
      

 
       elseif(   (b(i,j)  .gt.hc).and. 
     .           (hc.ge.b(i,j+1)).and.
     .       (abs(b(i,j)  -hc).le.ee).and.
     .       (abs(b(i,j+1)-hc).le.ee)     )then
        iu=iu+1

        cu(iu,1)=pv(i,j,1)+(pv(i,j+1,1)-pv(i,j,1))
     .  *(hc-b(i,j))/(b(i,j+1)-b(i,j))

c         cu(iu,1)=pv(i,j,1)


        cu(iu,2)=pv(i,j,2)+(pv(i,j+1,2)-pv(i,j,2))
     .  *(hc-b(i,j))/(b(i,j+1)-b(i,j))

         write(19,*)cu(iu,1),cu(iu,2),pv(i,j,1),pv(i,j+1,1),
     .  pv(i,j,2),pv(i,j+1,2)

        else
        endif
 
       
       
101    continue
100    continue


       write(2,*)"cmaker.f"
       write(2,*)"id,iu,np,nz"
       write(2,*)id,iu,np,nz 
       write(2,*)"hc",hc
  


       return
       end
