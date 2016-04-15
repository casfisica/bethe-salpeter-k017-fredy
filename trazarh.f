       subroutine trazarh(print,w1,z1,qv,pcv,ssp,ssm,svp,svm,F,traza1)
       implicit none
       double precision qq,qv(4)    
       double precision  
     . u0,u1,u2,u3,u4,u5,f0,f1,f2,f3,f4,f5,
     . z1,u,w1,w2,uh(4)
       double complex pcv(4),pc2,pc_q,traza1,ssp,ssm,svp,svm,
     . pc_uh,q_uh,i_,F(8,0:5),
     . F10,F11,F12,F13,F14,F15,F20,F21,F22,F23,F24,F25,
     . F30,F31,F32,F33,F34,F35,F40,F41,F42,F43,F44,F45

       integer i,j
       logical print
       w2= 1-w1
     

       i_= (0d0,1d0)


        f0 =  1d0 
        f1 =  1d0
        f2 =  1d0
        f3 =  1d0
        f4 =  1d0
        f5 =  1d0

       F10 = F(1,0)
       F11 = F(1,1)
       F12 = F(1,2)
       F13 = F(1,3)
       F14 = F(1,4)
       F15 = F(1,5)
       F20 = F(2,0)
       F21 = F(2,1)
       F22 = F(2,2)
       F23 = F(2,3)
       F24 = F(2,4)
       F25 = F(2,5)
       F30 = F(3,0)
       F31 = F(3,1)
       F32 = F(3,2)
       F33 = F(3,3)
       F34 = F(3,4)
       F35 = F(3,5)
       F40 = F(4,0)
       F41 = F(4,1)
       F42 = F(4,2)
       F43 = F(4,3)
       F44 = F(4,4)
       F45 = F(4,5)

      

       u0 = u(0,z1)
       u1 = u(1,z1) 
       u2 = u(2,z1)
       u3 = u(3,z1)
       u4 = u(4,z1)
       u5 = u(5,z1) 

       uh(1)= 0d0
       uh(2)= 0d0 
       uh(3)= 0d0
       uh(4)= 1d0 

       
       
       pc_q = 0d0
       qq   = 0d0
       pc2  = 0d0
       pc_uh= 0d0
       q_uh = 0d0

       do 100 i =1,4
       qq =qv(i)*qv(i)+qq
       pc2=pcv(i)*pcv(i)+pc2
       pc_q  =pcv(i)*qv(i)+pc_q
       pc_uh =pcv(i)*uh(i)+pc_uh
       q_uh  =qv(i)*uh(i)  +q_uh
100    continue


          traza1 = + 4.D0*i_*ssp*ssm*F15*u5
     &  + 4.D0*i_*ssp*ssm*F14*u4
     &  + 4.D0*i_*ssp*ssm*F13*u3
     &  + 4.D0*i_*ssp*ssm*F12*u2
     &  + 4.D0*i_*ssp*ssm*F11*u1
     &  + 4.D0*i_*ssp*ssm*F10*u0
     &  + 4.D0*i_*qq*svp*svm*F15*u5
     &  + 4.D0*i_*qq*svp*svm*F14*u4
     &  + 4.D0*i_*qq*svp*svm*F13*u3
     &  + 4.D0*i_*qq*svp*svm*F12*u2
     &  + 4.D0*i_*qq*svp*svm*F11*u1
     &  + 4.D0*i_*qq*svp*svm*F10*u0
     &  - 4.D0*i_*PC2*svp*svm*F15*u5*w1*w2
     &  - 4.D0*i_*PC2*svp*svm*F14*u4*w1*w2
     &  - 4.D0*i_*PC2*svp*svm*F13*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*svp*svm*F12*u2*w1*w2
     &  - 4.D0*i_*PC2*svp*svm*F11*u1*w1*w2
     &  - 4.D0*i_*PC2*svp*svm*F10*u0*w1*w2
     &  + 4.D0*i_*PC2*ssm*svp*F25*u5*w1
     &  + 4.D0*i_*PC2*ssm*svp*F24*u4*w1
     &  + 4.D0*i_*PC2*ssm*svp*F23*u3*w1
     &  + 4.D0*i_*PC2*ssm*svp*F22*u2*w1
     &  + 4.D0*i_*PC2*ssm*svp*F21*u1*w1
     &  + 4.D0*i_*PC2*ssm*svp*F20*u0*w1
     &  + 4.D0*i_*PC2*ssp*svm*F25*u5*w2
     &  + 4.D0*i_*PC2*ssp*svm*F24*u4*w2
     &  + 4.D0*i_*PC2*ssp*svm*F23*u3*w2
     &  + 4.D0*i_*PC2*ssp*svm*F22*u2*w2
     &  + 4.D0*i_*PC2*ssp*svm*F21*u1*w2
     &  + 4.D0*i_*PC2*ssp*svm*F20*u0*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*svm*F45*u5*w2
     &  - 4.D0*i_*PC2*qq*svp*svm*F45*u5*w1
     &  - 4.D0*i_*PC2*qq*svp*svm*F44*u4*w2
     &  - 4.D0*i_*PC2*qq*svp*svm*F44*u4*w1
     &  - 4.D0*i_*PC2*qq*svp*svm*F43*u3*w2
     &  - 4.D0*i_*PC2*qq*svp*svm*F43*u3*w1
     &  - 4.D0*i_*PC2*qq*svp*svm*F42*u2*w2
     &  - 4.D0*i_*PC2*qq*svp*svm*F42*u2*w1
     &  - 4.D0*i_*PC2*qq*svp*svm*F41*u1*w2
     &  - 4.D0*i_*PC2*qq*svp*svm*F41*u1*w1
     &  - 4.D0*i_*PC2*qq*svp*svm*F40*u0*w2
     &  - 4.D0*i_*PC2*qq*svp*svm*F40*u0*w1
     &  - 4.D0*PC_q*i_*svp*svm*F15*u5*w2
     &  + 4.D0*PC_q*i_*svp*svm*F15*u5*w1
     &  - 4.D0*PC_q*i_*svp*svm*F14*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svp*svm*F14*u4*w1
     &  - 4.D0*PC_q*i_*svp*svm*F13*u3*w2
     &  + 4.D0*PC_q*i_*svp*svm*F13*u3*w1
     &  - 4.D0*PC_q*i_*svp*svm*F12*u2*w2
     &  + 4.D0*PC_q*i_*svp*svm*F12*u2*w1
     &  - 4.D0*PC_q*i_*svp*svm*F11*u1*w2
     &  + 4.D0*PC_q*i_*svp*svm*F11*u1*w1
     &  - 4.D0*PC_q*i_*svp*svm*F10*u0*w2
     &  + 4.D0*PC_q*i_*svp*svm*F10*u0*w1
     &  + 4.D0*PC_q*i_*ssm*svp*F25*u5
     &  + 4.D0*PC_q*i_*ssm*svp*F24*u4
     &  + 4.D0*PC_q*i_*ssm*svp*F23*u3
     &  + 4.D0*PC_q*i_*ssm*svp*F22*u2
     &  + 4.D0*PC_q*i_*ssm*svp*F21*u1
     &  + 4.D0*PC_q*i_*ssm*svp*F20*u0
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssp*svm*F25*u5
     &  - 4.D0*PC_q*i_*ssp*svm*F24*u4
     &  - 4.D0*PC_q*i_*ssp*svm*F23*u3
     &  - 4.D0*PC_q*i_*ssp*svm*F22*u2
     &  - 4.D0*PC_q*i_*ssp*svm*F21*u1
     &  - 4.D0*PC_q*i_*ssp*svm*F20*u0
     &  + 4.D0*PC_q*i_*qq*ssm*svp*F35*u5
     &  + 4.D0*PC_q*i_*qq*ssm*svp*F34*u4
     &  + 4.D0*PC_q*i_*qq*ssm*svp*F33*u3
     &  + 4.D0*PC_q*i_*qq*ssm*svp*F32*u2
     &  + 4.D0*PC_q*i_*qq*ssm*svp*F31*u1
     &  + 4.D0*PC_q*i_*qq*ssm*svp*F30*u0
     &  - 4.D0*PC_q*i_*qq*ssp*svm*F35*u5
     &  - 4.D0*PC_q*i_*qq*ssp*svm*F34*u4
     &  - 4.D0*PC_q*i_*qq*ssp*svm*F33*u3
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssp*svm*F32*u2
     &  - 4.D0*PC_q*i_*qq*ssp*svm*F31*u1
     &  - 4.D0*PC_q*i_*qq*ssp*svm*F30*u0
     &  + 4.D0*PC_q**2*i_*svp*svm*F45*u5*w2
     &  + 4.D0*PC_q**2*i_*svp*svm*F45*u5*w1
     &  + 4.D0*PC_q**2*i_*svp*svm*F44*u4*w2
     &  + 4.D0*PC_q**2*i_*svp*svm*F44*u4*w1
     &  + 4.D0*PC_q**2*i_*svp*svm*F43*u3*w2
     &  + 4.D0*PC_q**2*i_*svp*svm*F43*u3*w1
     &  + 4.D0*PC_q**2*i_*svp*svm*F42*u2*w2
     &  + 4.D0*PC_q**2*i_*svp*svm*F42*u2*w1
     &  + 4.D0*PC_q**2*i_*svp*svm*F41*u1*w2
     &  + 4.D0*PC_q**2*i_*svp*svm*F41*u1*w1
     &  + 4.D0*PC_q**2*i_*svp*svm*F40*u0*w2
     &  + 4.D0*PC_q**2*i_*svp*svm*F40*u0*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*svp*F35*u5*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F34*u4*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F33*u3*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F32*u2*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F31*u1*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F30*u0*w1
     &  + 4.D0*PC_q**2*i_*ssp*svm*F35*u5*w2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F34*u4*w2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F33*u3*w2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F32*u2*w2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F31*u1*w2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F30*u0*w2
     &






       if(print.eqv..true.)then
       write(2,*)"trazaD-variables___________________________________"
       write(2,*)"qq",qq
       write(2,*)"F"  
       do  103 i =1,4
       write(2,*)(F(i,j),j=0,4)
103    continue
       write(2,*)"qv(1),qv(2),qv(3),qv(4)",qv(1),qv(2),qv(3),qv(4)
       write(2,*)"F10,F11,..."
       Write(2,*)      
     . F10,F11,F12,F13,F14,F15,F20,F21,F22,F23,F24,F25,
     . F30,F31,F32,F33,F34,F35,F40,F41,F42,F43,F44,F45
       write(2,*)"u0,u1,u2,u3,u4,u5",u0,u1,u2,u3,u4,u5
       write(2,*)"f0,f1,f2,f3,f4,f5",f0,f1,f2,f3,f4,f5
       write(2,*)"z1,w1,w2",z1,w1,w2
       write(2,*)"uh(1),uh(2),uh(3),uh(4)",uh(1),uh(2),uh(3),uh(4)
       write(2,*)"pcv(i)"
       write(2,*)"pcv(1),pcv(2),pcv(3),pcv(4)",
     . pcv(1),pcv(2),pcv(3),pcv(4)
       write(2,*)"pc2,pc_q,traza1",pc2,pc_q,traza1
       write(2,*)"ssp,ssm,svp,svm",ssp,ssm,svp,svm
       write(2,*)"pc_uh,q_uh",pc_uh,q_uh
       print = .false. 
       else
       endif



       return
       end
