       subroutine trazad(print,w1,z1,qv,pcv,ssp,ssm,svp,svm,F,traza1)
       implicit none
       double precision qq,qv(4)    
       double precision  
     . u0,u1,u2,u3,u4,u5,f0,f1,f2,f3,f4,f5,
     . z1,u,w1,w2,uh(4)
       double complex pcv(4),pc2,pc_q,traza1,ssp,ssm,svp,svm,
     . pc_uh,q_uh,F(4,0:5),
     . F10,F11,F12,F13,F14,F15,F20,F21,F22,F23,F24,F25,
     . F30,F31,F32,F33,F34,F35,F40,F41,F42,F43,F44,F45

       integer i,j
       logical print
       
        w2=1d0-w1



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





      traza1 = + 4.D0*PC_q*PC_uh*ssm*svp*F45*u5*w1
     &  + 4.D0*PC_q*PC_uh*ssm*svp*F44*u4*w1
     &  + 4.D0*PC_q*PC_uh*ssm*svp*F43*u3*w1
     &  + 4.D0*PC_q*PC_uh*ssm*svp*F42*u2*w1
     &  + 4.D0*PC_q*PC_uh*ssm*svp*F41*u1*w1
     &  + 4.D0*PC_q*PC_uh*ssm*svp*F40*u0*w1
     &  - 4.D0*PC_q*PC_uh*ssp*svm*F45*u5*w2
     &  - 4.D0*PC_q*PC_uh*ssp*svm*F44*u4*w2
     &  - 4.D0*PC_q*PC_uh*ssp*svm*F43*u3*w2
     &  - 4.D0*PC_q*PC_uh*ssp*svm*F42*u2*w2
     &  - 4.D0*PC_q*PC_uh*ssp*svm*F41*u1*w2
     &  - 4.D0*PC_q*PC_uh*ssp*svm*F40*u0*w2
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*u1*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*u1*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*u0*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*u0*w1
     &  - 8.D0*PC_q*q_uh*svp*svm*F25*u5
     &  - 8.D0*PC_q*q_uh*svp*svm*F24*u4
     &  - 8.D0*PC_q*q_uh*svp*svm*F23*u3
     &  - 8.D0*PC_q*q_uh*svp*svm*F22*u2
     &  - 8.D0*PC_q*q_uh*svp*svm*F21*u1
     &  - 8.D0*PC_q*q_uh*svp*svm*F20*u0
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*ssm*svp*F45*u5
     &  - 4.D0*PC_q*q_uh*ssm*svp*F44*u4
     &  - 4.D0*PC_q*q_uh*ssm*svp*F43*u3
     &  - 4.D0*PC_q*q_uh*ssm*svp*F42*u2
     &  - 4.D0*PC_q*q_uh*ssm*svp*F41*u1
     &  - 4.D0*PC_q*q_uh*ssm*svp*F40*u0
     &  - 4.D0*PC_q*q_uh*ssp*svm*F45*u5
     &  - 4.D0*PC_q*q_uh*ssp*svm*F44*u4
     &  - 4.D0*PC_q*q_uh*ssp*svm*F43*u3
     &  - 4.D0*PC_q*q_uh*ssp*svm*F42*u2
     &  - 4.D0*PC_q*q_uh*ssp*svm*F41*u1
     &  - 4.D0*PC_q*q_uh*ssp*svm*F40*u0
     &  - 4.D0*PC_q*q_uh*ssp*ssm*F35*u5
     &  - 4.D0*PC_q*q_uh*ssp*ssm*F34*u4
     &  - 4.D0*PC_q*q_uh*ssp*ssm*F33*u3
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*ssp*ssm*F32*u2
     &  - 4.D0*PC_q*q_uh*ssp*ssm*F31*u1
     &  - 4.D0*PC_q*q_uh*ssp*ssm*F30*u0
     &  - 4.D0*PC_q*q_uh*qq*svp*svm*F35*u5
     &  - 4.D0*PC_q*q_uh*qq*svp*svm*F34*u4
     &  - 4.D0*PC_q*q_uh*qq*svp*svm*F33*u3
     &  - 4.D0*PC_q*q_uh*qq*svp*svm*F32*u2
     &  - 4.D0*PC_q*q_uh*qq*svp*svm*F31*u1
     &  - 4.D0*PC_q*q_uh*qq*svp*svm*F30*u0
     &  - 4.D0*PC_q*q_uh*PC2*svp*svm*F35*u5*w1*w2
     &  - 4.D0*PC_q*q_uh*PC2*svp*svm*F34*u4*w1*w2
     &  - 4.D0*PC_q*q_uh*PC2*svp*svm*F33*u3*w1*w2
     &  - 4.D0*PC_q*q_uh*PC2*svp*svm*F32*u2*w1*w2
     &  - 4.D0*PC_q*q_uh*PC2*svp*svm*F31*u1*w1*w2
     &  - 4.D0*PC_q*q_uh*PC2*svp*svm*F30*u0*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q**2*PC_uh*svp*svm*F35*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*svp*svm*F34*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*svp*svm*F33*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*svp*svm*F32*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*svp*svm*F31*u1*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*svp*svm*F30*u0*w1*w2
     &  + 4.D0*PC_uh*ssm*svp*F15*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F14*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F12*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*u1*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*u0*w1
     &  + 4.D0*PC_uh*ssp*svm*F15*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F14*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_uh*ssp*svm*F12*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F11*u1*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*u0*w2
     &  - 4.D0*PC_uh*ssp*ssm*F25*u5
     &  - 4.D0*PC_uh*ssp*ssm*F24*u4
     &  - 4.D0*PC_uh*ssp*ssm*F23*u3
     &  - 4.D0*PC_uh*ssp*ssm*F22*u2
     &  - 4.D0*PC_uh*ssp*ssm*F21*u1
     &  - 4.D0*PC_uh*ssp*ssm*F20*u0
     &  + 4.D0*PC_uh*qq*svp*svm*F25*u5
     &  + 4.D0*PC_uh*qq*svp*svm*F24*u4
     &  + 4.D0*PC_uh*qq*svp*svm*F23*u3
     &  + 4.D0*PC_uh*qq*svp*svm*F22*u2
     &  + 4.D0*PC_uh*qq*svp*svm*F21*u1
     &  + 4.D0*PC_uh*qq*svp*svm*F20*u0
     &
      traza1 = traza1 + 4.D0*PC_uh*qq*ssm*svp*F45*u5
     &  + 4.D0*PC_uh*qq*ssm*svp*F44*u4
     &  + 4.D0*PC_uh*qq*ssm*svp*F43*u3
     &  + 4.D0*PC_uh*qq*ssm*svp*F42*u2
     &  + 4.D0*PC_uh*qq*ssm*svp*F41*u1
     &  + 4.D0*PC_uh*qq*ssm*svp*F40*u0
     &  + 4.D0*PC_uh*qq*ssp*svm*F45*u5
     &  + 4.D0*PC_uh*qq*ssp*svm*F44*u4
     &  + 4.D0*PC_uh*qq*ssp*svm*F43*u3
     &  + 4.D0*PC_uh*qq*ssp*svm*F42*u2
     &  + 4.D0*PC_uh*qq*ssp*svm*F41*u1
     &  + 4.D0*PC_uh*qq*ssp*svm*F40*u0
     &  + 4.D0*PC_uh*PC2*svp*svm*F25*u5*w1*w2
     &  + 4.D0*PC_uh*PC2*svp*svm*F24*u4*w1*w2
     &  + 4.D0*PC_uh*PC2*svp*svm*F23*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_uh*PC2*svp*svm*F22*u2*w1*w2
     &  + 4.D0*PC_uh*PC2*svp*svm*F21*u1*w1*w2
     &  + 4.D0*PC_uh*PC2*svp*svm*F20*u0*w1*w2
     &  + 4.D0*q_uh*ssm*svp*F15*u5
     &  + 4.D0*q_uh*ssm*svp*F14*u4
     &  + 4.D0*q_uh*ssm*svp*F13*u3
     &  + 4.D0*q_uh*ssm*svp*F12*u2
     &  + 4.D0*q_uh*ssm*svp*F11*u1
     &  + 4.D0*q_uh*ssm*svp*F10*u0
     &  - 4.D0*q_uh*ssp*svm*F15*u5
     &  - 4.D0*q_uh*ssp*svm*F14*u4
     &  - 4.D0*q_uh*ssp*svm*F13*u3
     &  - 4.D0*q_uh*ssp*svm*F12*u2
     &  - 4.D0*q_uh*ssp*svm*F11*u1
     &  - 4.D0*q_uh*ssp*svm*F10*u0
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*svp*svm*F25*u5*w2
     &  - 4.D0*q_uh*PC2*svp*svm*F25*u5*w1
     &  + 4.D0*q_uh*PC2*svp*svm*F24*u4*w2
     &  - 4.D0*q_uh*PC2*svp*svm*F24*u4*w1
     &  + 4.D0*q_uh*PC2*svp*svm*F23*u3*w2
     &  - 4.D0*q_uh*PC2*svp*svm*F23*u3*w1
     &  + 4.D0*q_uh*PC2*svp*svm*F22*u2*w2
     &  - 4.D0*q_uh*PC2*svp*svm*F22*u2*w1
     &  + 4.D0*q_uh*PC2*svp*svm*F21*u1*w2
     &  - 4.D0*q_uh*PC2*svp*svm*F21*u1*w1
     &  + 4.D0*q_uh*PC2*svp*svm*F20*u0*w2
     &  - 4.D0*q_uh*PC2*svp*svm*F20*u0*w1
     &  - 4.D0*q_uh*PC2*ssm*svp*F45*u5*w1
     &  - 4.D0*q_uh*PC2*ssm*svp*F44*u4*w1
     &  - 4.D0*q_uh*PC2*ssm*svp*F43*u3*w1
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*ssm*svp*F42*u2*w1
     &  - 4.D0*q_uh*PC2*ssm*svp*F41*u1*w1
     &  - 4.D0*q_uh*PC2*ssm*svp*F40*u0*w1
     &  + 4.D0*q_uh*PC2*ssp*svm*F45*u5*w2
     &  + 4.D0*q_uh*PC2*ssp*svm*F44*u4*w2
     &  + 4.D0*q_uh*PC2*ssp*svm*F43*u3*w2
     &  + 4.D0*q_uh*PC2*ssp*svm*F42*u2*w2
     &  + 4.D0*q_uh*PC2*ssp*svm*F41*u1*w2
     &  + 4.D0*q_uh*PC2*ssp*svm*F40*u0*w2
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
