       subroutine traza2(print,w1,z1,qv,pcv,ssp,ssm,svp,
     . svm,dssp,dssm,dsvp,dsvm,F,traza1)
       implicit none
       double precision qq,qv(4),uh(4)    
       double precision 
     . u0,u1,u2,u3,u4,u5,c0,c1,c2,c3,c4,c5,f1,f2,f3,f4,
     . z1,w1,w2,u
       double complex pcv(4),pc2,pc_q,traza1,pc_uh,i_,q_uh,
     . ssp,ssm,svp,svm,dssp,dssm,dsvp,dsvm,F(4,0:5),
     . F10,F11,F12,F13,F14,F15,F20,F21,F22,F23,F24,F25,
     . F30,F31,F32,F33,F34,F35,F40,F41,F42,F43,F44,F45,
     . F10r,F11r,F12r,F13r,F14r,F15r,F20r,F21r,F22r,F23r,F24r,F25r,
     . F30r,F31r,F32r,F33r,F34r,F35r,F40r,F41r,F42r,F43r,F44r,F45r,
     . F10i,F11i,F12i,F13i,F14i,F15i,F20i,F21i,F22i,F23i,F24i,F25i,
     . F30i,F31i,F32i,F33i,F34i,F35i,F40i,F41i,F42i,F43i,F44i,F45i

 
       integer i,j
       logical print
        i_ = (0d0,1d0)

        w2=1d0-w1


        f1 =  1d0
        f2 =  1d0
        f3 =  1d0
        f4 =  1d0
        
        c0 =  1d0
        c1 =  1d0
        c2 =  1d0
        c3 =  1d0
        c4 =  1d0
        c5 =  1d0

 

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

       F10r = dble(F(1,0))
       F11r = dble(F(1,1))
       F12r = dble(F(1,2))
       F13r = dble(F(1,3))
       F14r = dble(F(1,4))
       F15r = dble(F(1,5))
       F20r = dble(F(2,0))
       F21r = dble(F(2,1))
       F22r = dble(F(2,2))
       F23r = dble(F(2,3))
       F24r = dble(F(2,4))
       F25r = dble(F(2,5))
       F30r = dble(F(3,0))
       F31r = dble(F(3,1))
       F32r = dble(F(3,2))
       F33r = dble(F(3,3))
       F34r = dble(F(3,4))
       F35r = dble(F(3,5))
       F40r = dble(F(4,0))
       F41r = dble(F(4,1))
       F42r = dble(F(4,2))
       F43r = dble(F(4,3))
       F44r = dble(F(4,4))
       F45r = dble(F(4,5))

       F10i = dimag(F(1,0))
       F11i = dimag(F(1,1))
       F12i = dimag(F(1,2))
       F13i = dimag(F(1,3))
       F14i = dimag(F(1,4))
       F15i = dimag(F(1,5))
       F20i = dimag(F(2,0))
       F21i = dimag(F(2,1))
       F22i = dimag(F(2,2))
       F23i = dimag(F(2,3))
       F24i = dimag(F(2,4))
       F25i = dimag(F(2,5))
       F30i = dimag(F(3,0))
       F31i = dimag(F(3,1))
       F32i = dimag(F(3,2))
       F33i = dimag(F(3,3))
       F34i = dimag(F(3,4))
       F35i = dimag(F(3,5))
       F40i = dimag(F(4,0))
       F41i = dimag(F(4,1))
       F42i = dimag(F(4,2))
       F43i = dimag(F(4,3))
       F44i = dimag(F(4,4))
       F45i = dimag(F(4,5))

      

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
       pc_q =pcv(i)*qv(i)+pc_q
       pc_uh=pcv(i)*uh(i)+pc_uh
       q_uh =  qv(i)*uh(i)+q_uh
100    continue

            traza1 = - 4.D0*ssp*ssm*F15*F15r*u5**2*c5*f1
     &  - 4.D0*ssp*ssm*F15*F14r*u4*u5*c4*f1
     &  - 4.D0*ssp*ssm*F15*F13r*u3*u5*c3*f1
     &  - 4.D0*ssp*ssm*F15*F12r*u2*u5*c2*f1
     &  - 4.D0*ssp*ssm*F15*F11r*u1*u5*c1*f1
     &  - 4.D0*ssp*ssm*F15*F10r*u0*u5*c0*f1
     &  - 4.D0*ssp*ssm*F14*F15r*u4*u5*c5*f1
     &  - 4.D0*ssp*ssm*F14*F14r*u4**2*c4*f1
     &  - 4.D0*ssp*ssm*F14*F13r*u3*u4*c3*f1
     &  - 4.D0*ssp*ssm*F14*F12r*u2*u4*c2*f1
     &  - 4.D0*ssp*ssm*F14*F11r*u1*u4*c1*f1
     &  - 4.D0*ssp*ssm*F14*F10r*u0*u4*c0*f1
     &  - 4.D0*ssp*ssm*F13*F15r*u3*u5*c5*f1
     &  - 4.D0*ssp*ssm*F13*F14r*u3*u4*c4*f1
     &  - 4.D0*ssp*ssm*F13*F13r*u3**2*c3*f1
     &
      traza1 = traza1 - 4.D0*ssp*ssm*F13*F12r*u2*u3*c2*f1
     &  - 4.D0*ssp*ssm*F13*F11r*u1*u3*c1*f1
     &  - 4.D0*ssp*ssm*F13*F10r*u0*u3*c0*f1
     &  - 4.D0*ssp*ssm*F12*F15r*u2*u5*c5*f1
     &  - 4.D0*ssp*ssm*F12*F14r*u2*u4*c4*f1
     &  - 4.D0*ssp*ssm*F12*F13r*u2*u3*c3*f1
     &  - 4.D0*ssp*ssm*F12*F12r*u2**2*c2*f1
     &  - 4.D0*ssp*ssm*F12*F11r*u1*u2*c1*f1
     &  - 4.D0*ssp*ssm*F12*F10r*u0*u2*c0*f1
     &  - 4.D0*ssp*ssm*F11*F15r*u1*u5*c5*f1
     &  - 4.D0*ssp*ssm*F11*F14r*u1*u4*c4*f1
     &  - 4.D0*ssp*ssm*F11*F13r*u1*u3*c3*f1
     &  - 4.D0*ssp*ssm*F11*F12r*u1*u2*c2*f1
     &  - 4.D0*ssp*ssm*F11*F11r*u1**2*c1*f1
     &  - 4.D0*ssp*ssm*F11*F10r*u0*u1*c0*f1
     &
      traza1 = traza1 - 4.D0*ssp*ssm*F10*F15r*u0*u5*c5*f1
     &  - 4.D0*ssp*ssm*F10*F14r*u0*u4*c4*f1
     &  - 4.D0*ssp*ssm*F10*F13r*u0*u3*c3*f1
     &  - 4.D0*ssp*ssm*F10*F12r*u0*u2*c2*f1
     &  - 4.D0*ssp*ssm*F10*F11r*u0*u1*c1*f1
     &  - 4.D0*ssp*ssm*F10*F10r*u0**2*c0*f1
     &  - 4.D0*qq*svp*svm*F15*F15r*u5**2*c5*f1
     &  - 4.D0*qq*svp*svm*F15*F14r*u4*u5*c4*f1
     &  - 4.D0*qq*svp*svm*F15*F13r*u3*u5*c3*f1
     &  - 4.D0*qq*svp*svm*F15*F12r*u2*u5*c2*f1
     &  - 4.D0*qq*svp*svm*F15*F11r*u1*u5*c1*f1
     &  - 4.D0*qq*svp*svm*F15*F10r*u0*u5*c0*f1
     &  - 4.D0*qq*svp*svm*F14*F15r*u4*u5*c5*f1
     &  - 4.D0*qq*svp*svm*F14*F14r*u4**2*c4*f1
     &  - 4.D0*qq*svp*svm*F14*F13r*u3*u4*c3*f1
     &
      traza1 = traza1 - 4.D0*qq*svp*svm*F14*F12r*u2*u4*c2*f1
     &  - 4.D0*qq*svp*svm*F14*F11r*u1*u4*c1*f1
     &  - 4.D0*qq*svp*svm*F14*F10r*u0*u4*c0*f1
     &  - 4.D0*qq*svp*svm*F13*F15r*u3*u5*c5*f1
     &  - 4.D0*qq*svp*svm*F13*F14r*u3*u4*c4*f1
     &  - 4.D0*qq*svp*svm*F13*F13r*u3**2*c3*f1
     &  - 4.D0*qq*svp*svm*F13*F12r*u2*u3*c2*f1
     &  - 4.D0*qq*svp*svm*F13*F11r*u1*u3*c1*f1
     &  - 4.D0*qq*svp*svm*F13*F10r*u0*u3*c0*f1
     &  - 4.D0*qq*svp*svm*F12*F15r*u2*u5*c5*f1
     &  - 4.D0*qq*svp*svm*F12*F14r*u2*u4*c4*f1
     &  - 4.D0*qq*svp*svm*F12*F13r*u2*u3*c3*f1
     &  - 4.D0*qq*svp*svm*F12*F12r*u2**2*c2*f1
     &  - 4.D0*qq*svp*svm*F12*F11r*u1*u2*c1*f1
     &  - 4.D0*qq*svp*svm*F12*F10r*u0*u2*c0*f1
     &
      traza1 = traza1 - 4.D0*qq*svp*svm*F11*F15r*u1*u5*c5*f1
     &  - 4.D0*qq*svp*svm*F11*F14r*u1*u4*c4*f1
     &  - 4.D0*qq*svp*svm*F11*F13r*u1*u3*c3*f1
     &  - 4.D0*qq*svp*svm*F11*F12r*u1*u2*c2*f1
     &  - 4.D0*qq*svp*svm*F11*F11r*u1**2*c1*f1
     &  - 4.D0*qq*svp*svm*F11*F10r*u0*u1*c0*f1
     &  - 4.D0*qq*svp*svm*F10*F15r*u0*u5*c5*f1
     &  - 4.D0*qq*svp*svm*F10*F14r*u0*u4*c4*f1
     &  - 4.D0*qq*svp*svm*F10*F13r*u0*u3*c3*f1
     &  - 4.D0*qq*svp*svm*F10*F12r*u0*u2*c2*f1
     &  - 4.D0*qq*svp*svm*F10*F11r*u0*u1*c1*f1
     &  - 4.D0*qq*svp*svm*F10*F10r*u0**2*c0*f1
     &  + 4.D0*PC2*svp*svm*F15*F15r*u5**2*c5*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F15*F14r*u4*u5*c4*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F15*F13r*u3*u5*c3*f1*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svp*svm*F15*F12r*u2*u5*c2*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F15*F11r*u1*u5*c1*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F15*F10r*u0*u5*c0*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F14*F15r*u4*u5*c5*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F14*F14r*u4**2*c4*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F14*F13r*u3*u4*c3*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F14*F12r*u2*u4*c2*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F14*F11r*u1*u4*c1*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F14*F10r*u0*u4*c0*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F13*F15r*u3*u5*c5*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F13*F14r*u3*u4*c4*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F13*F13r*u3**2*c3*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F13*F12r*u2*u3*c2*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F13*F11r*u1*u3*c1*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F13*F10r*u0*u3*c0*f1*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svp*svm*F12*F15r*u2*u5*c5*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F12*F14r*u2*u4*c4*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F12*F13r*u2*u3*c3*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F12*F12r*u2**2*c2*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F12*F11r*u1*u2*c1*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F12*F10r*u0*u2*c0*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F11*F15r*u1*u5*c5*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F11*F14r*u1*u4*c4*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F11*F13r*u1*u3*c3*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F11*F12r*u1*u2*c2*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F11*F11r*u1**2*c1*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F11*F10r*u0*u1*c0*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F10*F15r*u0*u5*c5*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F10*F14r*u0*u4*c4*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F10*F13r*u0*u3*c3*f1*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svp*svm*F10*F12r*u0*u2*c2*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F10*F11r*u0*u1*c1*f1*w1*w2
     &  + 4.D0*PC2*svp*svm*F10*F10r*u0**2*c0*f1*w1*w2
     &  - 4.D0*PC2*ssm*svp*F25*F15r*u5**2*c5*f1*w1
     &  - 4.D0*PC2*ssm*svp*F25*F14r*u4*u5*c4*f1*w1
     &  - 4.D0*PC2*ssm*svp*F25*F13r*u3*u5*c3*f1*w1
     &  - 4.D0*PC2*ssm*svp*F25*F12r*u2*u5*c2*f1*w1
     &  - 4.D0*PC2*ssm*svp*F25*F11r*u1*u5*c1*f1*w1
     &  - 4.D0*PC2*ssm*svp*F25*F10r*u0*u5*c0*f1*w1
     &  - 4.D0*PC2*ssm*svp*F24*F15r*u4*u5*c5*f1*w1
     &  - 4.D0*PC2*ssm*svp*F24*F14r*u4**2*c4*f1*w1
     &  - 4.D0*PC2*ssm*svp*F24*F13r*u3*u4*c3*f1*w1
     &  - 4.D0*PC2*ssm*svp*F24*F12r*u2*u4*c2*f1*w1
     &  - 4.D0*PC2*ssm*svp*F24*F11r*u1*u4*c1*f1*w1
     &  - 4.D0*PC2*ssm*svp*F24*F10r*u0*u4*c0*f1*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*svp*F23*F15r*u3*u5*c5*f1*w1
     &  - 4.D0*PC2*ssm*svp*F23*F14r*u3*u4*c4*f1*w1
     &  - 4.D0*PC2*ssm*svp*F23*F13r*u3**2*c3*f1*w1
     &  - 4.D0*PC2*ssm*svp*F23*F12r*u2*u3*c2*f1*w1
     &  - 4.D0*PC2*ssm*svp*F23*F11r*u1*u3*c1*f1*w1
     &  - 4.D0*PC2*ssm*svp*F23*F10r*u0*u3*c0*f1*w1
     &  - 4.D0*PC2*ssm*svp*F22*F15r*u2*u5*c5*f1*w1
     &  - 4.D0*PC2*ssm*svp*F22*F14r*u2*u4*c4*f1*w1
     &  - 4.D0*PC2*ssm*svp*F22*F13r*u2*u3*c3*f1*w1
     &  - 4.D0*PC2*ssm*svp*F22*F12r*u2**2*c2*f1*w1
     &  - 4.D0*PC2*ssm*svp*F22*F11r*u1*u2*c1*f1*w1
     &  - 4.D0*PC2*ssm*svp*F22*F10r*u0*u2*c0*f1*w1
     &  - 4.D0*PC2*ssm*svp*F21*F15r*u1*u5*c5*f1*w1
     &  - 4.D0*PC2*ssm*svp*F21*F14r*u1*u4*c4*f1*w1
     &  - 4.D0*PC2*ssm*svp*F21*F13r*u1*u3*c3*f1*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*svp*F21*F12r*u1*u2*c2*f1*w1
     &  - 4.D0*PC2*ssm*svp*F21*F11r*u1**2*c1*f1*w1
     &  - 4.D0*PC2*ssm*svp*F21*F10r*u0*u1*c0*f1*w1
     &  - 4.D0*PC2*ssm*svp*F20*F15r*u0*u5*c5*f1*w1
     &  - 4.D0*PC2*ssm*svp*F20*F14r*u0*u4*c4*f1*w1
     &  - 4.D0*PC2*ssm*svp*F20*F13r*u0*u3*c3*f1*w1
     &  - 4.D0*PC2*ssm*svp*F20*F12r*u0*u2*c2*f1*w1
     &  - 4.D0*PC2*ssm*svp*F20*F11r*u0*u1*c1*f1*w1
     &  - 4.D0*PC2*ssm*svp*F20*F10r*u0**2*c0*f1*w1
     &  - 4.D0*PC2*ssm*svp*F15*F25r*u5**2*c5*f2*w1
     &  - 4.D0*PC2*ssm*svp*F15*F24r*u4*u5*c4*f2*w1
     &  - 4.D0*PC2*ssm*svp*F15*F23r*u3*u5*c3*f2*w1
     &  - 4.D0*PC2*ssm*svp*F15*F22r*u2*u5*c2*f2*w1
     &  - 4.D0*PC2*ssm*svp*F15*F21r*u1*u5*c1*f2*w1
     &  - 4.D0*PC2*ssm*svp*F15*F20r*u0*u5*c0*f2*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*svp*F14*F25r*u4*u5*c5*f2*w1
     &  - 4.D0*PC2*ssm*svp*F14*F24r*u4**2*c4*f2*w1
     &  - 4.D0*PC2*ssm*svp*F14*F23r*u3*u4*c3*f2*w1
     &  - 4.D0*PC2*ssm*svp*F14*F22r*u2*u4*c2*f2*w1
     &  - 4.D0*PC2*ssm*svp*F14*F21r*u1*u4*c1*f2*w1
     &  - 4.D0*PC2*ssm*svp*F14*F20r*u0*u4*c0*f2*w1
     &  - 4.D0*PC2*ssm*svp*F13*F25r*u3*u5*c5*f2*w1
     &  - 4.D0*PC2*ssm*svp*F13*F24r*u3*u4*c4*f2*w1
     &  - 4.D0*PC2*ssm*svp*F13*F23r*u3**2*c3*f2*w1
     &  - 4.D0*PC2*ssm*svp*F13*F22r*u2*u3*c2*f2*w1
     &  - 4.D0*PC2*ssm*svp*F13*F21r*u1*u3*c1*f2*w1
     &  - 4.D0*PC2*ssm*svp*F13*F20r*u0*u3*c0*f2*w1
     &  - 4.D0*PC2*ssm*svp*F12*F25r*u2*u5*c5*f2*w1
     &  - 4.D0*PC2*ssm*svp*F12*F24r*u2*u4*c4*f2*w1
     &  - 4.D0*PC2*ssm*svp*F12*F23r*u2*u3*c3*f2*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*svp*F12*F22r*u2**2*c2*f2*w1
     &  - 4.D0*PC2*ssm*svp*F12*F21r*u1*u2*c1*f2*w1
     &  - 4.D0*PC2*ssm*svp*F12*F20r*u0*u2*c0*f2*w1
     &  - 4.D0*PC2*ssm*svp*F11*F25r*u1*u5*c5*f2*w1
     &  - 4.D0*PC2*ssm*svp*F11*F24r*u1*u4*c4*f2*w1
     &  - 4.D0*PC2*ssm*svp*F11*F23r*u1*u3*c3*f2*w1
     &  - 4.D0*PC2*ssm*svp*F11*F22r*u1*u2*c2*f2*w1
     &  - 4.D0*PC2*ssm*svp*F11*F21r*u1**2*c1*f2*w1
     &  - 4.D0*PC2*ssm*svp*F11*F20r*u0*u1*c0*f2*w1
     &  - 4.D0*PC2*ssm*svp*F10*F25r*u0*u5*c5*f2*w1
     &  - 4.D0*PC2*ssm*svp*F10*F24r*u0*u4*c4*f2*w1
     &  - 4.D0*PC2*ssm*svp*F10*F23r*u0*u3*c3*f2*w1
     &  - 4.D0*PC2*ssm*svp*F10*F22r*u0*u2*c2*f2*w1
     &  - 4.D0*PC2*ssm*svp*F10*F21r*u0*u1*c1*f2*w1
     &  - 4.D0*PC2*ssm*svp*F10*F20r*u0**2*c0*f2*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssp*svm*F25*F15r*u5**2*c5*f1*w2
     &  - 4.D0*PC2*ssp*svm*F25*F14r*u4*u5*c4*f1*w2
     &  - 4.D0*PC2*ssp*svm*F25*F13r*u3*u5*c3*f1*w2
     &  - 4.D0*PC2*ssp*svm*F25*F12r*u2*u5*c2*f1*w2
     &  - 4.D0*PC2*ssp*svm*F25*F11r*u1*u5*c1*f1*w2
     &  - 4.D0*PC2*ssp*svm*F25*F10r*u0*u5*c0*f1*w2
     &  - 4.D0*PC2*ssp*svm*F24*F15r*u4*u5*c5*f1*w2
     &  - 4.D0*PC2*ssp*svm*F24*F14r*u4**2*c4*f1*w2
     &  - 4.D0*PC2*ssp*svm*F24*F13r*u3*u4*c3*f1*w2
     &  - 4.D0*PC2*ssp*svm*F24*F12r*u2*u4*c2*f1*w2
     &  - 4.D0*PC2*ssp*svm*F24*F11r*u1*u4*c1*f1*w2
     &  - 4.D0*PC2*ssp*svm*F24*F10r*u0*u4*c0*f1*w2
     &  - 4.D0*PC2*ssp*svm*F23*F15r*u3*u5*c5*f1*w2
     &  - 4.D0*PC2*ssp*svm*F23*F14r*u3*u4*c4*f1*w2
     &  - 4.D0*PC2*ssp*svm*F23*F13r*u3**2*c3*f1*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*svm*F23*F12r*u2*u3*c2*f1*w2
     &  - 4.D0*PC2*ssp*svm*F23*F11r*u1*u3*c1*f1*w2
     &  - 4.D0*PC2*ssp*svm*F23*F10r*u0*u3*c0*f1*w2
     &  - 4.D0*PC2*ssp*svm*F22*F15r*u2*u5*c5*f1*w2
     &  - 4.D0*PC2*ssp*svm*F22*F14r*u2*u4*c4*f1*w2
     &  - 4.D0*PC2*ssp*svm*F22*F13r*u2*u3*c3*f1*w2
     &  - 4.D0*PC2*ssp*svm*F22*F12r*u2**2*c2*f1*w2
     &  - 4.D0*PC2*ssp*svm*F22*F11r*u1*u2*c1*f1*w2
     &  - 4.D0*PC2*ssp*svm*F22*F10r*u0*u2*c0*f1*w2
     &  - 4.D0*PC2*ssp*svm*F21*F15r*u1*u5*c5*f1*w2
     &  - 4.D0*PC2*ssp*svm*F21*F14r*u1*u4*c4*f1*w2
     &  - 4.D0*PC2*ssp*svm*F21*F13r*u1*u3*c3*f1*w2
     &  - 4.D0*PC2*ssp*svm*F21*F12r*u1*u2*c2*f1*w2
     &  - 4.D0*PC2*ssp*svm*F21*F11r*u1**2*c1*f1*w2
     &  - 4.D0*PC2*ssp*svm*F21*F10r*u0*u1*c0*f1*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*svm*F20*F15r*u0*u5*c5*f1*w2
     &  - 4.D0*PC2*ssp*svm*F20*F14r*u0*u4*c4*f1*w2
     &  - 4.D0*PC2*ssp*svm*F20*F13r*u0*u3*c3*f1*w2
     &  - 4.D0*PC2*ssp*svm*F20*F12r*u0*u2*c2*f1*w2
     &  - 4.D0*PC2*ssp*svm*F20*F11r*u0*u1*c1*f1*w2
     &  - 4.D0*PC2*ssp*svm*F20*F10r*u0**2*c0*f1*w2
     &  - 4.D0*PC2*ssp*svm*F15*F25r*u5**2*c5*f2*w2
     &  - 4.D0*PC2*ssp*svm*F15*F24r*u4*u5*c4*f2*w2
     &  - 4.D0*PC2*ssp*svm*F15*F23r*u3*u5*c3*f2*w2
     &  - 4.D0*PC2*ssp*svm*F15*F22r*u2*u5*c2*f2*w2
     &  - 4.D0*PC2*ssp*svm*F15*F21r*u1*u5*c1*f2*w2
     &  - 4.D0*PC2*ssp*svm*F15*F20r*u0*u5*c0*f2*w2
     &  - 4.D0*PC2*ssp*svm*F14*F25r*u4*u5*c5*f2*w2
     &  - 4.D0*PC2*ssp*svm*F14*F24r*u4**2*c4*f2*w2
     &  - 4.D0*PC2*ssp*svm*F14*F23r*u3*u4*c3*f2*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*svm*F14*F22r*u2*u4*c2*f2*w2
     &  - 4.D0*PC2*ssp*svm*F14*F21r*u1*u4*c1*f2*w2
     &  - 4.D0*PC2*ssp*svm*F14*F20r*u0*u4*c0*f2*w2
     &  - 4.D0*PC2*ssp*svm*F13*F25r*u3*u5*c5*f2*w2
     &  - 4.D0*PC2*ssp*svm*F13*F24r*u3*u4*c4*f2*w2
     &  - 4.D0*PC2*ssp*svm*F13*F23r*u3**2*c3*f2*w2
     &  - 4.D0*PC2*ssp*svm*F13*F22r*u2*u3*c2*f2*w2
     &  - 4.D0*PC2*ssp*svm*F13*F21r*u1*u3*c1*f2*w2
     &  - 4.D0*PC2*ssp*svm*F13*F20r*u0*u3*c0*f2*w2
     &  - 4.D0*PC2*ssp*svm*F12*F25r*u2*u5*c5*f2*w2
     &  - 4.D0*PC2*ssp*svm*F12*F24r*u2*u4*c4*f2*w2
     &  - 4.D0*PC2*ssp*svm*F12*F23r*u2*u3*c3*f2*w2
     &  - 4.D0*PC2*ssp*svm*F12*F22r*u2**2*c2*f2*w2
     &  - 4.D0*PC2*ssp*svm*F12*F21r*u1*u2*c1*f2*w2
     &  - 4.D0*PC2*ssp*svm*F12*F20r*u0*u2*c0*f2*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*svm*F11*F25r*u1*u5*c5*f2*w2
     &  - 4.D0*PC2*ssp*svm*F11*F24r*u1*u4*c4*f2*w2
     &  - 4.D0*PC2*ssp*svm*F11*F23r*u1*u3*c3*f2*w2
     &  - 4.D0*PC2*ssp*svm*F11*F22r*u1*u2*c2*f2*w2
     &  - 4.D0*PC2*ssp*svm*F11*F21r*u1**2*c1*f2*w2
     &  - 4.D0*PC2*ssp*svm*F11*F20r*u0*u1*c0*f2*w2
     &  - 4.D0*PC2*ssp*svm*F10*F25r*u0*u5*c5*f2*w2
     &  - 4.D0*PC2*ssp*svm*F10*F24r*u0*u4*c4*f2*w2
     &  - 4.D0*PC2*ssp*svm*F10*F23r*u0*u3*c3*f2*w2
     &  - 4.D0*PC2*ssp*svm*F10*F22r*u0*u2*c2*f2*w2
     &  - 4.D0*PC2*ssp*svm*F10*F21r*u0*u1*c1*f2*w2
     &  - 4.D0*PC2*ssp*svm*F10*F20r*u0**2*c0*f2*w2
     &  + 4.D0*PC2*ssp*ssm*F25*F25r*u5**2*c5*f2
     &  + 4.D0*PC2*ssp*ssm*F25*F24r*u4*u5*c4*f2
     &  + 4.D0*PC2*ssp*ssm*F25*F23r*u3*u5*c3*f2
     &
      traza1 = traza1 + 4.D0*PC2*ssp*ssm*F25*F22r*u2*u5*c2*f2
     &  + 4.D0*PC2*ssp*ssm*F25*F21r*u1*u5*c1*f2
     &  + 4.D0*PC2*ssp*ssm*F25*F20r*u0*u5*c0*f2
     &  + 4.D0*PC2*ssp*ssm*F24*F25r*u4*u5*c5*f2
     &  + 4.D0*PC2*ssp*ssm*F24*F24r*u4**2*c4*f2
     &  + 4.D0*PC2*ssp*ssm*F24*F23r*u3*u4*c3*f2
     &  + 4.D0*PC2*ssp*ssm*F24*F22r*u2*u4*c2*f2
     &  + 4.D0*PC2*ssp*ssm*F24*F21r*u1*u4*c1*f2
     &  + 4.D0*PC2*ssp*ssm*F24*F20r*u0*u4*c0*f2
     &  + 4.D0*PC2*ssp*ssm*F23*F25r*u3*u5*c5*f2
     &  + 4.D0*PC2*ssp*ssm*F23*F24r*u3*u4*c4*f2
     &  + 4.D0*PC2*ssp*ssm*F23*F23r*u3**2*c3*f2
     &  + 4.D0*PC2*ssp*ssm*F23*F22r*u2*u3*c2*f2
     &  + 4.D0*PC2*ssp*ssm*F23*F21r*u1*u3*c1*f2
     &  + 4.D0*PC2*ssp*ssm*F23*F20r*u0*u3*c0*f2
     &
      traza1 = traza1 + 4.D0*PC2*ssp*ssm*F22*F25r*u2*u5*c5*f2
     &  + 4.D0*PC2*ssp*ssm*F22*F24r*u2*u4*c4*f2
     &  + 4.D0*PC2*ssp*ssm*F22*F23r*u2*u3*c3*f2
     &  + 4.D0*PC2*ssp*ssm*F22*F22r*u2**2*c2*f2
     &  + 4.D0*PC2*ssp*ssm*F22*F21r*u1*u2*c1*f2
     &  + 4.D0*PC2*ssp*ssm*F22*F20r*u0*u2*c0*f2
     &  + 4.D0*PC2*ssp*ssm*F21*F25r*u1*u5*c5*f2
     &  + 4.D0*PC2*ssp*ssm*F21*F24r*u1*u4*c4*f2
     &  + 4.D0*PC2*ssp*ssm*F21*F23r*u1*u3*c3*f2
     &  + 4.D0*PC2*ssp*ssm*F21*F22r*u1*u2*c2*f2
     &  + 4.D0*PC2*ssp*ssm*F21*F21r*u1**2*c1*f2
     &  + 4.D0*PC2*ssp*ssm*F21*F20r*u0*u1*c0*f2
     &  + 4.D0*PC2*ssp*ssm*F20*F25r*u0*u5*c5*f2
     &  + 4.D0*PC2*ssp*ssm*F20*F24r*u0*u4*c4*f2
     &  + 4.D0*PC2*ssp*ssm*F20*F23r*u0*u3*c3*f2
     &
      traza1 = traza1 + 4.D0*PC2*ssp*ssm*F20*F22r*u0*u2*c2*f2
     &  + 4.D0*PC2*ssp*ssm*F20*F21r*u0*u1*c1*f2
     &  + 4.D0*PC2*ssp*ssm*F20*F20r*u0**2*c0*f2
     &  + 4.D0*PC2*qq*svp*svm*F45*F15r*u5**2*c5*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F45*F15r*u5**2*c5*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F45*F14r*u4*u5*c4*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F45*F14r*u4*u5*c4*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F45*F13r*u3*u5*c3*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F45*F13r*u3*u5*c3*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F45*F12r*u2*u5*c2*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F45*F12r*u2*u5*c2*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F45*F11r*u1*u5*c1*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F45*F11r*u1*u5*c1*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F45*F10r*u0*u5*c0*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F45*F10r*u0*u5*c0*f1*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F44*F15r*u4*u5*c5*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F44*F15r*u4*u5*c5*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F44*F14r*u4**2*c4*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F44*F14r*u4**2*c4*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F44*F13r*u3*u4*c3*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F44*F13r*u3*u4*c3*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F44*F12r*u2*u4*c2*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F44*F12r*u2*u4*c2*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F44*F11r*u1*u4*c1*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F44*F11r*u1*u4*c1*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F44*F10r*u0*u4*c0*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F44*F10r*u0*u4*c0*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F43*F15r*u3*u5*c5*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F43*F15r*u3*u5*c5*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F43*F14r*u3*u4*c4*f1*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F43*F14r*u3*u4*c4*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F43*F13r*u3**2*c3*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F43*F13r*u3**2*c3*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F43*F12r*u2*u3*c2*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F43*F12r*u2*u3*c2*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F43*F11r*u1*u3*c1*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F43*F11r*u1*u3*c1*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F43*F10r*u0*u3*c0*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F43*F10r*u0*u3*c0*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F42*F15r*u2*u5*c5*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F42*F15r*u2*u5*c5*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F42*F14r*u2*u4*c4*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F42*F14r*u2*u4*c4*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F42*F13r*u2*u3*c3*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F42*F13r*u2*u3*c3*f1*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F42*F12r*u2**2*c2*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F42*F12r*u2**2*c2*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F42*F11r*u1*u2*c1*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F42*F11r*u1*u2*c1*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F42*F10r*u0*u2*c0*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F42*F10r*u0*u2*c0*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F41*F15r*u1*u5*c5*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F41*F15r*u1*u5*c5*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F41*F14r*u1*u4*c4*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F41*F14r*u1*u4*c4*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F41*F13r*u1*u3*c3*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F41*F13r*u1*u3*c3*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F41*F12r*u1*u2*c2*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F41*F12r*u1*u2*c2*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F41*F11r*u1**2*c1*f1*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F41*F11r*u1**2*c1*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F41*F10r*u0*u1*c0*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F41*F10r*u0*u1*c0*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F40*F15r*u0*u5*c5*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F40*F15r*u0*u5*c5*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F40*F14r*u0*u4*c4*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F40*F14r*u0*u4*c4*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F40*F13r*u0*u3*c3*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F40*F13r*u0*u3*c3*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F40*F12r*u0*u2*c2*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F40*F12r*u0*u2*c2*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F40*F11r*u0*u1*c1*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F40*F11r*u0*u1*c1*f1*w1
     &  + 4.D0*PC2*qq*svp*svm*F40*F10r*u0**2*c0*f1*w2
     &  + 4.D0*PC2*qq*svp*svm*F40*F10r*u0**2*c0*f1*w1
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*svm*F25*F25r*u5**2*c5*f2
     &  - 4.D0*PC2*qq*svp*svm*F25*F24r*u4*u5*c4*f2
     &  - 4.D0*PC2*qq*svp*svm*F25*F23r*u3*u5*c3*f2
     &  - 4.D0*PC2*qq*svp*svm*F25*F22r*u2*u5*c2*f2
     &  - 4.D0*PC2*qq*svp*svm*F25*F21r*u1*u5*c1*f2
     &  - 4.D0*PC2*qq*svp*svm*F25*F20r*u0*u5*c0*f2
     &  - 4.D0*PC2*qq*svp*svm*F24*F25r*u4*u5*c5*f2
     &  - 4.D0*PC2*qq*svp*svm*F24*F24r*u4**2*c4*f2
     &  - 4.D0*PC2*qq*svp*svm*F24*F23r*u3*u4*c3*f2
     &  - 4.D0*PC2*qq*svp*svm*F24*F22r*u2*u4*c2*f2
     &  - 4.D0*PC2*qq*svp*svm*F24*F21r*u1*u4*c1*f2
     &  - 4.D0*PC2*qq*svp*svm*F24*F20r*u0*u4*c0*f2
     &  - 4.D0*PC2*qq*svp*svm*F23*F25r*u3*u5*c5*f2
     &  - 4.D0*PC2*qq*svp*svm*F23*F24r*u3*u4*c4*f2
     &  - 4.D0*PC2*qq*svp*svm*F23*F23r*u3**2*c3*f2
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*svm*F23*F22r*u2*u3*c2*f2
     &  - 4.D0*PC2*qq*svp*svm*F23*F21r*u1*u3*c1*f2
     &  - 4.D0*PC2*qq*svp*svm*F23*F20r*u0*u3*c0*f2
     &  - 4.D0*PC2*qq*svp*svm*F22*F25r*u2*u5*c5*f2
     &  - 4.D0*PC2*qq*svp*svm*F22*F24r*u2*u4*c4*f2
     &  - 4.D0*PC2*qq*svp*svm*F22*F23r*u2*u3*c3*f2
     &  - 4.D0*PC2*qq*svp*svm*F22*F22r*u2**2*c2*f2
     &  - 4.D0*PC2*qq*svp*svm*F22*F21r*u1*u2*c1*f2
     &  - 4.D0*PC2*qq*svp*svm*F22*F20r*u0*u2*c0*f2
     &  - 4.D0*PC2*qq*svp*svm*F21*F25r*u1*u5*c5*f2
     &  - 4.D0*PC2*qq*svp*svm*F21*F24r*u1*u4*c4*f2
     &  - 4.D0*PC2*qq*svp*svm*F21*F23r*u1*u3*c3*f2
     &  - 4.D0*PC2*qq*svp*svm*F21*F22r*u1*u2*c2*f2
     &  - 4.D0*PC2*qq*svp*svm*F21*F21r*u1**2*c1*f2
     &  - 4.D0*PC2*qq*svp*svm*F21*F20r*u0*u1*c0*f2
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*svm*F20*F25r*u0*u5*c5*f2
     &  - 4.D0*PC2*qq*svp*svm*F20*F24r*u0*u4*c4*f2
     &  - 4.D0*PC2*qq*svp*svm*F20*F23r*u0*u3*c3*f2
     &  - 4.D0*PC2*qq*svp*svm*F20*F22r*u0*u2*c2*f2
     &  - 4.D0*PC2*qq*svp*svm*F20*F21r*u0*u1*c1*f2
     &  - 4.D0*PC2*qq*svp*svm*F20*F20r*u0**2*c0*f2
     &  + 4.D0*PC2*qq*svp*svm*F15*F45r*u5**2*c5*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F15*F45r*u5**2*c5*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F15*F44r*u4*u5*c4*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F15*F44r*u4*u5*c4*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F15*F43r*u3*u5*c3*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F15*F43r*u3*u5*c3*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F15*F42r*u2*u5*c2*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F15*F42r*u2*u5*c2*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F15*F41r*u1*u5*c1*f4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F15*F41r*u1*u5*c1*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F15*F40r*u0*u5*c0*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F15*F40r*u0*u5*c0*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F14*F45r*u4*u5*c5*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F14*F45r*u4*u5*c5*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F14*F44r*u4**2*c4*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F14*F44r*u4**2*c4*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F14*F43r*u3*u4*c3*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F14*F43r*u3*u4*c3*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F14*F42r*u2*u4*c2*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F14*F42r*u2*u4*c2*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F14*F41r*u1*u4*c1*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F14*F41r*u1*u4*c1*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F14*F40r*u0*u4*c0*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F14*F40r*u0*u4*c0*f4*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F13*F45r*u3*u5*c5*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F13*F45r*u3*u5*c5*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F13*F44r*u3*u4*c4*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F13*F44r*u3*u4*c4*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F13*F43r*u3**2*c3*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F13*F43r*u3**2*c3*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F13*F42r*u2*u3*c2*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F13*F42r*u2*u3*c2*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F13*F41r*u1*u3*c1*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F13*F41r*u1*u3*c1*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F13*F40r*u0*u3*c0*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F13*F40r*u0*u3*c0*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F12*F45r*u2*u5*c5*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F12*F45r*u2*u5*c5*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F12*F44r*u2*u4*c4*f4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F12*F44r*u2*u4*c4*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F12*F43r*u2*u3*c3*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F12*F43r*u2*u3*c3*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F12*F42r*u2**2*c2*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F12*F42r*u2**2*c2*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F12*F41r*u1*u2*c1*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F12*F41r*u1*u2*c1*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F12*F40r*u0*u2*c0*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F12*F40r*u0*u2*c0*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F11*F45r*u1*u5*c5*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F11*F45r*u1*u5*c5*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F11*F44r*u1*u4*c4*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F11*F44r*u1*u4*c4*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F11*F43r*u1*u3*c3*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F11*F43r*u1*u3*c3*f4*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F11*F42r*u1*u2*c2*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F11*F42r*u1*u2*c2*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F11*F41r*u1**2*c1*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F11*F41r*u1**2*c1*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F11*F40r*u0*u1*c0*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F11*F40r*u0*u1*c0*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F10*F45r*u0*u5*c5*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F10*F45r*u0*u5*c5*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F10*F44r*u0*u4*c4*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F10*F44r*u0*u4*c4*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F10*F43r*u0*u3*c3*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F10*F43r*u0*u3*c3*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F10*F42r*u0*u2*c2*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F10*F42r*u0*u2*c2*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F10*F41r*u0*u1*c1*f4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*svm*F10*F41r*u0*u1*c1*f4*w1
     &  + 4.D0*PC2*qq*svp*svm*F10*F40r*u0**2*c0*f4*w2
     &  + 4.D0*PC2*qq*svp*svm*F10*F40r*u0**2*c0*f4*w1
     &  - 4.D0*PC2*qq*ssm*svp*F45*F25r*u5**2*c5*f2
     &  - 4.D0*PC2*qq*ssm*svp*F45*F24r*u4*u5*c4*f2
     &  - 4.D0*PC2*qq*ssm*svp*F45*F23r*u3*u5*c3*f2
     &  - 4.D0*PC2*qq*ssm*svp*F45*F22r*u2*u5*c2*f2
     &  - 4.D0*PC2*qq*ssm*svp*F45*F21r*u1*u5*c1*f2
     &  - 4.D0*PC2*qq*ssm*svp*F45*F20r*u0*u5*c0*f2
     &  - 4.D0*PC2*qq*ssm*svp*F44*F25r*u4*u5*c5*f2
     &  - 4.D0*PC2*qq*ssm*svp*F44*F24r*u4**2*c4*f2
     &  - 4.D0*PC2*qq*ssm*svp*F44*F23r*u3*u4*c3*f2
     &  - 4.D0*PC2*qq*ssm*svp*F44*F22r*u2*u4*c2*f2
     &  - 4.D0*PC2*qq*ssm*svp*F44*F21r*u1*u4*c1*f2
     &  - 4.D0*PC2*qq*ssm*svp*F44*F20r*u0*u4*c0*f2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*svp*F43*F25r*u3*u5*c5*f2
     &  - 4.D0*PC2*qq*ssm*svp*F43*F24r*u3*u4*c4*f2
     &  - 4.D0*PC2*qq*ssm*svp*F43*F23r*u3**2*c3*f2
     &  - 4.D0*PC2*qq*ssm*svp*F43*F22r*u2*u3*c2*f2
     &  - 4.D0*PC2*qq*ssm*svp*F43*F21r*u1*u3*c1*f2
     &  - 4.D0*PC2*qq*ssm*svp*F43*F20r*u0*u3*c0*f2
     &  - 4.D0*PC2*qq*ssm*svp*F42*F25r*u2*u5*c5*f2
     &  - 4.D0*PC2*qq*ssm*svp*F42*F24r*u2*u4*c4*f2
     &  - 4.D0*PC2*qq*ssm*svp*F42*F23r*u2*u3*c3*f2
     &  - 4.D0*PC2*qq*ssm*svp*F42*F22r*u2**2*c2*f2
     &  - 4.D0*PC2*qq*ssm*svp*F42*F21r*u1*u2*c1*f2
     &  - 4.D0*PC2*qq*ssm*svp*F42*F20r*u0*u2*c0*f2
     &  - 4.D0*PC2*qq*ssm*svp*F41*F25r*u1*u5*c5*f2
     &  - 4.D0*PC2*qq*ssm*svp*F41*F24r*u1*u4*c4*f2
     &  - 4.D0*PC2*qq*ssm*svp*F41*F23r*u1*u3*c3*f2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*svp*F41*F22r*u1*u2*c2*f2
     &  - 4.D0*PC2*qq*ssm*svp*F41*F21r*u1**2*c1*f2
     &  - 4.D0*PC2*qq*ssm*svp*F41*F20r*u0*u1*c0*f2
     &  - 4.D0*PC2*qq*ssm*svp*F40*F25r*u0*u5*c5*f2
     &  - 4.D0*PC2*qq*ssm*svp*F40*F24r*u0*u4*c4*f2
     &  - 4.D0*PC2*qq*ssm*svp*F40*F23r*u0*u3*c3*f2
     &  - 4.D0*PC2*qq*ssm*svp*F40*F22r*u0*u2*c2*f2
     &  - 4.D0*PC2*qq*ssm*svp*F40*F21r*u0*u1*c1*f2
     &  - 4.D0*PC2*qq*ssm*svp*F40*F20r*u0**2*c0*f2
     &  - 4.D0*PC2*qq*ssm*svp*F25*F45r*u5**2*c5*f4
     &  - 4.D0*PC2*qq*ssm*svp*F25*F44r*u4*u5*c4*f4
     &  - 4.D0*PC2*qq*ssm*svp*F25*F43r*u3*u5*c3*f4
     &  - 4.D0*PC2*qq*ssm*svp*F25*F42r*u2*u5*c2*f4
     &  - 4.D0*PC2*qq*ssm*svp*F25*F41r*u1*u5*c1*f4
     &  - 4.D0*PC2*qq*ssm*svp*F25*F40r*u0*u5*c0*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*svp*F24*F45r*u4*u5*c5*f4
     &  - 4.D0*PC2*qq*ssm*svp*F24*F44r*u4**2*c4*f4
     &  - 4.D0*PC2*qq*ssm*svp*F24*F43r*u3*u4*c3*f4
     &  - 4.D0*PC2*qq*ssm*svp*F24*F42r*u2*u4*c2*f4
     &  - 4.D0*PC2*qq*ssm*svp*F24*F41r*u1*u4*c1*f4
     &  - 4.D0*PC2*qq*ssm*svp*F24*F40r*u0*u4*c0*f4
     &  - 4.D0*PC2*qq*ssm*svp*F23*F45r*u3*u5*c5*f4
     &  - 4.D0*PC2*qq*ssm*svp*F23*F44r*u3*u4*c4*f4
     &  - 4.D0*PC2*qq*ssm*svp*F23*F43r*u3**2*c3*f4
     &  - 4.D0*PC2*qq*ssm*svp*F23*F42r*u2*u3*c2*f4
     &  - 4.D0*PC2*qq*ssm*svp*F23*F41r*u1*u3*c1*f4
     &  - 4.D0*PC2*qq*ssm*svp*F23*F40r*u0*u3*c0*f4
     &  - 4.D0*PC2*qq*ssm*svp*F22*F45r*u2*u5*c5*f4
     &  - 4.D0*PC2*qq*ssm*svp*F22*F44r*u2*u4*c4*f4
     &  - 4.D0*PC2*qq*ssm*svp*F22*F43r*u2*u3*c3*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*svp*F22*F42r*u2**2*c2*f4
     &  - 4.D0*PC2*qq*ssm*svp*F22*F41r*u1*u2*c1*f4
     &  - 4.D0*PC2*qq*ssm*svp*F22*F40r*u0*u2*c0*f4
     &  - 4.D0*PC2*qq*ssm*svp*F21*F45r*u1*u5*c5*f4
     &  - 4.D0*PC2*qq*ssm*svp*F21*F44r*u1*u4*c4*f4
     &  - 4.D0*PC2*qq*ssm*svp*F21*F43r*u1*u3*c3*f4
     &  - 4.D0*PC2*qq*ssm*svp*F21*F42r*u1*u2*c2*f4
     &  - 4.D0*PC2*qq*ssm*svp*F21*F41r*u1**2*c1*f4
     &  - 4.D0*PC2*qq*ssm*svp*F21*F40r*u0*u1*c0*f4
     &  - 4.D0*PC2*qq*ssm*svp*F20*F45r*u0*u5*c5*f4
     &  - 4.D0*PC2*qq*ssm*svp*F20*F44r*u0*u4*c4*f4
     &  - 4.D0*PC2*qq*ssm*svp*F20*F43r*u0*u3*c3*f4
     &  - 4.D0*PC2*qq*ssm*svp*F20*F42r*u0*u2*c2*f4
     &  - 4.D0*PC2*qq*ssm*svp*F20*F41r*u0*u1*c1*f4
     &  - 4.D0*PC2*qq*ssm*svp*F20*F40r*u0**2*c0*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*svm*F45*F25r*u5**2*c5*f2
     &  - 4.D0*PC2*qq*ssp*svm*F45*F24r*u4*u5*c4*f2
     &  - 4.D0*PC2*qq*ssp*svm*F45*F23r*u3*u5*c3*f2
     &  - 4.D0*PC2*qq*ssp*svm*F45*F22r*u2*u5*c2*f2
     &  - 4.D0*PC2*qq*ssp*svm*F45*F21r*u1*u5*c1*f2
     &  - 4.D0*PC2*qq*ssp*svm*F45*F20r*u0*u5*c0*f2
     &  - 4.D0*PC2*qq*ssp*svm*F44*F25r*u4*u5*c5*f2
     &  - 4.D0*PC2*qq*ssp*svm*F44*F24r*u4**2*c4*f2
     &  - 4.D0*PC2*qq*ssp*svm*F44*F23r*u3*u4*c3*f2
     &  - 4.D0*PC2*qq*ssp*svm*F44*F22r*u2*u4*c2*f2
     &  - 4.D0*PC2*qq*ssp*svm*F44*F21r*u1*u4*c1*f2
     &  - 4.D0*PC2*qq*ssp*svm*F44*F20r*u0*u4*c0*f2
     &  - 4.D0*PC2*qq*ssp*svm*F43*F25r*u3*u5*c5*f2
     &  - 4.D0*PC2*qq*ssp*svm*F43*F24r*u3*u4*c4*f2
     &  - 4.D0*PC2*qq*ssp*svm*F43*F23r*u3**2*c3*f2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*svm*F43*F22r*u2*u3*c2*f2
     &  - 4.D0*PC2*qq*ssp*svm*F43*F21r*u1*u3*c1*f2
     &  - 4.D0*PC2*qq*ssp*svm*F43*F20r*u0*u3*c0*f2
     &  - 4.D0*PC2*qq*ssp*svm*F42*F25r*u2*u5*c5*f2
     &  - 4.D0*PC2*qq*ssp*svm*F42*F24r*u2*u4*c4*f2
     &  - 4.D0*PC2*qq*ssp*svm*F42*F23r*u2*u3*c3*f2
     &  - 4.D0*PC2*qq*ssp*svm*F42*F22r*u2**2*c2*f2
     &  - 4.D0*PC2*qq*ssp*svm*F42*F21r*u1*u2*c1*f2
     &  - 4.D0*PC2*qq*ssp*svm*F42*F20r*u0*u2*c0*f2
     &  - 4.D0*PC2*qq*ssp*svm*F41*F25r*u1*u5*c5*f2
     &  - 4.D0*PC2*qq*ssp*svm*F41*F24r*u1*u4*c4*f2
     &  - 4.D0*PC2*qq*ssp*svm*F41*F23r*u1*u3*c3*f2
     &  - 4.D0*PC2*qq*ssp*svm*F41*F22r*u1*u2*c2*f2
     &  - 4.D0*PC2*qq*ssp*svm*F41*F21r*u1**2*c1*f2
     &  - 4.D0*PC2*qq*ssp*svm*F41*F20r*u0*u1*c0*f2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*svm*F40*F25r*u0*u5*c5*f2
     &  - 4.D0*PC2*qq*ssp*svm*F40*F24r*u0*u4*c4*f2
     &  - 4.D0*PC2*qq*ssp*svm*F40*F23r*u0*u3*c3*f2
     &  - 4.D0*PC2*qq*ssp*svm*F40*F22r*u0*u2*c2*f2
     &  - 4.D0*PC2*qq*ssp*svm*F40*F21r*u0*u1*c1*f2
     &  - 4.D0*PC2*qq*ssp*svm*F40*F20r*u0**2*c0*f2
     &  - 4.D0*PC2*qq*ssp*svm*F25*F45r*u5**2*c5*f4
     &  - 4.D0*PC2*qq*ssp*svm*F25*F44r*u4*u5*c4*f4
     &  - 4.D0*PC2*qq*ssp*svm*F25*F43r*u3*u5*c3*f4
     &  - 4.D0*PC2*qq*ssp*svm*F25*F42r*u2*u5*c2*f4
     &  - 4.D0*PC2*qq*ssp*svm*F25*F41r*u1*u5*c1*f4
     &  - 4.D0*PC2*qq*ssp*svm*F25*F40r*u0*u5*c0*f4
     &  - 4.D0*PC2*qq*ssp*svm*F24*F45r*u4*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*svm*F24*F44r*u4**2*c4*f4
     &  - 4.D0*PC2*qq*ssp*svm*F24*F43r*u3*u4*c3*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*svm*F24*F42r*u2*u4*c2*f4
     &  - 4.D0*PC2*qq*ssp*svm*F24*F41r*u1*u4*c1*f4
     &  - 4.D0*PC2*qq*ssp*svm*F24*F40r*u0*u4*c0*f4
     &  - 4.D0*PC2*qq*ssp*svm*F23*F45r*u3*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*svm*F23*F44r*u3*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*svm*F23*F43r*u3**2*c3*f4
     &  - 4.D0*PC2*qq*ssp*svm*F23*F42r*u2*u3*c2*f4
     &  - 4.D0*PC2*qq*ssp*svm*F23*F41r*u1*u3*c1*f4
     &  - 4.D0*PC2*qq*ssp*svm*F23*F40r*u0*u3*c0*f4
     &  - 4.D0*PC2*qq*ssp*svm*F22*F45r*u2*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*svm*F22*F44r*u2*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*svm*F22*F43r*u2*u3*c3*f4
     &  - 4.D0*PC2*qq*ssp*svm*F22*F42r*u2**2*c2*f4
     &  - 4.D0*PC2*qq*ssp*svm*F22*F41r*u1*u2*c1*f4
     &  - 4.D0*PC2*qq*ssp*svm*F22*F40r*u0*u2*c0*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*svm*F21*F45r*u1*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*svm*F21*F44r*u1*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*svm*F21*F43r*u1*u3*c3*f4
     &  - 4.D0*PC2*qq*ssp*svm*F21*F42r*u1*u2*c2*f4
     &  - 4.D0*PC2*qq*ssp*svm*F21*F41r*u1**2*c1*f4
     &  - 4.D0*PC2*qq*ssp*svm*F21*F40r*u0*u1*c0*f4
     &  - 4.D0*PC2*qq*ssp*svm*F20*F45r*u0*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*svm*F20*F44r*u0*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*svm*F20*F43r*u0*u3*c3*f4
     &  - 4.D0*PC2*qq*ssp*svm*F20*F42r*u0*u2*c2*f4
     &  - 4.D0*PC2*qq*ssp*svm*F20*F41r*u0*u1*c1*f4
     &  - 4.D0*PC2*qq*ssp*svm*F20*F40r*u0**2*c0*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F45*F45r*u5**2*c5*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F45*F44r*u4*u5*c4*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F45*F43r*u3*u5*c3*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*ssm*F45*F42r*u2*u5*c2*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F45*F41r*u1*u5*c1*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F45*F40r*u0*u5*c0*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F44*F45r*u4*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F44*F44r*u4**2*c4*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F44*F43r*u3*u4*c3*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F44*F42r*u2*u4*c2*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F44*F41r*u1*u4*c1*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F44*F40r*u0*u4*c0*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F43*F45r*u3*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F43*F44r*u3*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F43*F43r*u3**2*c3*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F43*F42r*u2*u3*c2*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F43*F41r*u1*u3*c1*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F43*F40r*u0*u3*c0*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*ssm*F42*F45r*u2*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F42*F44r*u2*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F42*F43r*u2*u3*c3*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F42*F42r*u2**2*c2*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F42*F41r*u1*u2*c1*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F42*F40r*u0*u2*c0*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F41*F45r*u1*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F41*F44r*u1*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F41*F43r*u1*u3*c3*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F41*F42r*u1*u2*c2*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F41*F41r*u1**2*c1*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F41*F40r*u0*u1*c0*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F40*F45r*u0*u5*c5*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F40*F44r*u0*u4*c4*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F40*F43r*u0*u3*c3*f4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*ssm*F40*F42r*u0*u2*c2*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F40*F41r*u0*u1*c1*f4
     &  - 4.D0*PC2*qq*ssp*ssm*F40*F40r*u0**2*c0*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F45*F45r*u5**2*c5*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F45*F44r*u4*u5*c4*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F45*F43r*u3*u5*c3*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F45*F42r*u2*u5*c2*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F45*F41r*u1*u5*c1*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F45*F40r*u0*u5*c0*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F44*F45r*u4*u5*c5*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F44*F44r*u4**2*c4*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F44*F43r*u3*u4*c3*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F44*F42r*u2*u4*c2*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F44*F41r*u1*u4*c1*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F44*F40r*u0*u4*c0*f4
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svp*svm*F43*F45r*u3*u5*c5*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F43*F44r*u3*u4*c4*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F43*F43r*u3**2*c3*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F43*F42r*u2*u3*c2*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F43*F41r*u1*u3*c1*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F43*F40r*u0*u3*c0*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F42*F45r*u2*u5*c5*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F42*F44r*u2*u4*c4*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F42*F43r*u2*u3*c3*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F42*F42r*u2**2*c2*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F42*F41r*u1*u2*c1*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F42*F40r*u0*u2*c0*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F41*F45r*u1*u5*c5*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F41*F44r*u1*u4*c4*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F41*F43r*u1*u3*c3*f4
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svp*svm*F41*F42r*u1*u2*c2*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F41*F41r*u1**2*c1*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F41*F40r*u0*u1*c0*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F40*F45r*u0*u5*c5*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F40*F44r*u0*u4*c4*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F40*F43r*u0*u3*c3*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F40*F42r*u0*u2*c2*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F40*F41r*u0*u1*c1*f4
     &  + 4.D0*PC2*qq**2*svp*svm*F40*F40r*u0**2*c0*f4
     &  - 4.D0*PC2**2*svp*svm*F25*F25r*u5**2*c5*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F25*F24r*u4*u5*c4*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F25*F23r*u3*u5*c3*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F25*F22r*u2*u5*c2*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F25*F21r*u1*u5*c1*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F25*F20r*u0*u5*c0*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*svp*svm*F24*F25r*u4*u5*c5*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F24*F24r*u4**2*c4*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F24*F23r*u3*u4*c3*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F24*F22r*u2*u4*c2*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F24*F21r*u1*u4*c1*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F24*F20r*u0*u4*c0*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F23*F25r*u3*u5*c5*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F23*F24r*u3*u4*c4*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F23*F23r*u3**2*c3*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F23*F22r*u2*u3*c2*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F23*F21r*u1*u3*c1*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F23*F20r*u0*u3*c0*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F22*F25r*u2*u5*c5*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F22*F24r*u2*u4*c4*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F22*F23r*u2*u3*c3*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*svp*svm*F22*F22r*u2**2*c2*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F22*F21r*u1*u2*c1*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F22*F20r*u0*u2*c0*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F21*F25r*u1*u5*c5*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F21*F24r*u1*u4*c4*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F21*F23r*u1*u3*c3*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F21*F22r*u1*u2*c2*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F21*F21r*u1**2*c1*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F21*F20r*u0*u1*c0*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F20*F25r*u0*u5*c5*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F20*F24r*u0*u4*c4*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F20*F23r*u0*u3*c3*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F20*F22r*u0*u2*c2*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F20*F21r*u0*u1*c1*f2*w1*w2
     &  - 4.D0*PC2**2*svp*svm*F20*F20r*u0**2*c0*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svp*svm*F45*F45r*u5**2*c5*f4*w1*
     & w2
     &  - 4.D0*PC2**2*qq*svp*svm*F45*F44r*u4*u5*c4*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F45*F43r*u3*u5*c3*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F45*F42r*u2*u5*c2*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F45*F41r*u1*u5*c1*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F45*F40r*u0*u5*c0*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F44*F45r*u4*u5*c5*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F44*F44r*u4**2*c4*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F44*F43r*u3*u4*c3*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F44*F42r*u2*u4*c2*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F44*F41r*u1*u4*c1*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F44*F40r*u0*u4*c0*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F43*F45r*u3*u5*c5*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F43*F44r*u3*u4*c4*f4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svp*svm*F43*F43r*u3**2*c3*f4*w1*
     & w2
     &  - 4.D0*PC2**2*qq*svp*svm*F43*F42r*u2*u3*c2*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F43*F41r*u1*u3*c1*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F43*F40r*u0*u3*c0*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F42*F45r*u2*u5*c5*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F42*F44r*u2*u4*c4*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F42*F43r*u2*u3*c3*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F42*F42r*u2**2*c2*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F42*F41r*u1*u2*c1*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F42*F40r*u0*u2*c0*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F41*F45r*u1*u5*c5*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F41*F44r*u1*u4*c4*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F41*F43r*u1*u3*c3*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F41*F42r*u1*u2*c2*f4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svp*svm*F41*F41r*u1**2*c1*f4*w1*
     & w2
     &  - 4.D0*PC2**2*qq*svp*svm*F41*F40r*u0*u1*c0*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F40*F45r*u0*u5*c5*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F40*F44r*u0*u4*c4*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F40*F43r*u0*u3*c3*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F40*F42r*u0*u2*c2*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F40*F41r*u0*u1*c1*f4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*svm*F40*F40r*u0**2*c0*f4*w1*w2
     &  - 4.D0*i_*ssp*ssm*F15*F15i*u5**2*c5*f1
     &  - 4.D0*i_*ssp*ssm*F15*F14i*u4*u5*c4*f1
     &  - 4.D0*i_*ssp*ssm*F15*F13i*u3*u5*c3*f1
     &  - 4.D0*i_*ssp*ssm*F15*F12i*u2*u5*c2*f1
     &  - 4.D0*i_*ssp*ssm*F15*F11i*u1*u5*c1*f1
     &  - 4.D0*i_*ssp*ssm*F15*F10i*u0*u5*c0*f1
     &
      traza1 = traza1 - 4.D0*i_*ssp*ssm*F14*F15i*u4*u5*c5*f1
     &  - 4.D0*i_*ssp*ssm*F14*F14i*u4**2*c4*f1
     &  - 4.D0*i_*ssp*ssm*F14*F13i*u3*u4*c3*f1
     &  - 4.D0*i_*ssp*ssm*F14*F12i*u2*u4*c2*f1
     &  - 4.D0*i_*ssp*ssm*F14*F11i*u1*u4*c1*f1
     &  - 4.D0*i_*ssp*ssm*F14*F10i*u0*u4*c0*f1
     &  - 4.D0*i_*ssp*ssm*F13*F15i*u3*u5*c5*f1
     &  - 4.D0*i_*ssp*ssm*F13*F14i*u3*u4*c4*f1
     &  - 4.D0*i_*ssp*ssm*F13*F13i*u3**2*c3*f1
     &  - 4.D0*i_*ssp*ssm*F13*F12i*u2*u3*c2*f1
     &  - 4.D0*i_*ssp*ssm*F13*F11i*u1*u3*c1*f1
     &  - 4.D0*i_*ssp*ssm*F13*F10i*u0*u3*c0*f1
     &  - 4.D0*i_*ssp*ssm*F12*F15i*u2*u5*c5*f1
     &  - 4.D0*i_*ssp*ssm*F12*F14i*u2*u4*c4*f1
     &  - 4.D0*i_*ssp*ssm*F12*F13i*u2*u3*c3*f1
     &
      traza1 = traza1 - 4.D0*i_*ssp*ssm*F12*F12i*u2**2*c2*f1
     &  - 4.D0*i_*ssp*ssm*F12*F11i*u1*u2*c1*f1
     &  - 4.D0*i_*ssp*ssm*F12*F10i*u0*u2*c0*f1
     &  - 4.D0*i_*ssp*ssm*F11*F15i*u1*u5*c5*f1
     &  - 4.D0*i_*ssp*ssm*F11*F14i*u1*u4*c4*f1
     &  - 4.D0*i_*ssp*ssm*F11*F13i*u1*u3*c3*f1
     &  - 4.D0*i_*ssp*ssm*F11*F12i*u1*u2*c2*f1
     &  - 4.D0*i_*ssp*ssm*F11*F11i*u1**2*c1*f1
     &  - 4.D0*i_*ssp*ssm*F11*F10i*u0*u1*c0*f1
     &  - 4.D0*i_*ssp*ssm*F10*F15i*u0*u5*c5*f1
     &  - 4.D0*i_*ssp*ssm*F10*F14i*u0*u4*c4*f1
     &  - 4.D0*i_*ssp*ssm*F10*F13i*u0*u3*c3*f1
     &  - 4.D0*i_*ssp*ssm*F10*F12i*u0*u2*c2*f1
     &  - 4.D0*i_*ssp*ssm*F10*F11i*u0*u1*c1*f1
     &  - 4.D0*i_*ssp*ssm*F10*F10i*u0**2*c0*f1
     &
      traza1 = traza1 - 4.D0*i_*qq*svp*svm*F15*F15i*u5**2*c5*f1
     &  - 4.D0*i_*qq*svp*svm*F15*F14i*u4*u5*c4*f1
     &  - 4.D0*i_*qq*svp*svm*F15*F13i*u3*u5*c3*f1
     &  - 4.D0*i_*qq*svp*svm*F15*F12i*u2*u5*c2*f1
     &  - 4.D0*i_*qq*svp*svm*F15*F11i*u1*u5*c1*f1
     &  - 4.D0*i_*qq*svp*svm*F15*F10i*u0*u5*c0*f1
     &  - 4.D0*i_*qq*svp*svm*F14*F15i*u4*u5*c5*f1
     &  - 4.D0*i_*qq*svp*svm*F14*F14i*u4**2*c4*f1
     &  - 4.D0*i_*qq*svp*svm*F14*F13i*u3*u4*c3*f1
     &  - 4.D0*i_*qq*svp*svm*F14*F12i*u2*u4*c2*f1
     &  - 4.D0*i_*qq*svp*svm*F14*F11i*u1*u4*c1*f1
     &  - 4.D0*i_*qq*svp*svm*F14*F10i*u0*u4*c0*f1
     &  - 4.D0*i_*qq*svp*svm*F13*F15i*u3*u5*c5*f1
     &  - 4.D0*i_*qq*svp*svm*F13*F14i*u3*u4*c4*f1
     &  - 4.D0*i_*qq*svp*svm*F13*F13i*u3**2*c3*f1
     &
      traza1 = traza1 - 4.D0*i_*qq*svp*svm*F13*F12i*u2*u3*c2*f1
     &  - 4.D0*i_*qq*svp*svm*F13*F11i*u1*u3*c1*f1
     &  - 4.D0*i_*qq*svp*svm*F13*F10i*u0*u3*c0*f1
     &  - 4.D0*i_*qq*svp*svm*F12*F15i*u2*u5*c5*f1
     &  - 4.D0*i_*qq*svp*svm*F12*F14i*u2*u4*c4*f1
     &  - 4.D0*i_*qq*svp*svm*F12*F13i*u2*u3*c3*f1
     &  - 4.D0*i_*qq*svp*svm*F12*F12i*u2**2*c2*f1
     &  - 4.D0*i_*qq*svp*svm*F12*F11i*u1*u2*c1*f1
     &  - 4.D0*i_*qq*svp*svm*F12*F10i*u0*u2*c0*f1
     &  - 4.D0*i_*qq*svp*svm*F11*F15i*u1*u5*c5*f1
     &  - 4.D0*i_*qq*svp*svm*F11*F14i*u1*u4*c4*f1
     &  - 4.D0*i_*qq*svp*svm*F11*F13i*u1*u3*c3*f1
     &  - 4.D0*i_*qq*svp*svm*F11*F12i*u1*u2*c2*f1
     &  - 4.D0*i_*qq*svp*svm*F11*F11i*u1**2*c1*f1
     &  - 4.D0*i_*qq*svp*svm*F11*F10i*u0*u1*c0*f1
     &
      traza1 = traza1 - 4.D0*i_*qq*svp*svm*F10*F15i*u0*u5*c5*f1
     &  - 4.D0*i_*qq*svp*svm*F10*F14i*u0*u4*c4*f1
     &  - 4.D0*i_*qq*svp*svm*F10*F13i*u0*u3*c3*f1
     &  - 4.D0*i_*qq*svp*svm*F10*F12i*u0*u2*c2*f1
     &  - 4.D0*i_*qq*svp*svm*F10*F11i*u0*u1*c1*f1
     &  - 4.D0*i_*qq*svp*svm*F10*F10i*u0**2*c0*f1
     &  + 4.D0*i_*PC2*svp*svm*F15*F15i*u5**2*c5*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F15*F14i*u4*u5*c4*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F15*F13i*u3*u5*c3*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F15*F12i*u2*u5*c2*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F15*F11i*u1*u5*c1*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F15*F10i*u0*u5*c0*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F14*F15i*u4*u5*c5*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F14*F14i*u4**2*c4*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F14*F13i*u3*u4*c3*f1*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svp*svm*F14*F12i*u2*u4*c2*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F14*F11i*u1*u4*c1*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F14*F10i*u0*u4*c0*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F13*F15i*u3*u5*c5*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F13*F14i*u3*u4*c4*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F13*F13i*u3**2*c3*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F13*F12i*u2*u3*c2*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F13*F11i*u1*u3*c1*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F13*F10i*u0*u3*c0*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F12*F15i*u2*u5*c5*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F12*F14i*u2*u4*c4*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F12*F13i*u2*u3*c3*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F12*F12i*u2**2*c2*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F12*F11i*u1*u2*c1*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F12*F10i*u0*u2*c0*f1*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svp*svm*F11*F15i*u1*u5*c5*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F11*F14i*u1*u4*c4*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F11*F13i*u1*u3*c3*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F11*F12i*u1*u2*c2*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F11*F11i*u1**2*c1*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F11*F10i*u0*u1*c0*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F10*F15i*u0*u5*c5*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F10*F14i*u0*u4*c4*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F10*F13i*u0*u3*c3*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F10*F12i*u0*u2*c2*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F10*F11i*u0*u1*c1*f1*w1*w2
     &  + 4.D0*i_*PC2*svp*svm*F10*F10i*u0**2*c0*f1*w1*w2
     &  - 4.D0*i_*PC2*ssm*svp*F25*F15i*u5**2*c5*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F25*F14i*u4*u5*c4*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F25*F13i*u3*u5*c3*f1*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*svp*F25*F12i*u2*u5*c2*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F25*F11i*u1*u5*c1*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F25*F10i*u0*u5*c0*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F24*F15i*u4*u5*c5*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F24*F14i*u4**2*c4*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F24*F13i*u3*u4*c3*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F24*F12i*u2*u4*c2*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F24*F11i*u1*u4*c1*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F24*F10i*u0*u4*c0*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F23*F15i*u3*u5*c5*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F23*F14i*u3*u4*c4*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F23*F13i*u3**2*c3*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F23*F12i*u2*u3*c2*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F23*F11i*u1*u3*c1*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F23*F10i*u0*u3*c0*f1*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*svp*F22*F15i*u2*u5*c5*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F22*F14i*u2*u4*c4*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F22*F13i*u2*u3*c3*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F22*F12i*u2**2*c2*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F22*F11i*u1*u2*c1*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F22*F10i*u0*u2*c0*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F21*F15i*u1*u5*c5*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F21*F14i*u1*u4*c4*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F21*F13i*u1*u3*c3*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F21*F12i*u1*u2*c2*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F21*F11i*u1**2*c1*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F21*F10i*u0*u1*c0*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F20*F15i*u0*u5*c5*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F20*F14i*u0*u4*c4*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F20*F13i*u0*u3*c3*f1*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*svp*F20*F12i*u0*u2*c2*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F20*F11i*u0*u1*c1*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F20*F10i*u0**2*c0*f1*w1
     &  - 4.D0*i_*PC2*ssm*svp*F15*F25i*u5**2*c5*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F15*F24i*u4*u5*c4*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F15*F23i*u3*u5*c3*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F15*F22i*u2*u5*c2*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F15*F21i*u1*u5*c1*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F15*F20i*u0*u5*c0*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F14*F25i*u4*u5*c5*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F14*F24i*u4**2*c4*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F14*F23i*u3*u4*c3*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F14*F22i*u2*u4*c2*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F14*F21i*u1*u4*c1*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F14*F20i*u0*u4*c0*f2*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*svp*F13*F25i*u3*u5*c5*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F13*F24i*u3*u4*c4*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F13*F23i*u3**2*c3*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F13*F22i*u2*u3*c2*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F13*F21i*u1*u3*c1*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F13*F20i*u0*u3*c0*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F12*F25i*u2*u5*c5*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F12*F24i*u2*u4*c4*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F12*F23i*u2*u3*c3*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F12*F22i*u2**2*c2*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F12*F21i*u1*u2*c1*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F12*F20i*u0*u2*c0*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F11*F25i*u1*u5*c5*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F11*F24i*u1*u4*c4*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F11*F23i*u1*u3*c3*f2*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*svp*F11*F22i*u1*u2*c2*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F11*F21i*u1**2*c1*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F11*F20i*u0*u1*c0*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F10*F25i*u0*u5*c5*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F10*F24i*u0*u4*c4*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F10*F23i*u0*u3*c3*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F10*F22i*u0*u2*c2*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F10*F21i*u0*u1*c1*f2*w1
     &  - 4.D0*i_*PC2*ssm*svp*F10*F20i*u0**2*c0*f2*w1
     &  - 4.D0*i_*PC2*ssp*svm*F25*F15i*u5**2*c5*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F25*F14i*u4*u5*c4*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F25*F13i*u3*u5*c3*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F25*F12i*u2*u5*c2*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F25*F11i*u1*u5*c1*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F25*F10i*u0*u5*c0*f1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*svm*F24*F15i*u4*u5*c5*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F24*F14i*u4**2*c4*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F24*F13i*u3*u4*c3*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F24*F12i*u2*u4*c2*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F24*F11i*u1*u4*c1*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F24*F10i*u0*u4*c0*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F23*F15i*u3*u5*c5*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F23*F14i*u3*u4*c4*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F23*F13i*u3**2*c3*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F23*F12i*u2*u3*c2*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F23*F11i*u1*u3*c1*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F23*F10i*u0*u3*c0*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F22*F15i*u2*u5*c5*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F22*F14i*u2*u4*c4*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F22*F13i*u2*u3*c3*f1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*svm*F22*F12i*u2**2*c2*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F22*F11i*u1*u2*c1*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F22*F10i*u0*u2*c0*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F21*F15i*u1*u5*c5*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F21*F14i*u1*u4*c4*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F21*F13i*u1*u3*c3*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F21*F12i*u1*u2*c2*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F21*F11i*u1**2*c1*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F21*F10i*u0*u1*c0*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F20*F15i*u0*u5*c5*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F20*F14i*u0*u4*c4*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F20*F13i*u0*u3*c3*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F20*F12i*u0*u2*c2*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F20*F11i*u0*u1*c1*f1*w2
     &  - 4.D0*i_*PC2*ssp*svm*F20*F10i*u0**2*c0*f1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*svm*F15*F25i*u5**2*c5*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F15*F24i*u4*u5*c4*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F15*F23i*u3*u5*c3*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F15*F22i*u2*u5*c2*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F15*F21i*u1*u5*c1*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F15*F20i*u0*u5*c0*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F14*F25i*u4*u5*c5*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F14*F24i*u4**2*c4*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F14*F23i*u3*u4*c3*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F14*F22i*u2*u4*c2*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F14*F21i*u1*u4*c1*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F14*F20i*u0*u4*c0*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F13*F25i*u3*u5*c5*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F13*F24i*u3*u4*c4*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F13*F23i*u3**2*c3*f2*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*svm*F13*F22i*u2*u3*c2*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F13*F21i*u1*u3*c1*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F13*F20i*u0*u3*c0*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F12*F25i*u2*u5*c5*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F12*F24i*u2*u4*c4*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F12*F23i*u2*u3*c3*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F12*F22i*u2**2*c2*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F12*F21i*u1*u2*c1*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F12*F20i*u0*u2*c0*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F11*F25i*u1*u5*c5*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F11*F24i*u1*u4*c4*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F11*F23i*u1*u3*c3*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F11*F22i*u1*u2*c2*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F11*F21i*u1**2*c1*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F11*F20i*u0*u1*c0*f2*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*svm*F10*F25i*u0*u5*c5*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F10*F24i*u0*u4*c4*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F10*F23i*u0*u3*c3*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F10*F22i*u0*u2*c2*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F10*F21i*u0*u1*c1*f2*w2
     &  - 4.D0*i_*PC2*ssp*svm*F10*F20i*u0**2*c0*f2*w2
     &  + 4.D0*i_*PC2*ssp*ssm*F25*F25i*u5**2*c5*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F25*F24i*u4*u5*c4*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F25*F23i*u3*u5*c3*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F25*F22i*u2*u5*c2*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F25*F21i*u1*u5*c1*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F25*F20i*u0*u5*c0*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F24*F25i*u4*u5*c5*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F24*F24i*u4**2*c4*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F24*F23i*u3*u4*c3*f2
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssp*ssm*F24*F22i*u2*u4*c2*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F24*F21i*u1*u4*c1*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F24*F20i*u0*u4*c0*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F23*F25i*u3*u5*c5*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F23*F24i*u3*u4*c4*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F23*F23i*u3**2*c3*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F23*F22i*u2*u3*c2*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F23*F21i*u1*u3*c1*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F23*F20i*u0*u3*c0*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F22*F25i*u2*u5*c5*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F22*F24i*u2*u4*c4*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F22*F23i*u2*u3*c3*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F22*F22i*u2**2*c2*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F22*F21i*u1*u2*c1*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F22*F20i*u0*u2*c0*f2
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssp*ssm*F21*F25i*u1*u5*c5*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F21*F24i*u1*u4*c4*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F21*F23i*u1*u3*c3*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F21*F22i*u1*u2*c2*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F21*F21i*u1**2*c1*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F21*F20i*u0*u1*c0*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F20*F25i*u0*u5*c5*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F20*F24i*u0*u4*c4*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F20*F23i*u0*u3*c3*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F20*F22i*u0*u2*c2*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F20*F21i*u0*u1*c1*f2
     &  + 4.D0*i_*PC2*ssp*ssm*F20*F20i*u0**2*c0*f2
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F15i*u5**2*c5*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F15i*u5**2*c5*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F14i*u4*u5*c4*f1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F45*F14i*u4*u5*c4*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F13i*u3*u5*c3*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F13i*u3*u5*c3*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F12i*u2*u5*c2*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F12i*u2*u5*c2*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F11i*u1*u5*c1*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F11i*u1*u5*c1*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F10i*u0*u5*c0*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F45*F10i*u0*u5*c0*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F15i*u4*u5*c5*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F15i*u4*u5*c5*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F14i*u4**2*c4*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F14i*u4**2*c4*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F13i*u3*u4*c3*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F13i*u3*u4*c3*f1*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F44*F12i*u2*u4*c2*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F12i*u2*u4*c2*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F11i*u1*u4*c1*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F11i*u1*u4*c1*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F10i*u0*u4*c0*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F44*F10i*u0*u4*c0*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F15i*u3*u5*c5*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F15i*u3*u5*c5*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F14i*u3*u4*c4*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F14i*u3*u4*c4*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F13i*u3**2*c3*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F13i*u3**2*c3*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F12i*u2*u3*c2*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F12i*u2*u3*c2*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F11i*u1*u3*c1*f1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F43*F11i*u1*u3*c1*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F10i*u0*u3*c0*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F43*F10i*u0*u3*c0*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F15i*u2*u5*c5*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F15i*u2*u5*c5*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F14i*u2*u4*c4*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F14i*u2*u4*c4*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F13i*u2*u3*c3*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F13i*u2*u3*c3*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F12i*u2**2*c2*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F12i*u2**2*c2*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F11i*u1*u2*c1*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F11i*u1*u2*c1*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F10i*u0*u2*c0*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F42*F10i*u0*u2*c0*f1*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F41*F15i*u1*u5*c5*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F15i*u1*u5*c5*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F14i*u1*u4*c4*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F14i*u1*u4*c4*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F13i*u1*u3*c3*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F13i*u1*u3*c3*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F12i*u1*u2*c2*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F12i*u1*u2*c2*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F11i*u1**2*c1*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F11i*u1**2*c1*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F10i*u0*u1*c0*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F41*F10i*u0*u1*c0*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F15i*u0*u5*c5*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F15i*u0*u5*c5*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F14i*u0*u4*c4*f1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F40*F14i*u0*u4*c4*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F13i*u0*u3*c3*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F13i*u0*u3*c3*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F12i*u0*u2*c2*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F12i*u0*u2*c2*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F11i*u0*u1*c1*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F11i*u0*u1*c1*f1*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F10i*u0**2*c0*f1*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F40*F10i*u0**2*c0*f1*w1
     &  - 4.D0*i_*PC2*qq*svp*svm*F25*F25i*u5**2*c5*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F25*F24i*u4*u5*c4*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F25*F23i*u3*u5*c3*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F25*F22i*u2*u5*c2*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F25*F21i*u1*u5*c1*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F25*F20i*u0*u5*c0*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*svm*F24*F25i*u4*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F24*F24i*u4**2*c4*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F24*F23i*u3*u4*c3*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F24*F22i*u2*u4*c2*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F24*F21i*u1*u4*c1*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F24*F20i*u0*u4*c0*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F23*F25i*u3*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F23*F24i*u3*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F23*F23i*u3**2*c3*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F23*F22i*u2*u3*c2*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F23*F21i*u1*u3*c1*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F23*F20i*u0*u3*c0*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F22*F25i*u2*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F22*F24i*u2*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F22*F23i*u2*u3*c3*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*svm*F22*F22i*u2**2*c2*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F22*F21i*u1*u2*c1*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F22*F20i*u0*u2*c0*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F21*F25i*u1*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F21*F24i*u1*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F21*F23i*u1*u3*c3*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F21*F22i*u1*u2*c2*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F21*F21i*u1**2*c1*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F21*F20i*u0*u1*c0*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F20*F25i*u0*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F20*F24i*u0*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F20*F23i*u0*u3*c3*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F20*F22i*u0*u2*c2*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F20*F21i*u0*u1*c1*f2
     &  - 4.D0*i_*PC2*qq*svp*svm*F20*F20i*u0**2*c0*f2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F15*F45i*u5**2*c5*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F45i*u5**2*c5*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F44i*u4*u5*c4*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F44i*u4*u5*c4*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F43i*u3*u5*c3*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F43i*u3*u5*c3*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F42i*u2*u5*c2*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F42i*u2*u5*c2*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F41i*u1*u5*c1*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F41i*u1*u5*c1*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F40i*u0*u5*c0*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F15*F40i*u0*u5*c0*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F45i*u4*u5*c5*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F45i*u4*u5*c5*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F44i*u4**2*c4*f4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F14*F44i*u4**2*c4*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F43i*u3*u4*c3*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F43i*u3*u4*c3*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F42i*u2*u4*c2*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F42i*u2*u4*c2*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F41i*u1*u4*c1*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F41i*u1*u4*c1*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F40i*u0*u4*c0*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F14*F40i*u0*u4*c0*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F45i*u3*u5*c5*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F45i*u3*u5*c5*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F44i*u3*u4*c4*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F44i*u3*u4*c4*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F43i*u3**2*c3*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F43i*u3**2*c3*f4*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F13*F42i*u2*u3*c2*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F42i*u2*u3*c2*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F41i*u1*u3*c1*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F41i*u1*u3*c1*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F40i*u0*u3*c0*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F13*F40i*u0*u3*c0*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F45i*u2*u5*c5*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F45i*u2*u5*c5*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F44i*u2*u4*c4*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F44i*u2*u4*c4*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F43i*u2*u3*c3*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F43i*u2*u3*c3*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F42i*u2**2*c2*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F42i*u2**2*c2*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F41i*u1*u2*c1*f4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F12*F41i*u1*u2*c1*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F40i*u0*u2*c0*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F12*F40i*u0*u2*c0*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F45i*u1*u5*c5*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F45i*u1*u5*c5*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F44i*u1*u4*c4*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F44i*u1*u4*c4*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F43i*u1*u3*c3*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F43i*u1*u3*c3*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F42i*u1*u2*c2*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F42i*u1*u2*c2*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F41i*u1**2*c1*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F41i*u1**2*c1*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F40i*u0*u1*c0*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F11*F40i*u0*u1*c0*f4*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*svm*F10*F45i*u0*u5*c5*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F45i*u0*u5*c5*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F44i*u0*u4*c4*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F44i*u0*u4*c4*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F43i*u0*u3*c3*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F43i*u0*u3*c3*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F42i*u0*u2*c2*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F42i*u0*u2*c2*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F41i*u0*u1*c1*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F41i*u0*u1*c1*f4*w1
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F40i*u0**2*c0*f4*w2
     &  + 4.D0*i_*PC2*qq*svp*svm*F10*F40i*u0**2*c0*f4*w1
     &  - 4.D0*i_*PC2*qq*ssm*svp*F45*F25i*u5**2*c5*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F45*F24i*u4*u5*c4*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F45*F23i*u3*u5*c3*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*svp*F45*F22i*u2*u5*c2*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F45*F21i*u1*u5*c1*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F45*F20i*u0*u5*c0*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F44*F25i*u4*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F44*F24i*u4**2*c4*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F44*F23i*u3*u4*c3*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F44*F22i*u2*u4*c2*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F44*F21i*u1*u4*c1*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F44*F20i*u0*u4*c0*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F43*F25i*u3*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F43*F24i*u3*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F43*F23i*u3**2*c3*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F43*F22i*u2*u3*c2*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F43*F21i*u1*u3*c1*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F43*F20i*u0*u3*c0*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*svp*F42*F25i*u2*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F42*F24i*u2*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F42*F23i*u2*u3*c3*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F42*F22i*u2**2*c2*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F42*F21i*u1*u2*c1*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F42*F20i*u0*u2*c0*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F41*F25i*u1*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F41*F24i*u1*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F41*F23i*u1*u3*c3*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F41*F22i*u1*u2*c2*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F41*F21i*u1**2*c1*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F41*F20i*u0*u1*c0*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F40*F25i*u0*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F40*F24i*u0*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F40*F23i*u0*u3*c3*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*svp*F40*F22i*u0*u2*c2*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F40*F21i*u0*u1*c1*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F40*F20i*u0**2*c0*f2
     &  - 4.D0*i_*PC2*qq*ssm*svp*F25*F45i*u5**2*c5*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F25*F44i*u4*u5*c4*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F25*F43i*u3*u5*c3*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F25*F42i*u2*u5*c2*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F25*F41i*u1*u5*c1*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F25*F40i*u0*u5*c0*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F24*F45i*u4*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F24*F44i*u4**2*c4*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F24*F43i*u3*u4*c3*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F24*F42i*u2*u4*c2*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F24*F41i*u1*u4*c1*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F24*F40i*u0*u4*c0*f4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*svp*F23*F45i*u3*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F23*F44i*u3*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F23*F43i*u3**2*c3*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F23*F42i*u2*u3*c2*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F23*F41i*u1*u3*c1*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F23*F40i*u0*u3*c0*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F22*F45i*u2*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F22*F44i*u2*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F22*F43i*u2*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F22*F42i*u2**2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F22*F41i*u1*u2*c1*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F22*F40i*u0*u2*c0*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F21*F45i*u1*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F21*F44i*u1*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F21*F43i*u1*u3*c3*f4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*svp*F21*F42i*u1*u2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F21*F41i*u1**2*c1*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F21*F40i*u0*u1*c0*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F20*F45i*u0*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F20*F44i*u0*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F20*F43i*u0*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F20*F42i*u0*u2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F20*F41i*u0*u1*c1*f4
     &  - 4.D0*i_*PC2*qq*ssm*svp*F20*F40i*u0**2*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F45*F25i*u5**2*c5*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F45*F24i*u4*u5*c4*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F45*F23i*u3*u5*c3*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F45*F22i*u2*u5*c2*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F45*F21i*u1*u5*c1*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F45*F20i*u0*u5*c0*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*svm*F44*F25i*u4*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F44*F24i*u4**2*c4*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F44*F23i*u3*u4*c3*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F44*F22i*u2*u4*c2*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F44*F21i*u1*u4*c1*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F44*F20i*u0*u4*c0*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F43*F25i*u3*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F43*F24i*u3*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F43*F23i*u3**2*c3*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F43*F22i*u2*u3*c2*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F43*F21i*u1*u3*c1*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F43*F20i*u0*u3*c0*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F42*F25i*u2*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F42*F24i*u2*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F42*F23i*u2*u3*c3*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*svm*F42*F22i*u2**2*c2*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F42*F21i*u1*u2*c1*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F42*F20i*u0*u2*c0*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F41*F25i*u1*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F41*F24i*u1*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F41*F23i*u1*u3*c3*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F41*F22i*u1*u2*c2*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F41*F21i*u1**2*c1*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F41*F20i*u0*u1*c0*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F40*F25i*u0*u5*c5*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F40*F24i*u0*u4*c4*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F40*F23i*u0*u3*c3*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F40*F22i*u0*u2*c2*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F40*F21i*u0*u1*c1*f2
     &  - 4.D0*i_*PC2*qq*ssp*svm*F40*F20i*u0**2*c0*f2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*svm*F25*F45i*u5**2*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F25*F44i*u4*u5*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F25*F43i*u3*u5*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F25*F42i*u2*u5*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F25*F41i*u1*u5*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F25*F40i*u0*u5*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F24*F45i*u4*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F24*F44i*u4**2*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F24*F43i*u3*u4*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F24*F42i*u2*u4*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F24*F41i*u1*u4*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F24*F40i*u0*u4*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F23*F45i*u3*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F23*F44i*u3*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F23*F43i*u3**2*c3*f4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*svm*F23*F42i*u2*u3*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F23*F41i*u1*u3*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F23*F40i*u0*u3*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F22*F45i*u2*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F22*F44i*u2*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F22*F43i*u2*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F22*F42i*u2**2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F22*F41i*u1*u2*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F22*F40i*u0*u2*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F21*F45i*u1*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F21*F44i*u1*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F21*F43i*u1*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F21*F42i*u1*u2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F21*F41i*u1**2*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F21*F40i*u0*u1*c0*f4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*svm*F20*F45i*u0*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F20*F44i*u0*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F20*F43i*u0*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F20*F42i*u0*u2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F20*F41i*u0*u1*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*svm*F20*F40i*u0**2*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F45*F45i*u5**2*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F45*F44i*u4*u5*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F45*F43i*u3*u5*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F45*F42i*u2*u5*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F45*F41i*u1*u5*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F45*F40i*u0*u5*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F44*F45i*u4*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F44*F44i*u4**2*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F44*F43i*u3*u4*c3*f4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*ssm*F44*F42i*u2*u4*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F44*F41i*u1*u4*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F44*F40i*u0*u4*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F43*F45i*u3*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F43*F44i*u3*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F43*F43i*u3**2*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F43*F42i*u2*u3*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F43*F41i*u1*u3*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F43*F40i*u0*u3*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F42*F45i*u2*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F42*F44i*u2*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F42*F43i*u2*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F42*F42i*u2**2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F42*F41i*u1*u2*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F42*F40i*u0*u2*c0*f4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*ssm*F41*F45i*u1*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F41*F44i*u1*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F41*F43i*u1*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F41*F42i*u1*u2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F41*F41i*u1**2*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F41*F40i*u0*u1*c0*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F40*F45i*u0*u5*c5*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F40*F44i*u0*u4*c4*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F40*F43i*u0*u3*c3*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F40*F42i*u0*u2*c2*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F40*F41i*u0*u1*c1*f4
     &  - 4.D0*i_*PC2*qq*ssp*ssm*F40*F40i*u0**2*c0*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F45*F45i*u5**2*c5*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F45*F44i*u4*u5*c4*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F45*F43i*u3*u5*c3*f4
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svp*svm*F45*F42i*u2*u5*c2*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F45*F41i*u1*u5*c1*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F45*F40i*u0*u5*c0*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F44*F45i*u4*u5*c5*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F44*F44i*u4**2*c4*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F44*F43i*u3*u4*c3*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F44*F42i*u2*u4*c2*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F44*F41i*u1*u4*c1*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F44*F40i*u0*u4*c0*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F43*F45i*u3*u5*c5*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F43*F44i*u3*u4*c4*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F43*F43i*u3**2*c3*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F43*F42i*u2*u3*c2*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F43*F41i*u1*u3*c1*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F43*F40i*u0*u3*c0*f4
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svp*svm*F42*F45i*u2*u5*c5*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F42*F44i*u2*u4*c4*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F42*F43i*u2*u3*c3*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F42*F42i*u2**2*c2*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F42*F41i*u1*u2*c1*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F42*F40i*u0*u2*c0*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F41*F45i*u1*u5*c5*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F41*F44i*u1*u4*c4*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F41*F43i*u1*u3*c3*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F41*F42i*u1*u2*c2*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F41*F41i*u1**2*c1*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F41*F40i*u0*u1*c0*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F40*F45i*u0*u5*c5*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F40*F44i*u0*u4*c4*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F40*F43i*u0*u3*c3*f4
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svp*svm*F40*F42i*u0*u2*c2*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F40*F41i*u0*u1*c1*f4
     &  + 4.D0*i_*PC2*qq**2*svp*svm*F40*F40i*u0**2*c0*f4
     &  - 4.D0*i_*PC2**2*svp*svm*F25*F25i*u5**2*c5*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F25*F24i*u4*u5*c4*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F25*F23i*u3*u5*c3*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F25*F22i*u2*u5*c2*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F25*F21i*u1*u5*c1*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F25*F20i*u0*u5*c0*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F24*F25i*u4*u5*c5*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F24*F24i*u4**2*c4*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F24*F23i*u3*u4*c3*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F24*F22i*u2*u4*c2*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F24*F21i*u1*u4*c1*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F24*F20i*u0*u4*c0*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svp*svm*F23*F25i*u3*u5*c5*f2*w1*
     & w2
     &  - 4.D0*i_*PC2**2*svp*svm*F23*F24i*u3*u4*c4*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F23*F23i*u3**2*c3*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F23*F22i*u2*u3*c2*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F23*F21i*u1*u3*c1*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F23*F20i*u0*u3*c0*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F22*F25i*u2*u5*c5*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F22*F24i*u2*u4*c4*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F22*F23i*u2*u3*c3*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F22*F22i*u2**2*c2*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F22*F21i*u1*u2*c1*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F22*F20i*u0*u2*c0*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F21*F25i*u1*u5*c5*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F21*F24i*u1*u4*c4*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svp*svm*F21*F23i*u1*u3*c3*f2*w1*
     & w2
     &  - 4.D0*i_*PC2**2*svp*svm*F21*F22i*u1*u2*c2*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F21*F21i*u1**2*c1*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F21*F20i*u0*u1*c0*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F20*F25i*u0*u5*c5*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F20*F24i*u0*u4*c4*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F20*F23i*u0*u3*c3*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F20*F22i*u0*u2*c2*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F20*F21i*u0*u1*c1*f2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*svm*F20*F20i*u0**2*c0*f2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F45*F45i*u5**2*c5*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F45*F44i*u4*u5*c4*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F45*F43i*u3*u5*c3*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F45*F42i*u2*u5*c2*f4*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svp*svm*F45*F41i*u1*u5*c1*f4*
     & w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F45*F40i*u0*u5*c0*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F44*F45i*u4*u5*c5*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F44*F44i*u4**2*c4*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F44*F43i*u3*u4*c3*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F44*F42i*u2*u4*c2*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F44*F41i*u1*u4*c1*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F44*F40i*u0*u4*c0*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F43*F45i*u3*u5*c5*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F43*F44i*u3*u4*c4*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F43*F43i*u3**2*c3*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F43*F42i*u2*u3*c2*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F43*F41i*u1*u3*c1*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F43*F40i*u0*u3*c0*f4*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svp*svm*F42*F45i*u2*u5*c5*f4*
     & w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F42*F44i*u2*u4*c4*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F42*F43i*u2*u3*c3*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F42*F42i*u2**2*c2*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F42*F41i*u1*u2*c1*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F42*F40i*u0*u2*c0*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F41*F45i*u1*u5*c5*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F41*F44i*u1*u4*c4*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F41*F43i*u1*u3*c3*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F41*F42i*u1*u2*c2*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F41*F41i*u1**2*c1*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F41*F40i*u0*u1*c0*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F40*F45i*u0*u5*c5*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F40*F44i*u0*u4*c4*f4*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svp*svm*F40*F43i*u0*u3*c3*f4*
     & w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F40*F42i*u0*u2*c2*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F40*F41i*u0*u1*c1*f4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*svm*F40*F40i*u0**2*c0*f4*w1*w2
     &  + 4.D0*PC_q*svp*svm*F15*F15r*u5**2*c5*f1*w2
     &  - 4.D0*PC_q*svp*svm*F15*F15r*u5**2*c5*f1*w1
     &  + 4.D0*PC_q*svp*svm*F15*F14r*u4*u5*c4*f1*w2
     &  - 4.D0*PC_q*svp*svm*F15*F14r*u4*u5*c4*f1*w1
     &  + 4.D0*PC_q*svp*svm*F15*F13r*u3*u5*c3*f1*w2
     &  - 4.D0*PC_q*svp*svm*F15*F13r*u3*u5*c3*f1*w1
     &  + 4.D0*PC_q*svp*svm*F15*F12r*u2*u5*c2*f1*w2
     &  - 4.D0*PC_q*svp*svm*F15*F12r*u2*u5*c2*f1*w1
     &  + 4.D0*PC_q*svp*svm*F15*F11r*u1*u5*c1*f1*w2
     &  - 4.D0*PC_q*svp*svm*F15*F11r*u1*u5*c1*f1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svp*svm*F15*F10r*u0*u5*c0*f1*w2
     &  - 4.D0*PC_q*svp*svm*F15*F10r*u0*u5*c0*f1*w1
     &  + 4.D0*PC_q*svp*svm*F14*F15r*u4*u5*c5*f1*w2
     &  - 4.D0*PC_q*svp*svm*F14*F15r*u4*u5*c5*f1*w1
     &  + 4.D0*PC_q*svp*svm*F14*F14r*u4**2*c4*f1*w2
     &  - 4.D0*PC_q*svp*svm*F14*F14r*u4**2*c4*f1*w1
     &  + 4.D0*PC_q*svp*svm*F14*F13r*u3*u4*c3*f1*w2
     &  - 4.D0*PC_q*svp*svm*F14*F13r*u3*u4*c3*f1*w1
     &  + 4.D0*PC_q*svp*svm*F14*F12r*u2*u4*c2*f1*w2
     &  - 4.D0*PC_q*svp*svm*F14*F12r*u2*u4*c2*f1*w1
     &  + 4.D0*PC_q*svp*svm*F14*F11r*u1*u4*c1*f1*w2
     &  - 4.D0*PC_q*svp*svm*F14*F11r*u1*u4*c1*f1*w1
     &  + 4.D0*PC_q*svp*svm*F14*F10r*u0*u4*c0*f1*w2
     &  - 4.D0*PC_q*svp*svm*F14*F10r*u0*u4*c0*f1*w1
     &  + 4.D0*PC_q*svp*svm*F13*F15r*u3*u5*c5*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*svm*F13*F15r*u3*u5*c5*f1*w1
     &  + 4.D0*PC_q*svp*svm*F13*F14r*u3*u4*c4*f1*w2
     &  - 4.D0*PC_q*svp*svm*F13*F14r*u3*u4*c4*f1*w1
     &  + 4.D0*PC_q*svp*svm*F13*F13r*u3**2*c3*f1*w2
     &  - 4.D0*PC_q*svp*svm*F13*F13r*u3**2*c3*f1*w1
     &  + 4.D0*PC_q*svp*svm*F13*F12r*u2*u3*c2*f1*w2
     &  - 4.D0*PC_q*svp*svm*F13*F12r*u2*u3*c2*f1*w1
     &  + 4.D0*PC_q*svp*svm*F13*F11r*u1*u3*c1*f1*w2
     &  - 4.D0*PC_q*svp*svm*F13*F11r*u1*u3*c1*f1*w1
     &  + 4.D0*PC_q*svp*svm*F13*F10r*u0*u3*c0*f1*w2
     &  - 4.D0*PC_q*svp*svm*F13*F10r*u0*u3*c0*f1*w1
     &  + 4.D0*PC_q*svp*svm*F12*F15r*u2*u5*c5*f1*w2
     &  - 4.D0*PC_q*svp*svm*F12*F15r*u2*u5*c5*f1*w1
     &  + 4.D0*PC_q*svp*svm*F12*F14r*u2*u4*c4*f1*w2
     &  - 4.D0*PC_q*svp*svm*F12*F14r*u2*u4*c4*f1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svp*svm*F12*F13r*u2*u3*c3*f1*w2
     &  - 4.D0*PC_q*svp*svm*F12*F13r*u2*u3*c3*f1*w1
     &  + 4.D0*PC_q*svp*svm*F12*F12r*u2**2*c2*f1*w2
     &  - 4.D0*PC_q*svp*svm*F12*F12r*u2**2*c2*f1*w1
     &  + 4.D0*PC_q*svp*svm*F12*F11r*u1*u2*c1*f1*w2
     &  - 4.D0*PC_q*svp*svm*F12*F11r*u1*u2*c1*f1*w1
     &  + 4.D0*PC_q*svp*svm*F12*F10r*u0*u2*c0*f1*w2
     &  - 4.D0*PC_q*svp*svm*F12*F10r*u0*u2*c0*f1*w1
     &  + 4.D0*PC_q*svp*svm*F11*F15r*u1*u5*c5*f1*w2
     &  - 4.D0*PC_q*svp*svm*F11*F15r*u1*u5*c5*f1*w1
     &  + 4.D0*PC_q*svp*svm*F11*F14r*u1*u4*c4*f1*w2
     &  - 4.D0*PC_q*svp*svm*F11*F14r*u1*u4*c4*f1*w1
     &  + 4.D0*PC_q*svp*svm*F11*F13r*u1*u3*c3*f1*w2
     &  - 4.D0*PC_q*svp*svm*F11*F13r*u1*u3*c3*f1*w1
     &  + 4.D0*PC_q*svp*svm*F11*F12r*u1*u2*c2*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*svm*F11*F12r*u1*u2*c2*f1*w1
     &  + 4.D0*PC_q*svp*svm*F11*F11r*u1**2*c1*f1*w2
     &  - 4.D0*PC_q*svp*svm*F11*F11r*u1**2*c1*f1*w1
     &  + 4.D0*PC_q*svp*svm*F11*F10r*u0*u1*c0*f1*w2
     &  - 4.D0*PC_q*svp*svm*F11*F10r*u0*u1*c0*f1*w1
     &  + 4.D0*PC_q*svp*svm*F10*F15r*u0*u5*c5*f1*w2
     &  - 4.D0*PC_q*svp*svm*F10*F15r*u0*u5*c5*f1*w1
     &  + 4.D0*PC_q*svp*svm*F10*F14r*u0*u4*c4*f1*w2
     &  - 4.D0*PC_q*svp*svm*F10*F14r*u0*u4*c4*f1*w1
     &  + 4.D0*PC_q*svp*svm*F10*F13r*u0*u3*c3*f1*w2
     &  - 4.D0*PC_q*svp*svm*F10*F13r*u0*u3*c3*f1*w1
     &  + 4.D0*PC_q*svp*svm*F10*F12r*u0*u2*c2*f1*w2
     &  - 4.D0*PC_q*svp*svm*F10*F12r*u0*u2*c2*f1*w1
     &  + 4.D0*PC_q*svp*svm*F10*F11r*u0*u1*c1*f1*w2
     &  - 4.D0*PC_q*svp*svm*F10*F11r*u0*u1*c1*f1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svp*svm*F10*F10r*u0**2*c0*f1*w2
     &  - 4.D0*PC_q*svp*svm*F10*F10r*u0**2*c0*f1*w1
     &  - 4.D0*PC_q*ssm*svp*F25*F15r*u5**2*c5*f1
     &  - 4.D0*PC_q*ssm*svp*F25*F14r*u4*u5*c4*f1
     &  - 4.D0*PC_q*ssm*svp*F25*F13r*u3*u5*c3*f1
     &  - 4.D0*PC_q*ssm*svp*F25*F12r*u2*u5*c2*f1
     &  - 4.D0*PC_q*ssm*svp*F25*F11r*u1*u5*c1*f1
     &  - 4.D0*PC_q*ssm*svp*F25*F10r*u0*u5*c0*f1
     &  - 4.D0*PC_q*ssm*svp*F24*F15r*u4*u5*c5*f1
     &  - 4.D0*PC_q*ssm*svp*F24*F14r*u4**2*c4*f1
     &  - 4.D0*PC_q*ssm*svp*F24*F13r*u3*u4*c3*f1
     &  - 4.D0*PC_q*ssm*svp*F24*F12r*u2*u4*c2*f1
     &  - 4.D0*PC_q*ssm*svp*F24*F11r*u1*u4*c1*f1
     &  - 4.D0*PC_q*ssm*svp*F24*F10r*u0*u4*c0*f1
     &  - 4.D0*PC_q*ssm*svp*F23*F15r*u3*u5*c5*f1
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*svp*F23*F14r*u3*u4*c4*f1
     &  - 4.D0*PC_q*ssm*svp*F23*F13r*u3**2*c3*f1
     &  - 4.D0*PC_q*ssm*svp*F23*F12r*u2*u3*c2*f1
     &  - 4.D0*PC_q*ssm*svp*F23*F11r*u1*u3*c1*f1
     &  - 4.D0*PC_q*ssm*svp*F23*F10r*u0*u3*c0*f1
     &  - 4.D0*PC_q*ssm*svp*F22*F15r*u2*u5*c5*f1
     &  - 4.D0*PC_q*ssm*svp*F22*F14r*u2*u4*c4*f1
     &  - 4.D0*PC_q*ssm*svp*F22*F13r*u2*u3*c3*f1
     &  - 4.D0*PC_q*ssm*svp*F22*F12r*u2**2*c2*f1
     &  - 4.D0*PC_q*ssm*svp*F22*F11r*u1*u2*c1*f1
     &  - 4.D0*PC_q*ssm*svp*F22*F10r*u0*u2*c0*f1
     &  - 4.D0*PC_q*ssm*svp*F21*F15r*u1*u5*c5*f1
     &  - 4.D0*PC_q*ssm*svp*F21*F14r*u1*u4*c4*f1
     &  - 4.D0*PC_q*ssm*svp*F21*F13r*u1*u3*c3*f1
     &  - 4.D0*PC_q*ssm*svp*F21*F12r*u1*u2*c2*f1
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*svp*F21*F11r*u1**2*c1*f1
     &  - 4.D0*PC_q*ssm*svp*F21*F10r*u0*u1*c0*f1
     &  - 4.D0*PC_q*ssm*svp*F20*F15r*u0*u5*c5*f1
     &  - 4.D0*PC_q*ssm*svp*F20*F14r*u0*u4*c4*f1
     &  - 4.D0*PC_q*ssm*svp*F20*F13r*u0*u3*c3*f1
     &  - 4.D0*PC_q*ssm*svp*F20*F12r*u0*u2*c2*f1
     &  - 4.D0*PC_q*ssm*svp*F20*F11r*u0*u1*c1*f1
     &  - 4.D0*PC_q*ssm*svp*F20*F10r*u0**2*c0*f1
     &  - 4.D0*PC_q*ssm*svp*F15*F25r*u5**2*c5*f2
     &  - 4.D0*PC_q*ssm*svp*F15*F24r*u4*u5*c4*f2
     &  - 4.D0*PC_q*ssm*svp*F15*F23r*u3*u5*c3*f2
     &  - 4.D0*PC_q*ssm*svp*F15*F22r*u2*u5*c2*f2
     &  - 4.D0*PC_q*ssm*svp*F15*F21r*u1*u5*c1*f2
     &  - 4.D0*PC_q*ssm*svp*F15*F20r*u0*u5*c0*f2
     &  - 4.D0*PC_q*ssm*svp*F14*F25r*u4*u5*c5*f2
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*svp*F14*F24r*u4**2*c4*f2
     &  - 4.D0*PC_q*ssm*svp*F14*F23r*u3*u4*c3*f2
     &  - 4.D0*PC_q*ssm*svp*F14*F22r*u2*u4*c2*f2
     &  - 4.D0*PC_q*ssm*svp*F14*F21r*u1*u4*c1*f2
     &  - 4.D0*PC_q*ssm*svp*F14*F20r*u0*u4*c0*f2
     &  - 4.D0*PC_q*ssm*svp*F13*F25r*u3*u5*c5*f2
     &  - 4.D0*PC_q*ssm*svp*F13*F24r*u3*u4*c4*f2
     &  - 4.D0*PC_q*ssm*svp*F13*F23r*u3**2*c3*f2
     &  - 4.D0*PC_q*ssm*svp*F13*F22r*u2*u3*c2*f2
     &  - 4.D0*PC_q*ssm*svp*F13*F21r*u1*u3*c1*f2
     &  - 4.D0*PC_q*ssm*svp*F13*F20r*u0*u3*c0*f2
     &  - 4.D0*PC_q*ssm*svp*F12*F25r*u2*u5*c5*f2
     &  - 4.D0*PC_q*ssm*svp*F12*F24r*u2*u4*c4*f2
     &  - 4.D0*PC_q*ssm*svp*F12*F23r*u2*u3*c3*f2
     &  - 4.D0*PC_q*ssm*svp*F12*F22r*u2**2*c2*f2
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*svp*F12*F21r*u1*u2*c1*f2
     &  - 4.D0*PC_q*ssm*svp*F12*F20r*u0*u2*c0*f2
     &  - 4.D0*PC_q*ssm*svp*F11*F25r*u1*u5*c5*f2
     &  - 4.D0*PC_q*ssm*svp*F11*F24r*u1*u4*c4*f2
     &  - 4.D0*PC_q*ssm*svp*F11*F23r*u1*u3*c3*f2
     &  - 4.D0*PC_q*ssm*svp*F11*F22r*u1*u2*c2*f2
     &  - 4.D0*PC_q*ssm*svp*F11*F21r*u1**2*c1*f2
     &  - 4.D0*PC_q*ssm*svp*F11*F20r*u0*u1*c0*f2
     &  - 4.D0*PC_q*ssm*svp*F10*F25r*u0*u5*c5*f2
     &  - 4.D0*PC_q*ssm*svp*F10*F24r*u0*u4*c4*f2
     &  - 4.D0*PC_q*ssm*svp*F10*F23r*u0*u3*c3*f2
     &  - 4.D0*PC_q*ssm*svp*F10*F22r*u0*u2*c2*f2
     &  - 4.D0*PC_q*ssm*svp*F10*F21r*u0*u1*c1*f2
     &  - 4.D0*PC_q*ssm*svp*F10*F20r*u0**2*c0*f2
     &  + 4.D0*PC_q*ssp*svm*F25*F15r*u5**2*c5*f1
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*svm*F25*F14r*u4*u5*c4*f1
     &  + 4.D0*PC_q*ssp*svm*F25*F13r*u3*u5*c3*f1
     &  + 4.D0*PC_q*ssp*svm*F25*F12r*u2*u5*c2*f1
     &  + 4.D0*PC_q*ssp*svm*F25*F11r*u1*u5*c1*f1
     &  + 4.D0*PC_q*ssp*svm*F25*F10r*u0*u5*c0*f1
     &  + 4.D0*PC_q*ssp*svm*F24*F15r*u4*u5*c5*f1
     &  + 4.D0*PC_q*ssp*svm*F24*F14r*u4**2*c4*f1
     &  + 4.D0*PC_q*ssp*svm*F24*F13r*u3*u4*c3*f1
     &  + 4.D0*PC_q*ssp*svm*F24*F12r*u2*u4*c2*f1
     &  + 4.D0*PC_q*ssp*svm*F24*F11r*u1*u4*c1*f1
     &  + 4.D0*PC_q*ssp*svm*F24*F10r*u0*u4*c0*f1
     &  + 4.D0*PC_q*ssp*svm*F23*F15r*u3*u5*c5*f1
     &  + 4.D0*PC_q*ssp*svm*F23*F14r*u3*u4*c4*f1
     &  + 4.D0*PC_q*ssp*svm*F23*F13r*u3**2*c3*f1
     &  + 4.D0*PC_q*ssp*svm*F23*F12r*u2*u3*c2*f1
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*svm*F23*F11r*u1*u3*c1*f1
     &  + 4.D0*PC_q*ssp*svm*F23*F10r*u0*u3*c0*f1
     &  + 4.D0*PC_q*ssp*svm*F22*F15r*u2*u5*c5*f1
     &  + 4.D0*PC_q*ssp*svm*F22*F14r*u2*u4*c4*f1
     &  + 4.D0*PC_q*ssp*svm*F22*F13r*u2*u3*c3*f1
     &  + 4.D0*PC_q*ssp*svm*F22*F12r*u2**2*c2*f1
     &  + 4.D0*PC_q*ssp*svm*F22*F11r*u1*u2*c1*f1
     &  + 4.D0*PC_q*ssp*svm*F22*F10r*u0*u2*c0*f1
     &  + 4.D0*PC_q*ssp*svm*F21*F15r*u1*u5*c5*f1
     &  + 4.D0*PC_q*ssp*svm*F21*F14r*u1*u4*c4*f1
     &  + 4.D0*PC_q*ssp*svm*F21*F13r*u1*u3*c3*f1
     &  + 4.D0*PC_q*ssp*svm*F21*F12r*u1*u2*c2*f1
     &  + 4.D0*PC_q*ssp*svm*F21*F11r*u1**2*c1*f1
     &  + 4.D0*PC_q*ssp*svm*F21*F10r*u0*u1*c0*f1
     &  + 4.D0*PC_q*ssp*svm*F20*F15r*u0*u5*c5*f1
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*svm*F20*F14r*u0*u4*c4*f1
     &  + 4.D0*PC_q*ssp*svm*F20*F13r*u0*u3*c3*f1
     &  + 4.D0*PC_q*ssp*svm*F20*F12r*u0*u2*c2*f1
     &  + 4.D0*PC_q*ssp*svm*F20*F11r*u0*u1*c1*f1
     &  + 4.D0*PC_q*ssp*svm*F20*F10r*u0**2*c0*f1
     &  + 4.D0*PC_q*ssp*svm*F15*F25r*u5**2*c5*f2
     &  + 4.D0*PC_q*ssp*svm*F15*F24r*u4*u5*c4*f2
     &  + 4.D0*PC_q*ssp*svm*F15*F23r*u3*u5*c3*f2
     &  + 4.D0*PC_q*ssp*svm*F15*F22r*u2*u5*c2*f2
     &  + 4.D0*PC_q*ssp*svm*F15*F21r*u1*u5*c1*f2
     &  + 4.D0*PC_q*ssp*svm*F15*F20r*u0*u5*c0*f2
     &  + 4.D0*PC_q*ssp*svm*F14*F25r*u4*u5*c5*f2
     &  + 4.D0*PC_q*ssp*svm*F14*F24r*u4**2*c4*f2
     &  + 4.D0*PC_q*ssp*svm*F14*F23r*u3*u4*c3*f2
     &  + 4.D0*PC_q*ssp*svm*F14*F22r*u2*u4*c2*f2
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*svm*F14*F21r*u1*u4*c1*f2
     &  + 4.D0*PC_q*ssp*svm*F14*F20r*u0*u4*c0*f2
     &  + 4.D0*PC_q*ssp*svm*F13*F25r*u3*u5*c5*f2
     &  + 4.D0*PC_q*ssp*svm*F13*F24r*u3*u4*c4*f2
     &  + 4.D0*PC_q*ssp*svm*F13*F23r*u3**2*c3*f2
     &  + 4.D0*PC_q*ssp*svm*F13*F22r*u2*u3*c2*f2
     &  + 4.D0*PC_q*ssp*svm*F13*F21r*u1*u3*c1*f2
     &  + 4.D0*PC_q*ssp*svm*F13*F20r*u0*u3*c0*f2
     &  + 4.D0*PC_q*ssp*svm*F12*F25r*u2*u5*c5*f2
     &  + 4.D0*PC_q*ssp*svm*F12*F24r*u2*u4*c4*f2
     &  + 4.D0*PC_q*ssp*svm*F12*F23r*u2*u3*c3*f2
     &  + 4.D0*PC_q*ssp*svm*F12*F22r*u2**2*c2*f2
     &  + 4.D0*PC_q*ssp*svm*F12*F21r*u1*u2*c1*f2
     &  + 4.D0*PC_q*ssp*svm*F12*F20r*u0*u2*c0*f2
     &  + 4.D0*PC_q*ssp*svm*F11*F25r*u1*u5*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*svm*F11*F24r*u1*u4*c4*f2
     &  + 4.D0*PC_q*ssp*svm*F11*F23r*u1*u3*c3*f2
     &  + 4.D0*PC_q*ssp*svm*F11*F22r*u1*u2*c2*f2
     &  + 4.D0*PC_q*ssp*svm*F11*F21r*u1**2*c1*f2
     &  + 4.D0*PC_q*ssp*svm*F11*F20r*u0*u1*c0*f2
     &  + 4.D0*PC_q*ssp*svm*F10*F25r*u0*u5*c5*f2
     &  + 4.D0*PC_q*ssp*svm*F10*F24r*u0*u4*c4*f2
     &  + 4.D0*PC_q*ssp*svm*F10*F23r*u0*u3*c3*f2
     &  + 4.D0*PC_q*ssp*svm*F10*F22r*u0*u2*c2*f2
     &  + 4.D0*PC_q*ssp*svm*F10*F21r*u0*u1*c1*f2
     &  + 4.D0*PC_q*ssp*svm*F10*F20r*u0**2*c0*f2
     &  - 4.D0*PC_q*qq*ssm*svp*F35*F15r*u5**2*c5*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F35*F14r*u4*u5*c4*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F35*F13r*u3*u5*c3*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F35*F12r*u2*u5*c2*f1
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*svp*F35*F11r*u1*u5*c1*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F35*F10r*u0*u5*c0*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F34*F15r*u4*u5*c5*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F34*F14r*u4**2*c4*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F34*F13r*u3*u4*c3*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F34*F12r*u2*u4*c2*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F34*F11r*u1*u4*c1*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F34*F10r*u0*u4*c0*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F33*F15r*u3*u5*c5*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F33*F14r*u3*u4*c4*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F33*F13r*u3**2*c3*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F33*F12r*u2*u3*c2*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F33*F11r*u1*u3*c1*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F33*F10r*u0*u3*c0*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F32*F15r*u2*u5*c5*f1
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*svp*F32*F14r*u2*u4*c4*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F32*F13r*u2*u3*c3*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F32*F12r*u2**2*c2*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F32*F11r*u1*u2*c1*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F32*F10r*u0*u2*c0*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F31*F15r*u1*u5*c5*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F31*F14r*u1*u4*c4*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F31*F13r*u1*u3*c3*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F31*F12r*u1*u2*c2*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F31*F11r*u1**2*c1*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F31*F10r*u0*u1*c0*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F30*F15r*u0*u5*c5*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F30*F14r*u0*u4*c4*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F30*F13r*u0*u3*c3*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F30*F12r*u0*u2*c2*f1
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*svp*F30*F11r*u0*u1*c1*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F30*F10r*u0**2*c0*f1
     &  - 4.D0*PC_q*qq*ssm*svp*F15*F35r*u5**2*c5*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F15*F34r*u4*u5*c4*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F15*F33r*u3*u5*c3*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F15*F32r*u2*u5*c2*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F15*F31r*u1*u5*c1*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F15*F30r*u0*u5*c0*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F14*F35r*u4*u5*c5*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F14*F34r*u4**2*c4*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F14*F33r*u3*u4*c3*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F14*F32r*u2*u4*c2*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F14*F31r*u1*u4*c1*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F14*F30r*u0*u4*c0*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F13*F35r*u3*u5*c5*f3
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*svp*F13*F34r*u3*u4*c4*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F13*F33r*u3**2*c3*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F13*F32r*u2*u3*c2*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F13*F31r*u1*u3*c1*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F13*F30r*u0*u3*c0*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F12*F35r*u2*u5*c5*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F12*F34r*u2*u4*c4*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F12*F33r*u2*u3*c3*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F12*F32r*u2**2*c2*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F12*F31r*u1*u2*c1*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F12*F30r*u0*u2*c0*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F11*F35r*u1*u5*c5*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F11*F34r*u1*u4*c4*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F11*F33r*u1*u3*c3*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F11*F32r*u1*u2*c2*f3
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*svp*F11*F31r*u1**2*c1*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F11*F30r*u0*u1*c0*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F10*F35r*u0*u5*c5*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F10*F34r*u0*u4*c4*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F10*F33r*u0*u3*c3*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F10*F32r*u0*u2*c2*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F10*F31r*u0*u1*c1*f3
     &  - 4.D0*PC_q*qq*ssm*svp*F10*F30r*u0**2*c0*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F35*F15r*u5**2*c5*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F35*F14r*u4*u5*c4*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F35*F13r*u3*u5*c3*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F35*F12r*u2*u5*c2*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F35*F11r*u1*u5*c1*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F35*F10r*u0*u5*c0*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F34*F15r*u4*u5*c5*f1
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*svm*F34*F14r*u4**2*c4*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F34*F13r*u3*u4*c3*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F34*F12r*u2*u4*c2*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F34*F11r*u1*u4*c1*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F34*F10r*u0*u4*c0*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F33*F15r*u3*u5*c5*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F33*F14r*u3*u4*c4*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F33*F13r*u3**2*c3*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F33*F12r*u2*u3*c2*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F33*F11r*u1*u3*c1*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F33*F10r*u0*u3*c0*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F32*F15r*u2*u5*c5*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F32*F14r*u2*u4*c4*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F32*F13r*u2*u3*c3*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F32*F12r*u2**2*c2*f1
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*svm*F32*F11r*u1*u2*c1*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F32*F10r*u0*u2*c0*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F31*F15r*u1*u5*c5*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F31*F14r*u1*u4*c4*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F31*F13r*u1*u3*c3*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F31*F12r*u1*u2*c2*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F31*F11r*u1**2*c1*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F31*F10r*u0*u1*c0*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F30*F15r*u0*u5*c5*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F30*F14r*u0*u4*c4*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F30*F13r*u0*u3*c3*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F30*F12r*u0*u2*c2*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F30*F11r*u0*u1*c1*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F30*F10r*u0**2*c0*f1
     &  + 4.D0*PC_q*qq*ssp*svm*F15*F35r*u5**2*c5*f3
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*svm*F15*F34r*u4*u5*c4*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F15*F33r*u3*u5*c3*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F15*F32r*u2*u5*c2*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F15*F31r*u1*u5*c1*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F15*F30r*u0*u5*c0*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F14*F35r*u4*u5*c5*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F14*F34r*u4**2*c4*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F14*F33r*u3*u4*c3*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F14*F32r*u2*u4*c2*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F14*F31r*u1*u4*c1*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F14*F30r*u0*u4*c0*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F13*F35r*u3*u5*c5*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F13*F34r*u3*u4*c4*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F13*F33r*u3**2*c3*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F13*F32r*u2*u3*c2*f3
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*svm*F13*F31r*u1*u3*c1*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F13*F30r*u0*u3*c0*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F12*F35r*u2*u5*c5*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F12*F34r*u2*u4*c4*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F12*F33r*u2*u3*c3*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F12*F32r*u2**2*c2*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F12*F31r*u1*u2*c1*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F12*F30r*u0*u2*c0*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F11*F35r*u1*u5*c5*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F11*F34r*u1*u4*c4*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F11*F33r*u1*u3*c3*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F11*F32r*u1*u2*c2*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F11*F31r*u1**2*c1*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F11*F30r*u0*u1*c0*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F10*F35r*u0*u5*c5*f3
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*svm*F10*F34r*u0*u4*c4*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F10*F33r*u0*u3*c3*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F10*F32r*u0*u2*c2*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F10*F31r*u0*u1*c1*f3
     &  + 4.D0*PC_q*qq*ssp*svm*F10*F30r*u0**2*c0*f3
     &  - 4.D0*PC_q*PC2*svp*svm*F25*F25r*u5**2*c5*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F25*F25r*u5**2*c5*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F25*F24r*u4*u5*c4*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F25*F24r*u4*u5*c4*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F25*F23r*u3*u5*c3*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F25*F23r*u3*u5*c3*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F25*F22r*u2*u5*c2*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F25*F22r*u2*u5*c2*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F25*F21r*u1*u5*c1*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F25*F21r*u1*u5*c1*f2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svp*svm*F25*F20r*u0*u5*c0*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F25*F20r*u0*u5*c0*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F24*F25r*u4*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F24*F25r*u4*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F24*F24r*u4**2*c4*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F24*F24r*u4**2*c4*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F24*F23r*u3*u4*c3*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F24*F23r*u3*u4*c3*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F24*F22r*u2*u4*c2*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F24*F22r*u2*u4*c2*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F24*F21r*u1*u4*c1*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F24*F21r*u1*u4*c1*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F24*F20r*u0*u4*c0*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F24*F20r*u0*u4*c0*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F23*F25r*u3*u5*c5*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svp*svm*F23*F25r*u3*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F23*F24r*u3*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F23*F24r*u3*u4*c4*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F23*F23r*u3**2*c3*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F23*F23r*u3**2*c3*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F23*F22r*u2*u3*c2*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F23*F22r*u2*u3*c2*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F23*F21r*u1*u3*c1*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F23*F21r*u1*u3*c1*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F23*F20r*u0*u3*c0*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F23*F20r*u0*u3*c0*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F22*F25r*u2*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F22*F25r*u2*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F22*F24r*u2*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F22*F24r*u2*u4*c4*f2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svp*svm*F22*F23r*u2*u3*c3*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F22*F23r*u2*u3*c3*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F22*F22r*u2**2*c2*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F22*F22r*u2**2*c2*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F22*F21r*u1*u2*c1*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F22*F21r*u1*u2*c1*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F22*F20r*u0*u2*c0*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F22*F20r*u0*u2*c0*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F21*F25r*u1*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F21*F25r*u1*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F21*F24r*u1*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F21*F24r*u1*u4*c4*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F21*F23r*u1*u3*c3*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F21*F23r*u1*u3*c3*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F21*F22r*u1*u2*c2*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svp*svm*F21*F22r*u1*u2*c2*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F21*F21r*u1**2*c1*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F21*F21r*u1**2*c1*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F21*F20r*u0*u1*c0*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F21*F20r*u0*u1*c0*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F20*F25r*u0*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F20*F25r*u0*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F20*F24r*u0*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F20*F24r*u0*u4*c4*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F20*F23r*u0*u3*c3*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F20*F23r*u0*u3*c3*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F20*F22r*u0*u2*c2*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F20*F22r*u0*u2*c2*f2*w1
     &  - 4.D0*PC_q*PC2*svp*svm*F20*F21r*u0*u1*c1*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F20*F21r*u0*u1*c1*f2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svp*svm*F20*F20r*u0**2*c0*f2*w2
     &  + 4.D0*PC_q*PC2*svp*svm*F20*F20r*u0**2*c0*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F45*F45r*u5**2*c5*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F45*F45r*u5**2*c5*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F45*F44r*u4*u5*c4*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F45*F44r*u4*u5*c4*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F45*F43r*u3*u5*c3*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F45*F43r*u3*u5*c3*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F45*F42r*u2*u5*c2*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F45*F42r*u2*u5*c2*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F45*F41r*u1*u5*c1*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F45*F41r*u1*u5*c1*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F45*F40r*u0*u5*c0*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F45*F40r*u0*u5*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F44*F45r*u4*u5*c5*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F44*F45r*u4*u5*c5*f4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F44*F44r*u4**2*c4*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F44*F44r*u4**2*c4*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F44*F43r*u3*u4*c3*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F44*F43r*u3*u4*c3*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F44*F42r*u2*u4*c2*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F44*F42r*u2*u4*c2*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F44*F41r*u1*u4*c1*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F44*F41r*u1*u4*c1*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F44*F40r*u0*u4*c0*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F44*F40r*u0*u4*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F43*F45r*u3*u5*c5*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F43*F45r*u3*u5*c5*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F43*F44r*u3*u4*c4*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F43*F44r*u3*u4*c4*f4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F43*F43r*u3**2*c3*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F43*F43r*u3**2*c3*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F43*F42r*u2*u3*c2*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F43*F42r*u2*u3*c2*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F43*F41r*u1*u3*c1*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F43*F41r*u1*u3*c1*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F43*F40r*u0*u3*c0*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F43*F40r*u0*u3*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F42*F45r*u2*u5*c5*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F42*F45r*u2*u5*c5*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F42*F44r*u2*u4*c4*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F42*F44r*u2*u4*c4*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F42*F43r*u2*u3*c3*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F42*F43r*u2*u3*c3*f4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F42*F42r*u2**2*c2*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F42*F42r*u2**2*c2*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F42*F41r*u1*u2*c1*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F42*F41r*u1*u2*c1*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F42*F40r*u0*u2*c0*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F42*F40r*u0*u2*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F41*F45r*u1*u5*c5*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F41*F45r*u1*u5*c5*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F41*F44r*u1*u4*c4*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F41*F44r*u1*u4*c4*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F41*F43r*u1*u3*c3*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F41*F43r*u1*u3*c3*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F41*F42r*u1*u2*c2*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F41*F42r*u1*u2*c2*f4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F41*F41r*u1**2*c1*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F41*F41r*u1**2*c1*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F41*F40r*u0*u1*c0*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F41*F40r*u0*u1*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F40*F45r*u0*u5*c5*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F40*F45r*u0*u5*c5*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F40*F44r*u0*u4*c4*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F40*F44r*u0*u4*c4*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F40*F43r*u0*u3*c3*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F40*F43r*u0*u3*c3*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F40*F42r*u0*u2*c2*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F40*F42r*u0*u2*c2*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F40*F41r*u0*u1*c1*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F40*F41r*u0*u1*c1*f4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F40*F40r*u0**2*c0*f4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F40*F40r*u0**2*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F35*F25r*u5**2*c5*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F35*F25r*u5**2*c5*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F35*F24r*u4*u5*c4*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F35*F24r*u4*u5*c4*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F35*F23r*u3*u5*c3*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F35*F23r*u3*u5*c3*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F35*F22r*u2*u5*c2*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F35*F22r*u2*u5*c2*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F35*F21r*u1*u5*c1*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F35*F21r*u1*u5*c1*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F35*F20r*u0*u5*c0*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F35*F20r*u0*u5*c0*f2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F34*F25r*u4*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F34*F25r*u4*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F34*F24r*u4**2*c4*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F34*F24r*u4**2*c4*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F34*F23r*u3*u4*c3*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F34*F23r*u3*u4*c3*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F34*F22r*u2*u4*c2*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F34*F22r*u2*u4*c2*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F34*F21r*u1*u4*c1*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F34*F21r*u1*u4*c1*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F34*F20r*u0*u4*c0*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F34*F20r*u0*u4*c0*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F33*F25r*u3*u5*c5*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F33*F25r*u3*u5*c5*f2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F33*F24r*u3*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F33*F24r*u3*u4*c4*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F33*F23r*u3**2*c3*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F33*F23r*u3**2*c3*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F33*F22r*u2*u3*c2*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F33*F22r*u2*u3*c2*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F33*F21r*u1*u3*c1*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F33*F21r*u1*u3*c1*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F33*F20r*u0*u3*c0*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F33*F20r*u0*u3*c0*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F32*F25r*u2*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F32*F25r*u2*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F32*F24r*u2*u4*c4*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F32*F24r*u2*u4*c4*f2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F32*F23r*u2*u3*c3*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F32*F23r*u2*u3*c3*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F32*F22r*u2**2*c2*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F32*F22r*u2**2*c2*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F32*F21r*u1*u2*c1*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F32*F21r*u1*u2*c1*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F32*F20r*u0*u2*c0*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F32*F20r*u0*u2*c0*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F31*F25r*u1*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F31*F25r*u1*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F31*F24r*u1*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F31*F24r*u1*u4*c4*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F31*F23r*u1*u3*c3*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F31*F23r*u1*u3*c3*f2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F31*F22r*u1*u2*c2*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F31*F22r*u1*u2*c2*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F31*F21r*u1**2*c1*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F31*F21r*u1**2*c1*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F31*F20r*u0*u1*c0*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F31*F20r*u0*u1*c0*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F30*F25r*u0*u5*c5*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F30*F25r*u0*u5*c5*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F30*F24r*u0*u4*c4*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F30*F24r*u0*u4*c4*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F30*F23r*u0*u3*c3*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F30*F23r*u0*u3*c3*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F30*F22r*u0*u2*c2*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F30*F22r*u0*u2*c2*f2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F30*F21r*u0*u1*c1*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F30*F21r*u0*u1*c1*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F30*F20r*u0**2*c0*f2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F30*F20r*u0**2*c0*f2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F25*F35r*u5**2*c5*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F25*F35r*u5**2*c5*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F25*F34r*u4*u5*c4*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F25*F34r*u4*u5*c4*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F25*F33r*u3*u5*c3*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F25*F33r*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F25*F32r*u2*u5*c2*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F25*F32r*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F25*F31r*u1*u5*c1*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F25*F31r*u1*u5*c1*f3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F25*F30r*u0*u5*c0*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F25*F30r*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F24*F35r*u4*u5*c5*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F24*F35r*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F24*F34r*u4**2*c4*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F24*F34r*u4**2*c4*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F24*F33r*u3*u4*c3*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F24*F33r*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F24*F32r*u2*u4*c2*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F24*F32r*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F24*F31r*u1*u4*c1*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F24*F31r*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F24*F30r*u0*u4*c0*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F24*F30r*u0*u4*c0*f3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F23*F35r*u3*u5*c5*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F23*F35r*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F23*F34r*u3*u4*c4*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F23*F34r*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F23*F33r*u3**2*c3*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F23*F33r*u3**2*c3*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F23*F32r*u2*u3*c2*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F23*F32r*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F23*F31r*u1*u3*c1*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F23*F31r*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F23*F30r*u0*u3*c0*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F23*F30r*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F22*F35r*u2*u5*c5*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F22*F35r*u2*u5*c5*f3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F22*F34r*u2*u4*c4*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F22*F34r*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F22*F33r*u2*u3*c3*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F22*F33r*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F22*F32r*u2**2*c2*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F22*F32r*u2**2*c2*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F22*F31r*u1*u2*c1*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F22*F31r*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F22*F30r*u0*u2*c0*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F22*F30r*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F21*F35r*u1*u5*c5*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F21*F35r*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F21*F34r*u1*u4*c4*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F21*F34r*u1*u4*c4*f3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F21*F33r*u1*u3*c3*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F21*F33r*u1*u3*c3*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F21*F32r*u1*u2*c2*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F21*F32r*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F21*F31r*u1**2*c1*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F21*F31r*u1**2*c1*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F21*F30r*u0*u1*c0*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F21*F30r*u0*u1*c0*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F20*F35r*u0*u5*c5*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F20*F35r*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F20*F34r*u0*u4*c4*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F20*F34r*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F20*F33r*u0*u3*c3*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*svm*F20*F33r*u0*u3*c3*f3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F20*F32r*u0*u2*c2*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F20*F32r*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F20*F31r*u0*u1*c1*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F20*F31r*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*svm*F20*F30r*u0**2*c0*f3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*svm*F20*F30r*u0**2*c0*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F45*F35r*u5**2*c5*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F45*F34r*u4*u5*c4*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F45*F33r*u3*u5*c3*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F45*F32r*u2*u5*c2*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F45*F31r*u1*u5*c1*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F45*F30r*u0*u5*c0*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F44*F35r*u4*u5*c5*f3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*svp*F44*F34r*u4**2*c4*f3*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F44*F33r*u3*u4*c3*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F44*F32r*u2*u4*c2*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F44*F31r*u1*u4*c1*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F44*F30r*u0*u4*c0*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F43*F35r*u3*u5*c5*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F43*F34r*u3*u4*c4*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F43*F33r*u3**2*c3*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F43*F32r*u2*u3*c2*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F43*F31r*u1*u3*c1*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F43*F30r*u0*u3*c0*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F42*F35r*u2*u5*c5*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F42*F34r*u2*u4*c4*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F42*F33r*u2*u3*c3*f3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*svp*F42*F32r*u2**2*c2*f3*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F42*F31r*u1*u2*c1*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F42*F30r*u0*u2*c0*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F41*F35r*u1*u5*c5*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F41*F34r*u1*u4*c4*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F41*F33r*u1*u3*c3*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F41*F32r*u1*u2*c2*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F41*F31r*u1**2*c1*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F41*F30r*u0*u1*c0*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F40*F35r*u0*u5*c5*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F40*F34r*u0*u4*c4*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F40*F33r*u0*u3*c3*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F40*F32r*u0*u2*c2*f3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F40*F31r*u0*u1*c1*f3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*svp*F40*F30r*u0**2*c0*f3*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F35*F45r*u5**2*c5*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F35*F44r*u4*u5*c4*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F35*F43r*u3*u5*c3*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F35*F42r*u2*u5*c2*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F35*F41r*u1*u5*c1*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F35*F40r*u0*u5*c0*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F34*F45r*u4*u5*c5*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F34*F44r*u4**2*c4*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F34*F43r*u3*u4*c3*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F34*F42r*u2*u4*c2*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F34*F41r*u1*u4*c1*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F34*F40r*u0*u4*c0*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F33*F45r*u3*u5*c5*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*svp*F33*F44r*u3*u4*c4*f4*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F33*F43r*u3**2*c3*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F33*F42r*u2*u3*c2*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F33*F41r*u1*u3*c1*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F33*F40r*u0*u3*c0*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F32*F45r*u2*u5*c5*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F32*F44r*u2*u4*c4*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F32*F43r*u2*u3*c3*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F32*F42r*u2**2*c2*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F32*F41r*u1*u2*c1*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F32*F40r*u0*u2*c0*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F31*F45r*u1*u5*c5*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F31*F44r*u1*u4*c4*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F31*F43r*u1*u3*c3*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*svp*F31*F42r*u1*u2*c2*f4*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F31*F41r*u1**2*c1*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F31*F40r*u0*u1*c0*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F30*F45r*u0*u5*c5*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F30*F44r*u0*u4*c4*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F30*F43r*u0*u3*c3*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F30*F42r*u0*u2*c2*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F30*F41r*u0*u1*c1*f4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*svp*F30*F40r*u0**2*c0*f4*w1
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F45*F35r*u5**2*c5*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F45*F34r*u4*u5*c4*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F45*F33r*u3*u5*c3*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F45*F32r*u2*u5*c2*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F45*F31r*u1*u5*c1*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*svm*F45*F30r*u0*u5*c0*f3*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F44*F35r*u4*u5*c5*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F44*F34r*u4**2*c4*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F44*F33r*u3*u4*c3*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F44*F32r*u2*u4*c2*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F44*F31r*u1*u4*c1*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F44*F30r*u0*u4*c0*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F43*F35r*u3*u5*c5*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F43*F34r*u3*u4*c4*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F43*F33r*u3**2*c3*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F43*F32r*u2*u3*c2*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F43*F31r*u1*u3*c1*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F43*F30r*u0*u3*c0*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F42*F35r*u2*u5*c5*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*svm*F42*F34r*u2*u4*c4*f3*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F42*F33r*u2*u3*c3*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F42*F32r*u2**2*c2*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F42*F31r*u1*u2*c1*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F42*F30r*u0*u2*c0*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F41*F35r*u1*u5*c5*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F41*F34r*u1*u4*c4*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F41*F33r*u1*u3*c3*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F41*F32r*u1*u2*c2*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F41*F31r*u1**2*c1*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F41*F30r*u0*u1*c0*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F40*F35r*u0*u5*c5*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F40*F34r*u0*u4*c4*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F40*F33r*u0*u3*c3*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*svm*F40*F32r*u0*u2*c2*f3*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F40*F31r*u0*u1*c1*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F40*F30r*u0**2*c0*f3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F35*F45r*u5**2*c5*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F35*F44r*u4*u5*c4*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F35*F43r*u3*u5*c3*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F35*F42r*u2*u5*c2*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F35*F41r*u1*u5*c1*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F35*F40r*u0*u5*c0*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F34*F45r*u4*u5*c5*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F34*F44r*u4**2*c4*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F34*F43r*u3*u4*c3*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F34*F42r*u2*u4*c2*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F34*F41r*u1*u4*c1*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*svm*F34*F40r*u0*u4*c0*f4*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F33*F45r*u3*u5*c5*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F33*F44r*u3*u4*c4*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F33*F43r*u3**2*c3*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F33*F42r*u2*u3*c2*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F33*F41r*u1*u3*c1*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F33*F40r*u0*u3*c0*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F32*F45r*u2*u5*c5*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F32*F44r*u2*u4*c4*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F32*F43r*u2*u3*c3*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F32*F42r*u2**2*c2*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F32*F41r*u1*u2*c1*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F32*F40r*u0*u2*c0*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F31*F45r*u1*u5*c5*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*svm*F31*F44r*u1*u4*c4*f4*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F31*F43r*u1*u3*c3*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F31*F42r*u1*u2*c2*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F31*F41r*u1**2*c1*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F31*F40r*u0*u1*c0*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F30*F45r*u0*u5*c5*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F30*F44r*u0*u4*c4*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F30*F43r*u0*u3*c3*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F30*F42r*u0*u2*c2*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F30*F41r*u0*u1*c1*f4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*svm*F30*F40r*u0**2*c0*f4*w2
     &  + 4.D0*PC_q*i_*svp*svm*F15*F15i*u5**2*c5*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F15*F15i*u5**2*c5*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F15*F14i*u4*u5*c4*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*svm*F15*F14i*u4*u5*c4*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F15*F13i*u3*u5*c3*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F15*F13i*u3*u5*c3*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F15*F12i*u2*u5*c2*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F15*F12i*u2*u5*c2*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F15*F11i*u1*u5*c1*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F15*F11i*u1*u5*c1*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F15*F10i*u0*u5*c0*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F15*F10i*u0*u5*c0*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F14*F15i*u4*u5*c5*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F14*F15i*u4*u5*c5*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F14*F14i*u4**2*c4*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F14*F14i*u4**2*c4*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F14*F13i*u3*u4*c3*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F14*F13i*u3*u4*c3*f1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svp*svm*F14*F12i*u2*u4*c2*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F14*F12i*u2*u4*c2*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F14*F11i*u1*u4*c1*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F14*F11i*u1*u4*c1*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F14*F10i*u0*u4*c0*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F14*F10i*u0*u4*c0*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F13*F15i*u3*u5*c5*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F13*F15i*u3*u5*c5*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F13*F14i*u3*u4*c4*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F13*F14i*u3*u4*c4*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F13*F13i*u3**2*c3*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F13*F13i*u3**2*c3*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F13*F12i*u2*u3*c2*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F13*F12i*u2*u3*c2*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F13*F11i*u1*u3*c1*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*svm*F13*F11i*u1*u3*c1*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F13*F10i*u0*u3*c0*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F13*F10i*u0*u3*c0*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F12*F15i*u2*u5*c5*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F12*F15i*u2*u5*c5*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F12*F14i*u2*u4*c4*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F12*F14i*u2*u4*c4*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F12*F13i*u2*u3*c3*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F12*F13i*u2*u3*c3*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F12*F12i*u2**2*c2*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F12*F12i*u2**2*c2*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F12*F11i*u1*u2*c1*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F12*F11i*u1*u2*c1*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F12*F10i*u0*u2*c0*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F12*F10i*u0*u2*c0*f1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svp*svm*F11*F15i*u1*u5*c5*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F11*F15i*u1*u5*c5*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F11*F14i*u1*u4*c4*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F11*F14i*u1*u4*c4*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F11*F13i*u1*u3*c3*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F11*F13i*u1*u3*c3*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F11*F12i*u1*u2*c2*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F11*F12i*u1*u2*c2*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F11*F11i*u1**2*c1*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F11*F11i*u1**2*c1*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F11*F10i*u0*u1*c0*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F11*F10i*u0*u1*c0*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F10*F15i*u0*u5*c5*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F10*F15i*u0*u5*c5*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F10*F14i*u0*u4*c4*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*svm*F10*F14i*u0*u4*c4*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F10*F13i*u0*u3*c3*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F10*F13i*u0*u3*c3*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F10*F12i*u0*u2*c2*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F10*F12i*u0*u2*c2*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F10*F11i*u0*u1*c1*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F10*F11i*u0*u1*c1*f1*w1
     &  + 4.D0*PC_q*i_*svp*svm*F10*F10i*u0**2*c0*f1*w2
     &  - 4.D0*PC_q*i_*svp*svm*F10*F10i*u0**2*c0*f1*w1
     &  - 4.D0*PC_q*i_*ssm*svp*F25*F15i*u5**2*c5*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F25*F14i*u4*u5*c4*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F25*F13i*u3*u5*c3*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F25*F12i*u2*u5*c2*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F25*F11i*u1*u5*c1*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F25*F10i*u0*u5*c0*f1
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*svp*F24*F15i*u4*u5*c5*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F24*F14i*u4**2*c4*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F24*F13i*u3*u4*c3*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F24*F12i*u2*u4*c2*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F24*F11i*u1*u4*c1*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F24*F10i*u0*u4*c0*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F23*F15i*u3*u5*c5*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F23*F14i*u3*u4*c4*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F23*F13i*u3**2*c3*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F23*F12i*u2*u3*c2*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F23*F11i*u1*u3*c1*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F23*F10i*u0*u3*c0*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F22*F15i*u2*u5*c5*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F22*F14i*u2*u4*c4*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F22*F13i*u2*u3*c3*f1
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*svp*F22*F12i*u2**2*c2*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F22*F11i*u1*u2*c1*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F22*F10i*u0*u2*c0*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F21*F15i*u1*u5*c5*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F21*F14i*u1*u4*c4*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F21*F13i*u1*u3*c3*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F21*F12i*u1*u2*c2*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F21*F11i*u1**2*c1*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F21*F10i*u0*u1*c0*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F20*F15i*u0*u5*c5*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F20*F14i*u0*u4*c4*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F20*F13i*u0*u3*c3*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F20*F12i*u0*u2*c2*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F20*F11i*u0*u1*c1*f1
     &  - 4.D0*PC_q*i_*ssm*svp*F20*F10i*u0**2*c0*f1
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*svp*F15*F25i*u5**2*c5*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F15*F24i*u4*u5*c4*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F15*F23i*u3*u5*c3*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F15*F22i*u2*u5*c2*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F15*F21i*u1*u5*c1*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F15*F20i*u0*u5*c0*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F14*F25i*u4*u5*c5*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F14*F24i*u4**2*c4*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F14*F23i*u3*u4*c3*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F14*F22i*u2*u4*c2*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F14*F21i*u1*u4*c1*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F14*F20i*u0*u4*c0*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F13*F25i*u3*u5*c5*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F13*F24i*u3*u4*c4*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F13*F23i*u3**2*c3*f2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*svp*F13*F22i*u2*u3*c2*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F13*F21i*u1*u3*c1*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F13*F20i*u0*u3*c0*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F12*F25i*u2*u5*c5*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F12*F24i*u2*u4*c4*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F12*F23i*u2*u3*c3*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F12*F22i*u2**2*c2*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F12*F21i*u1*u2*c1*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F12*F20i*u0*u2*c0*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F11*F25i*u1*u5*c5*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F11*F24i*u1*u4*c4*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F11*F23i*u1*u3*c3*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F11*F22i*u1*u2*c2*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F11*F21i*u1**2*c1*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F11*F20i*u0*u1*c0*f2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*svp*F10*F25i*u0*u5*c5*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F10*F24i*u0*u4*c4*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F10*F23i*u0*u3*c3*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F10*F22i*u0*u2*c2*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F10*F21i*u0*u1*c1*f2
     &  - 4.D0*PC_q*i_*ssm*svp*F10*F20i*u0**2*c0*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F25*F15i*u5**2*c5*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F25*F14i*u4*u5*c4*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F25*F13i*u3*u5*c3*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F25*F12i*u2*u5*c2*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F25*F11i*u1*u5*c1*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F25*F10i*u0*u5*c0*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F24*F15i*u4*u5*c5*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F24*F14i*u4**2*c4*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F24*F13i*u3*u4*c3*f1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*svm*F24*F12i*u2*u4*c2*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F24*F11i*u1*u4*c1*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F24*F10i*u0*u4*c0*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F23*F15i*u3*u5*c5*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F23*F14i*u3*u4*c4*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F23*F13i*u3**2*c3*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F23*F12i*u2*u3*c2*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F23*F11i*u1*u3*c1*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F23*F10i*u0*u3*c0*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F22*F15i*u2*u5*c5*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F22*F14i*u2*u4*c4*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F22*F13i*u2*u3*c3*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F22*F12i*u2**2*c2*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F22*F11i*u1*u2*c1*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F22*F10i*u0*u2*c0*f1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*svm*F21*F15i*u1*u5*c5*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F21*F14i*u1*u4*c4*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F21*F13i*u1*u3*c3*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F21*F12i*u1*u2*c2*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F21*F11i*u1**2*c1*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F21*F10i*u0*u1*c0*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F20*F15i*u0*u5*c5*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F20*F14i*u0*u4*c4*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F20*F13i*u0*u3*c3*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F20*F12i*u0*u2*c2*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F20*F11i*u0*u1*c1*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F20*F10i*u0**2*c0*f1
     &  + 4.D0*PC_q*i_*ssp*svm*F15*F25i*u5**2*c5*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F15*F24i*u4*u5*c4*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F15*F23i*u3*u5*c3*f2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*svm*F15*F22i*u2*u5*c2*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F15*F21i*u1*u5*c1*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F15*F20i*u0*u5*c0*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F14*F25i*u4*u5*c5*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F14*F24i*u4**2*c4*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F14*F23i*u3*u4*c3*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F14*F22i*u2*u4*c2*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F14*F21i*u1*u4*c1*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F14*F20i*u0*u4*c0*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F13*F25i*u3*u5*c5*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F13*F24i*u3*u4*c4*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F13*F23i*u3**2*c3*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F13*F22i*u2*u3*c2*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F13*F21i*u1*u3*c1*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F13*F20i*u0*u3*c0*f2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*svm*F12*F25i*u2*u5*c5*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F12*F24i*u2*u4*c4*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F12*F23i*u2*u3*c3*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F12*F22i*u2**2*c2*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F12*F21i*u1*u2*c1*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F12*F20i*u0*u2*c0*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F11*F25i*u1*u5*c5*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F11*F24i*u1*u4*c4*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F11*F23i*u1*u3*c3*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F11*F22i*u1*u2*c2*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F11*F21i*u1**2*c1*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F11*F20i*u0*u1*c0*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F10*F25i*u0*u5*c5*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F10*F24i*u0*u4*c4*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F10*F23i*u0*u3*c3*f2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*svm*F10*F22i*u0*u2*c2*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F10*F21i*u0*u1*c1*f2
     &  + 4.D0*PC_q*i_*ssp*svm*F10*F20i*u0**2*c0*f2
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F35*F15i*u5**2*c5*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F35*F14i*u4*u5*c4*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F35*F13i*u3*u5*c3*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F35*F12i*u2*u5*c2*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F35*F11i*u1*u5*c1*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F35*F10i*u0*u5*c0*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F34*F15i*u4*u5*c5*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F34*F14i*u4**2*c4*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F34*F13i*u3*u4*c3*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F34*F12i*u2*u4*c2*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F34*F11i*u1*u4*c1*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F34*F10i*u0*u4*c0*f1
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*svp*F33*F15i*u3*u5*c5*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F33*F14i*u3*u4*c4*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F33*F13i*u3**2*c3*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F33*F12i*u2*u3*c2*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F33*F11i*u1*u3*c1*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F33*F10i*u0*u3*c0*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F32*F15i*u2*u5*c5*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F32*F14i*u2*u4*c4*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F32*F13i*u2*u3*c3*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F32*F12i*u2**2*c2*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F32*F11i*u1*u2*c1*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F32*F10i*u0*u2*c0*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F31*F15i*u1*u5*c5*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F31*F14i*u1*u4*c4*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F31*F13i*u1*u3*c3*f1
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*svp*F31*F12i*u1*u2*c2*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F31*F11i*u1**2*c1*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F31*F10i*u0*u1*c0*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F30*F15i*u0*u5*c5*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F30*F14i*u0*u4*c4*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F30*F13i*u0*u3*c3*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F30*F12i*u0*u2*c2*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F30*F11i*u0*u1*c1*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F30*F10i*u0**2*c0*f1
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F15*F35i*u5**2*c5*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F15*F34i*u4*u5*c4*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F15*F33i*u3*u5*c3*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F15*F32i*u2*u5*c2*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F15*F31i*u1*u5*c1*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F15*F30i*u0*u5*c0*f3
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*svp*F14*F35i*u4*u5*c5*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F14*F34i*u4**2*c4*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F14*F33i*u3*u4*c3*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F14*F32i*u2*u4*c2*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F14*F31i*u1*u4*c1*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F14*F30i*u0*u4*c0*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F13*F35i*u3*u5*c5*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F13*F34i*u3*u4*c4*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F13*F33i*u3**2*c3*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F13*F32i*u2*u3*c2*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F13*F31i*u1*u3*c1*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F13*F30i*u0*u3*c0*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F12*F35i*u2*u5*c5*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F12*F34i*u2*u4*c4*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F12*F33i*u2*u3*c3*f3
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*svp*F12*F32i*u2**2*c2*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F12*F31i*u1*u2*c1*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F12*F30i*u0*u2*c0*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F11*F35i*u1*u5*c5*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F11*F34i*u1*u4*c4*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F11*F33i*u1*u3*c3*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F11*F32i*u1*u2*c2*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F11*F31i*u1**2*c1*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F11*F30i*u0*u1*c0*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F10*F35i*u0*u5*c5*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F10*F34i*u0*u4*c4*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F10*F33i*u0*u3*c3*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F10*F32i*u0*u2*c2*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F10*F31i*u0*u1*c1*f3
     &  - 4.D0*PC_q*i_*qq*ssm*svp*F10*F30i*u0**2*c0*f3
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*svm*F35*F15i*u5**2*c5*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F35*F14i*u4*u5*c4*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F35*F13i*u3*u5*c3*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F35*F12i*u2*u5*c2*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F35*F11i*u1*u5*c1*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F35*F10i*u0*u5*c0*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F34*F15i*u4*u5*c5*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F34*F14i*u4**2*c4*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F34*F13i*u3*u4*c3*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F34*F12i*u2*u4*c2*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F34*F11i*u1*u4*c1*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F34*F10i*u0*u4*c0*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F33*F15i*u3*u5*c5*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F33*F14i*u3*u4*c4*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F33*F13i*u3**2*c3*f1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*svm*F33*F12i*u2*u3*c2*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F33*F11i*u1*u3*c1*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F33*F10i*u0*u3*c0*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F32*F15i*u2*u5*c5*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F32*F14i*u2*u4*c4*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F32*F13i*u2*u3*c3*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F32*F12i*u2**2*c2*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F32*F11i*u1*u2*c1*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F32*F10i*u0*u2*c0*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F31*F15i*u1*u5*c5*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F31*F14i*u1*u4*c4*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F31*F13i*u1*u3*c3*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F31*F12i*u1*u2*c2*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F31*F11i*u1**2*c1*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F31*F10i*u0*u1*c0*f1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*svm*F30*F15i*u0*u5*c5*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F30*F14i*u0*u4*c4*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F30*F13i*u0*u3*c3*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F30*F12i*u0*u2*c2*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F30*F11i*u0*u1*c1*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F30*F10i*u0**2*c0*f1
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F15*F35i*u5**2*c5*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F15*F34i*u4*u5*c4*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F15*F33i*u3*u5*c3*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F15*F32i*u2*u5*c2*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F15*F31i*u1*u5*c1*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F15*F30i*u0*u5*c0*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F14*F35i*u4*u5*c5*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F14*F34i*u4**2*c4*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F14*F33i*u3*u4*c3*f3
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*svm*F14*F32i*u2*u4*c2*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F14*F31i*u1*u4*c1*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F14*F30i*u0*u4*c0*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F13*F35i*u3*u5*c5*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F13*F34i*u3*u4*c4*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F13*F33i*u3**2*c3*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F13*F32i*u2*u3*c2*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F13*F31i*u1*u3*c1*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F13*F30i*u0*u3*c0*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F12*F35i*u2*u5*c5*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F12*F34i*u2*u4*c4*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F12*F33i*u2*u3*c3*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F12*F32i*u2**2*c2*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F12*F31i*u1*u2*c1*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F12*F30i*u0*u2*c0*f3
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*svm*F11*F35i*u1*u5*c5*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F11*F34i*u1*u4*c4*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F11*F33i*u1*u3*c3*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F11*F32i*u1*u2*c2*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F11*F31i*u1**2*c1*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F11*F30i*u0*u1*c0*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F10*F35i*u0*u5*c5*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F10*F34i*u0*u4*c4*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F10*F33i*u0*u3*c3*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F10*F32i*u0*u2*c2*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F10*F31i*u0*u1*c1*f3
     &  + 4.D0*PC_q*i_*qq*ssp*svm*F10*F30i*u0**2*c0*f3
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F25*F25i*u5**2*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F25*F25i*u5**2*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F25*F24i*u4*u5*c4*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*svm*F25*F24i*u4*u5*c4*f2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F25*F23i*u3*u5*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F25*F23i*u3*u5*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F25*F22i*u2*u5*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F25*F22i*u2*u5*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F25*F21i*u1*u5*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F25*F21i*u1*u5*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F25*F20i*u0*u5*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F25*F20i*u0*u5*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F24*F25i*u4*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F24*F25i*u4*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F24*F24i*u4**2*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F24*F24i*u4**2*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F24*F23i*u3*u4*c3*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*svm*F24*F23i*u3*u4*c3*f2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F24*F22i*u2*u4*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F24*F22i*u2*u4*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F24*F21i*u1*u4*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F24*F21i*u1*u4*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F24*F20i*u0*u4*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F24*F20i*u0*u4*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F23*F25i*u3*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F23*F25i*u3*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F23*F24i*u3*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F23*F24i*u3*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F23*F23i*u3**2*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F23*F23i*u3**2*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F23*F22i*u2*u3*c2*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*svm*F23*F22i*u2*u3*c2*f2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F23*F21i*u1*u3*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F23*F21i*u1*u3*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F23*F20i*u0*u3*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F23*F20i*u0*u3*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F22*F25i*u2*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F22*F25i*u2*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F22*F24i*u2*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F22*F24i*u2*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F22*F23i*u2*u3*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F22*F23i*u2*u3*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F22*F22i*u2**2*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F22*F22i*u2**2*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F22*F21i*u1*u2*c1*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*svm*F22*F21i*u1*u2*c1*f2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F22*F20i*u0*u2*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F22*F20i*u0*u2*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F21*F25i*u1*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F21*F25i*u1*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F21*F24i*u1*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F21*F24i*u1*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F21*F23i*u1*u3*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F21*F23i*u1*u3*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F21*F22i*u1*u2*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F21*F22i*u1*u2*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F21*F21i*u1**2*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F21*F21i*u1**2*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F21*F20i*u0*u1*c0*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*svm*F21*F20i*u0*u1*c0*f2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F20*F25i*u0*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F20*F25i*u0*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F20*F24i*u0*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F20*F24i*u0*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F20*F23i*u0*u3*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F20*F23i*u0*u3*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F20*F22i*u0*u2*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F20*F22i*u0*u2*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F20*F21i*u0*u1*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F20*F21i*u0*u1*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*svm*F20*F20i*u0**2*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*svm*F20*F20i*u0**2*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F45i*u5**2*c5*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F45i*u5**2*c5*
     & f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F44i*u4*u5*c4*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F44i*u4*u5*c4*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F43i*u3*u5*c3*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F43i*u3*u5*c3*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F42i*u2*u5*c2*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F42i*u2*u5*c2*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F41i*u1*u5*c1*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F41i*u1*u5*c1*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F40i*u0*u5*c0*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F45*F40i*u0*u5*c0*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F45i*u4*u5*c5*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F45i*u4*u5*c5*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F44i*u4**2*c4*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F44i*u4**2*c4*
     & f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F43i*u3*u4*c3*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F43i*u3*u4*c3*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F42i*u2*u4*c2*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F42i*u2*u4*c2*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F41i*u1*u4*c1*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F41i*u1*u4*c1*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F40i*u0*u4*c0*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F44*F40i*u0*u4*c0*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F45i*u3*u5*c5*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F45i*u3*u5*c5*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F44i*u3*u4*c4*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F44i*u3*u4*c4*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F43i*u3**2*c3*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F43i*u3**2*c3*
     & f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F42i*u2*u3*c2*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F42i*u2*u3*c2*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F41i*u1*u3*c1*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F41i*u1*u3*c1*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F40i*u0*u3*c0*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F43*F40i*u0*u3*c0*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F45i*u2*u5*c5*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F45i*u2*u5*c5*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F44i*u2*u4*c4*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F44i*u2*u4*c4*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F43i*u2*u3*c3*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F43i*u2*u3*c3*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F42i*u2**2*c2*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F42i*u2**2*c2*
     & f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F41i*u1*u2*c1*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F41i*u1*u2*c1*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F40i*u0*u2*c0*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F42*F40i*u0*u2*c0*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F45i*u1*u5*c5*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F45i*u1*u5*c5*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F44i*u1*u4*c4*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F44i*u1*u4*c4*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F43i*u1*u3*c3*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F43i*u1*u3*c3*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F42i*u1*u2*c2*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F42i*u1*u2*c2*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F41i*u1**2*c1*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F41i*u1**2*c1*
     & f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F40i*u0*u1*c0*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F41*F40i*u0*u1*c0*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F45i*u0*u5*c5*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F45i*u0*u5*c5*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F44i*u0*u4*c4*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F44i*u0*u4*c4*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F43i*u0*u3*c3*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F43i*u0*u3*c3*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F42i*u0*u2*c2*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F42i*u0*u2*c2*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F41i*u0*u1*c1*f4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F41i*u0*u1*c1*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F40i*u0**2*c0*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F40*F40i*u0**2*c0*
     & f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F25i*u5**2*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F25i*u5**2*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F24i*u4*u5*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F24i*u4*u5*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F23i*u3*u5*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F23i*u3*u5*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F22i*u2*u5*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F22i*u2*u5*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F21i*u1*u5*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F21i*u1*u5*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F20i*u0*u5*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F35*F20i*u0*u5*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F25i*u4*u5*c5*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F25i*u4*u5*c5*
     & f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F24i*u4**2*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F24i*u4**2*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F23i*u3*u4*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F23i*u3*u4*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F22i*u2*u4*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F22i*u2*u4*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F21i*u1*u4*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F21i*u1*u4*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F20i*u0*u4*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F34*F20i*u0*u4*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F25i*u3*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F25i*u3*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F24i*u3*u4*c4*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F24i*u3*u4*c4*
     & f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F23i*u3**2*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F23i*u3**2*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F22i*u2*u3*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F22i*u2*u3*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F21i*u1*u3*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F21i*u1*u3*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F20i*u0*u3*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F33*F20i*u0*u3*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F25i*u2*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F25i*u2*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F24i*u2*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F24i*u2*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F23i*u2*u3*c3*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F23i*u2*u3*c3*
     & f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F22i*u2**2*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F22i*u2**2*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F21i*u1*u2*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F21i*u1*u2*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F20i*u0*u2*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F32*F20i*u0*u2*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F25i*u1*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F25i*u1*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F24i*u1*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F24i*u1*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F23i*u1*u3*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F23i*u1*u3*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F22i*u1*u2*c2*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F22i*u1*u2*c2*
     & f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F21i*u1**2*c1*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F21i*u1**2*c1*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F20i*u0*u1*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F31*F20i*u0*u1*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F25i*u0*u5*c5*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F25i*u0*u5*c5*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F24i*u0*u4*c4*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F24i*u0*u4*c4*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F23i*u0*u3*c3*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F23i*u0*u3*c3*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F22i*u0*u2*c2*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F22i*u0*u2*c2*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F21i*u0*u1*c1*f2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F21i*u0*u1*c1*
     & f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F20i*u0**2*c0*f2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F30*F20i*u0**2*c0*f2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F35i*u5**2*c5*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F35i*u5**2*c5*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F34i*u4*u5*c4*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F34i*u4*u5*c4*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F33i*u3*u5*c3*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F33i*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F32i*u2*u5*c2*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F32i*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F31i*u1*u5*c1*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F31i*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F30i*u0*u5*c0*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F25*F30i*u0*u5*c0*
     & f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F35i*u4*u5*c5*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F35i*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F34i*u4**2*c4*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F34i*u4**2*c4*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F33i*u3*u4*c3*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F33i*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F32i*u2*u4*c2*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F32i*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F31i*u1*u4*c1*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F31i*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F30i*u0*u4*c0*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F24*F30i*u0*u4*c0*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F35i*u3*u5*c5*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F35i*u3*u5*c5*
     & f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F34i*u3*u4*c4*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F34i*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F33i*u3**2*c3*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F33i*u3**2*c3*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F32i*u2*u3*c2*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F32i*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F31i*u1*u3*c1*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F31i*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F30i*u0*u3*c0*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F23*F30i*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F35i*u2*u5*c5*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F35i*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F34i*u2*u4*c4*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F34i*u2*u4*c4*
     & f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F33i*u2*u3*c3*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F33i*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F32i*u2**2*c2*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F32i*u2**2*c2*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F31i*u1*u2*c1*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F31i*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F30i*u0*u2*c0*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F22*F30i*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F35i*u1*u5*c5*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F35i*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F34i*u1*u4*c4*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F34i*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F33i*u1*u3*c3*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F33i*u1*u3*c3*
     & f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F32i*u1*u2*c2*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F32i*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F31i*u1**2*c1*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F31i*u1**2*c1*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F30i*u0*u1*c0*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F21*F30i*u0*u1*c0*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F35i*u0*u5*c5*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F35i*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F34i*u0*u4*c4*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F34i*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F33i*u0*u3*c3*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F33i*u0*u3*c3*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F32i*u0*u2*c2*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F32i*u0*u2*c2*
     & f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F31i*u0*u1*c1*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F31i*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F30i*u0**2*c0*f3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*svm*F20*F30i*u0**2*c0*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F45*F35i*u5**2*c5*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F45*F34i*u4*u5*c4*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F45*F33i*u3*u5*c3*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F45*F32i*u2*u5*c2*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F45*F31i*u1*u5*c1*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F45*F30i*u0*u5*c0*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F44*F35i*u4*u5*c5*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F44*F34i*u4**2*c4*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F44*F33i*u3*u4*c3*f3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F44*F32i*u2*u4*c2*
     & f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F44*F31i*u1*u4*c1*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F44*F30i*u0*u4*c0*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F43*F35i*u3*u5*c5*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F43*F34i*u3*u4*c4*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F43*F33i*u3**2*c3*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F43*F32i*u2*u3*c2*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F43*F31i*u1*u3*c1*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F43*F30i*u0*u3*c0*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F42*F35i*u2*u5*c5*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F42*F34i*u2*u4*c4*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F42*F33i*u2*u3*c3*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F42*F32i*u2**2*c2*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F42*F31i*u1*u2*c1*f3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F42*F30i*u0*u2*c0*
     & f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F41*F35i*u1*u5*c5*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F41*F34i*u1*u4*c4*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F41*F33i*u1*u3*c3*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F41*F32i*u1*u2*c2*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F41*F31i*u1**2*c1*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F41*F30i*u0*u1*c0*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F40*F35i*u0*u5*c5*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F40*F34i*u0*u4*c4*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F40*F33i*u0*u3*c3*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F40*F32i*u0*u2*c2*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F40*F31i*u0*u1*c1*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F40*F30i*u0**2*c0*f3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F35*F45i*u5**2*c5*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F35*F44i*u4*u5*c4*
     & f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F35*F43i*u3*u5*c3*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F35*F42i*u2*u5*c2*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F35*F41i*u1*u5*c1*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F35*F40i*u0*u5*c0*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F34*F45i*u4*u5*c5*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F34*F44i*u4**2*c4*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F34*F43i*u3*u4*c3*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F34*F42i*u2*u4*c2*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F34*F41i*u1*u4*c1*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F34*F40i*u0*u4*c0*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F33*F45i*u3*u5*c5*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F33*F44i*u3*u4*c4*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F33*F43i*u3**2*c3*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F33*F42i*u2*u3*c2*
     & f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F33*F41i*u1*u3*c1*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F33*F40i*u0*u3*c0*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F32*F45i*u2*u5*c5*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F32*F44i*u2*u4*c4*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F32*F43i*u2*u3*c3*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F32*F42i*u2**2*c2*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F32*F41i*u1*u2*c1*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F32*F40i*u0*u2*c0*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F31*F45i*u1*u5*c5*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F31*F44i*u1*u4*c4*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F31*F43i*u1*u3*c3*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F31*F42i*u1*u2*c2*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F31*F41i*u1**2*c1*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F31*F40i*u0*u1*c0*
     & f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F30*F45i*u0*u5*c5*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F30*F44i*u0*u4*c4*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F30*F43i*u0*u3*c3*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F30*F42i*u0*u2*c2*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F30*F41i*u0*u1*c1*f4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*svp*F30*F40i*u0**2*c0*f4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F45*F35i*u5**2*c5*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F45*F34i*u4*u5*c4*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F45*F33i*u3*u5*c3*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F45*F32i*u2*u5*c2*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F45*F31i*u1*u5*c1*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F45*F30i*u0*u5*c0*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F44*F35i*u4*u5*c5*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F44*F34i*u4**2*c4*
     & f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F44*F33i*u3*u4*c3*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F44*F32i*u2*u4*c2*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F44*F31i*u1*u4*c1*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F44*F30i*u0*u4*c0*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F43*F35i*u3*u5*c5*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F43*F34i*u3*u4*c4*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F43*F33i*u3**2*c3*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F43*F32i*u2*u3*c2*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F43*F31i*u1*u3*c1*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F43*F30i*u0*u3*c0*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F42*F35i*u2*u5*c5*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F42*F34i*u2*u4*c4*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F42*F33i*u2*u3*c3*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F42*F32i*u2**2*c2*
     & f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F42*F31i*u1*u2*c1*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F42*F30i*u0*u2*c0*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F41*F35i*u1*u5*c5*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F41*F34i*u1*u4*c4*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F41*F33i*u1*u3*c3*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F41*F32i*u1*u2*c2*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F41*F31i*u1**2*c1*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F41*F30i*u0*u1*c0*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F40*F35i*u0*u5*c5*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F40*F34i*u0*u4*c4*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F40*F33i*u0*u3*c3*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F40*F32i*u0*u2*c2*f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F40*F31i*u0*u1*c1*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F40*F30i*u0**2*c0*
     & f3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F35*F45i*u5**2*c5*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F35*F44i*u4*u5*c4*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F35*F43i*u3*u5*c3*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F35*F42i*u2*u5*c2*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F35*F41i*u1*u5*c1*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F35*F40i*u0*u5*c0*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F34*F45i*u4*u5*c5*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F34*F44i*u4**2*c4*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F34*F43i*u3*u4*c3*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F34*F42i*u2*u4*c2*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F34*F41i*u1*u4*c1*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F34*F40i*u0*u4*c0*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F33*F45i*u3*u5*c5*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F33*F44i*u3*u4*c4*
     & f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F33*F43i*u3**2*c3*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F33*F42i*u2*u3*c2*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F33*F41i*u1*u3*c1*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F33*F40i*u0*u3*c0*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F32*F45i*u2*u5*c5*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F32*F44i*u2*u4*c4*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F32*F43i*u2*u3*c3*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F32*F42i*u2**2*c2*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F32*F41i*u1*u2*c1*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F32*F40i*u0*u2*c0*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F31*F45i*u1*u5*c5*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F31*F44i*u1*u4*c4*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F31*F43i*u1*u3*c3*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F31*F42i*u1*u2*c2*
     & f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F31*F41i*u1**2*c1*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F31*F40i*u0*u1*c0*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F30*F45i*u0*u5*c5*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F30*F44i*u0*u4*c4*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F30*F43i*u0*u3*c3*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F30*F42i*u0*u2*c2*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F30*F41i*u0*u1*c1*f4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*svm*F30*F40i*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F45*F15r*u5**2*c5*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F45*F15r*u5**2*c5*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F45*F14r*u4*u5*c4*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F45*F14r*u4*u5*c4*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F45*F13r*u3*u5*c3*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F45*F13r*u3*u5*c3*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F45*F12r*u2*u5*c2*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F45*F12r*u2*u5*c2*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F45*F11r*u1*u5*c1*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F45*F11r*u1*u5*c1*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F45*F10r*u0*u5*c0*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F45*F10r*u0*u5*c0*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F44*F15r*u4*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F44*F15r*u4*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F44*F14r*u4**2*c4*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F44*F14r*u4**2*c4*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F44*F13r*u3*u4*c3*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F44*F13r*u3*u4*c3*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F44*F12r*u2*u4*c2*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F44*F12r*u2*u4*c2*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F44*F11r*u1*u4*c1*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F44*F11r*u1*u4*c1*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F44*F10r*u0*u4*c0*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F44*F10r*u0*u4*c0*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F43*F15r*u3*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F43*F15r*u3*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F43*F14r*u3*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F43*F14r*u3*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F43*F13r*u3**2*c3*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F43*F13r*u3**2*c3*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F43*F12r*u2*u3*c2*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F43*F12r*u2*u3*c2*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F43*F11r*u1*u3*c1*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F43*F11r*u1*u3*c1*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F43*F10r*u0*u3*c0*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F43*F10r*u0*u3*c0*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F42*F15r*u2*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F42*F15r*u2*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F42*F14r*u2*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F42*F14r*u2*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F42*F13r*u2*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F42*F13r*u2*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F42*F12r*u2**2*c2*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F42*F12r*u2**2*c2*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F42*F11r*u1*u2*c1*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F42*F11r*u1*u2*c1*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F42*F10r*u0*u2*c0*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F42*F10r*u0*u2*c0*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F41*F15r*u1*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F41*F15r*u1*u5*c5*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F41*F14r*u1*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F41*F14r*u1*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F41*F13r*u1*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F41*F13r*u1*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F41*F12r*u1*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F41*F12r*u1*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F41*F11r*u1**2*c1*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F41*F11r*u1**2*c1*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F41*F10r*u0*u1*c0*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F41*F10r*u0*u1*c0*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F40*F15r*u0*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F40*F15r*u0*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F40*F14r*u0*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F40*F14r*u0*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F40*F13r*u0*u3*c3*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F40*F13r*u0*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F40*F12r*u0*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F40*F12r*u0*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F40*F11r*u0*u1*c1*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F40*F11r*u0*u1*c1*f1*w1
     &  - 4.D0*PC_q**2*svp*svm*F40*F10r*u0**2*c0*f1*w2
     &  - 4.D0*PC_q**2*svp*svm*F40*F10r*u0**2*c0*f1*w1
     &  + 8.D0*PC_q**2*svp*svm*F25*F25r*u5**2*c5*f2
     &  + 8.D0*PC_q**2*svp*svm*F25*F24r*u4*u5*c4*f2
     &  + 8.D0*PC_q**2*svp*svm*F25*F23r*u3*u5*c3*f2
     &  + 8.D0*PC_q**2*svp*svm*F25*F22r*u2*u5*c2*f2
     &  + 8.D0*PC_q**2*svp*svm*F25*F21r*u1*u5*c1*f2
     &  + 8.D0*PC_q**2*svp*svm*F25*F20r*u0*u5*c0*f2
     &  + 8.D0*PC_q**2*svp*svm*F24*F25r*u4*u5*c5*f2
     &  + 8.D0*PC_q**2*svp*svm*F24*F24r*u4**2*c4*f2
     &
      traza1 = traza1 + 8.D0*PC_q**2*svp*svm*F24*F23r*u3*u4*c3*f2
     &  + 8.D0*PC_q**2*svp*svm*F24*F22r*u2*u4*c2*f2
     &  + 8.D0*PC_q**2*svp*svm*F24*F21r*u1*u4*c1*f2
     &  + 8.D0*PC_q**2*svp*svm*F24*F20r*u0*u4*c0*f2
     &  + 8.D0*PC_q**2*svp*svm*F23*F25r*u3*u5*c5*f2
     &  + 8.D0*PC_q**2*svp*svm*F23*F24r*u3*u4*c4*f2
     &  + 8.D0*PC_q**2*svp*svm*F23*F23r*u3**2*c3*f2
     &  + 8.D0*PC_q**2*svp*svm*F23*F22r*u2*u3*c2*f2
     &  + 8.D0*PC_q**2*svp*svm*F23*F21r*u1*u3*c1*f2
     &  + 8.D0*PC_q**2*svp*svm*F23*F20r*u0*u3*c0*f2
     &  + 8.D0*PC_q**2*svp*svm*F22*F25r*u2*u5*c5*f2
     &  + 8.D0*PC_q**2*svp*svm*F22*F24r*u2*u4*c4*f2
     &  + 8.D0*PC_q**2*svp*svm*F22*F23r*u2*u3*c3*f2
     &  + 8.D0*PC_q**2*svp*svm*F22*F22r*u2**2*c2*f2
     &  + 8.D0*PC_q**2*svp*svm*F22*F21r*u1*u2*c1*f2
     &
      traza1 = traza1 + 8.D0*PC_q**2*svp*svm*F22*F20r*u0*u2*c0*f2
     &  + 8.D0*PC_q**2*svp*svm*F21*F25r*u1*u5*c5*f2
     &  + 8.D0*PC_q**2*svp*svm*F21*F24r*u1*u4*c4*f2
     &  + 8.D0*PC_q**2*svp*svm*F21*F23r*u1*u3*c3*f2
     &  + 8.D0*PC_q**2*svp*svm*F21*F22r*u1*u2*c2*f2
     &  + 8.D0*PC_q**2*svp*svm*F21*F21r*u1**2*c1*f2
     &  + 8.D0*PC_q**2*svp*svm*F21*F20r*u0*u1*c0*f2
     &  + 8.D0*PC_q**2*svp*svm*F20*F25r*u0*u5*c5*f2
     &  + 8.D0*PC_q**2*svp*svm*F20*F24r*u0*u4*c4*f2
     &  + 8.D0*PC_q**2*svp*svm*F20*F23r*u0*u3*c3*f2
     &  + 8.D0*PC_q**2*svp*svm*F20*F22r*u0*u2*c2*f2
     &  + 8.D0*PC_q**2*svp*svm*F20*F21r*u0*u1*c1*f2
     &  + 8.D0*PC_q**2*svp*svm*F20*F20r*u0**2*c0*f2
     &  - 4.D0*PC_q**2*svp*svm*F15*F45r*u5**2*c5*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F15*F45r*u5**2*c5*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F15*F44r*u4*u5*c4*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F15*F44r*u4*u5*c4*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F15*F43r*u3*u5*c3*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F15*F43r*u3*u5*c3*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F15*F42r*u2*u5*c2*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F15*F42r*u2*u5*c2*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F15*F41r*u1*u5*c1*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F15*F41r*u1*u5*c1*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F15*F40r*u0*u5*c0*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F15*F40r*u0*u5*c0*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F14*F45r*u4*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F14*F45r*u4*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F14*F44r*u4**2*c4*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F14*F44r*u4**2*c4*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F14*F43r*u3*u4*c3*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F14*F43r*u3*u4*c3*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F14*F42r*u2*u4*c2*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F14*F42r*u2*u4*c2*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F14*F41r*u1*u4*c1*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F14*F41r*u1*u4*c1*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F14*F40r*u0*u4*c0*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F14*F40r*u0*u4*c0*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F13*F45r*u3*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F13*F45r*u3*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F13*F44r*u3*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F13*F44r*u3*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F13*F43r*u3**2*c3*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F13*F43r*u3**2*c3*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F13*F42r*u2*u3*c2*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F13*F42r*u2*u3*c2*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F13*F41r*u1*u3*c1*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F13*F41r*u1*u3*c1*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F13*F40r*u0*u3*c0*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F13*F40r*u0*u3*c0*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F12*F45r*u2*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F12*F45r*u2*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F12*F44r*u2*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F12*F44r*u2*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F12*F43r*u2*u3*c3*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F12*F43r*u2*u3*c3*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F12*F42r*u2**2*c2*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F12*F42r*u2**2*c2*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F12*F41r*u1*u2*c1*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F12*F41r*u1*u2*c1*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F12*F40r*u0*u2*c0*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F12*F40r*u0*u2*c0*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F11*F45r*u1*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F11*F45r*u1*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F11*F44r*u1*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F11*F44r*u1*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F11*F43r*u1*u3*c3*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F11*F43r*u1*u3*c3*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F11*F42r*u1*u2*c2*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F11*F42r*u1*u2*c2*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F11*F41r*u1**2*c1*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F11*F41r*u1**2*c1*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F11*F40r*u0*u1*c0*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F11*F40r*u0*u1*c0*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F10*F45r*u0*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F10*F45r*u0*u5*c5*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*svm*F10*F44r*u0*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F10*F44r*u0*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F10*F43r*u0*u3*c3*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F10*F43r*u0*u3*c3*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F10*F42r*u0*u2*c2*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F10*F42r*u0*u2*c2*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F10*F41r*u0*u1*c1*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F10*F41r*u0*u1*c1*f4*w1
     &  - 4.D0*PC_q**2*svp*svm*F10*F40r*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**2*svp*svm*F10*F40r*u0**2*c0*f4*w1
     &  + 4.D0*PC_q**2*ssm*svp*F45*F25r*u5**2*c5*f2
     &  + 4.D0*PC_q**2*ssm*svp*F45*F24r*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*ssm*svp*F45*F23r*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*ssm*svp*F45*F22r*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*ssm*svp*F45*F21r*u1*u5*c1*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*svp*F45*F20r*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*ssm*svp*F44*F25r*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*ssm*svp*F44*F24r*u4**2*c4*f2
     &  + 4.D0*PC_q**2*ssm*svp*F44*F23r*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*ssm*svp*F44*F22r*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*ssm*svp*F44*F21r*u1*u4*c1*f2
     &  + 4.D0*PC_q**2*ssm*svp*F44*F20r*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*ssm*svp*F43*F25r*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*ssm*svp*F43*F24r*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*ssm*svp*F43*F23r*u3**2*c3*f2
     &  + 4.D0*PC_q**2*ssm*svp*F43*F22r*u2*u3*c2*f2
     &  + 4.D0*PC_q**2*ssm*svp*F43*F21r*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*ssm*svp*F43*F20r*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*ssm*svp*F42*F25r*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*ssm*svp*F42*F24r*u2*u4*c4*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*svp*F42*F23r*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*ssm*svp*F42*F22r*u2**2*c2*f2
     &  + 4.D0*PC_q**2*ssm*svp*F42*F21r*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*ssm*svp*F42*F20r*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*ssm*svp*F41*F25r*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*ssm*svp*F41*F24r*u1*u4*c4*f2
     &  + 4.D0*PC_q**2*ssm*svp*F41*F23r*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*ssm*svp*F41*F22r*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*ssm*svp*F41*F21r*u1**2*c1*f2
     &  + 4.D0*PC_q**2*ssm*svp*F41*F20r*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*ssm*svp*F40*F25r*u0*u5*c5*f2
     &  + 4.D0*PC_q**2*ssm*svp*F40*F24r*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*ssm*svp*F40*F23r*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*ssm*svp*F40*F22r*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*ssm*svp*F40*F21r*u0*u1*c1*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*svp*F40*F20r*u0**2*c0*f2
     &  - 4.D0*PC_q**2*ssm*svp*F35*F15r*u5**2*c5*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F35*F14r*u4*u5*c4*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F35*F13r*u3*u5*c3*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F35*F12r*u2*u5*c2*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F35*F11r*u1*u5*c1*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F35*F10r*u0*u5*c0*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F34*F15r*u4*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F34*F14r*u4**2*c4*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F34*F13r*u3*u4*c3*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F34*F12r*u2*u4*c2*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F34*F11r*u1*u4*c1*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F34*F10r*u0*u4*c0*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F33*F15r*u3*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F33*F14r*u3*u4*c4*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*svp*F33*F13r*u3**2*c3*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F33*F12r*u2*u3*c2*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F33*F11r*u1*u3*c1*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F33*F10r*u0*u3*c0*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F32*F15r*u2*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F32*F14r*u2*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F32*F13r*u2*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F32*F12r*u2**2*c2*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F32*F11r*u1*u2*c1*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F32*F10r*u0*u2*c0*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F31*F15r*u1*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F31*F14r*u1*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F31*F13r*u1*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F31*F12r*u1*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F31*F11r*u1**2*c1*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*svp*F31*F10r*u0*u1*c0*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F30*F15r*u0*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F30*F14r*u0*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F30*F13r*u0*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F30*F12r*u0*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F30*F11r*u0*u1*c1*f1*w1
     &  - 4.D0*PC_q**2*ssm*svp*F30*F10r*u0**2*c0*f1*w1
     &  + 4.D0*PC_q**2*ssm*svp*F25*F45r*u5**2*c5*f4
     &  + 4.D0*PC_q**2*ssm*svp*F25*F44r*u4*u5*c4*f4
     &  + 4.D0*PC_q**2*ssm*svp*F25*F43r*u3*u5*c3*f4
     &  + 4.D0*PC_q**2*ssm*svp*F25*F42r*u2*u5*c2*f4
     &  + 4.D0*PC_q**2*ssm*svp*F25*F41r*u1*u5*c1*f4
     &  + 4.D0*PC_q**2*ssm*svp*F25*F40r*u0*u5*c0*f4
     &  + 4.D0*PC_q**2*ssm*svp*F24*F45r*u4*u5*c5*f4
     &  + 4.D0*PC_q**2*ssm*svp*F24*F44r*u4**2*c4*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*svp*F24*F43r*u3*u4*c3*f4
     &  + 4.D0*PC_q**2*ssm*svp*F24*F42r*u2*u4*c2*f4
     &  + 4.D0*PC_q**2*ssm*svp*F24*F41r*u1*u4*c1*f4
     &  + 4.D0*PC_q**2*ssm*svp*F24*F40r*u0*u4*c0*f4
     &  + 4.D0*PC_q**2*ssm*svp*F23*F45r*u3*u5*c5*f4
     &  + 4.D0*PC_q**2*ssm*svp*F23*F44r*u3*u4*c4*f4
     &  + 4.D0*PC_q**2*ssm*svp*F23*F43r*u3**2*c3*f4
     &  + 4.D0*PC_q**2*ssm*svp*F23*F42r*u2*u3*c2*f4
     &  + 4.D0*PC_q**2*ssm*svp*F23*F41r*u1*u3*c1*f4
     &  + 4.D0*PC_q**2*ssm*svp*F23*F40r*u0*u3*c0*f4
     &  + 4.D0*PC_q**2*ssm*svp*F22*F45r*u2*u5*c5*f4
     &  + 4.D0*PC_q**2*ssm*svp*F22*F44r*u2*u4*c4*f4
     &  + 4.D0*PC_q**2*ssm*svp*F22*F43r*u2*u3*c3*f4
     &  + 4.D0*PC_q**2*ssm*svp*F22*F42r*u2**2*c2*f4
     &  + 4.D0*PC_q**2*ssm*svp*F22*F41r*u1*u2*c1*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*svp*F22*F40r*u0*u2*c0*f4
     &  + 4.D0*PC_q**2*ssm*svp*F21*F45r*u1*u5*c5*f4
     &  + 4.D0*PC_q**2*ssm*svp*F21*F44r*u1*u4*c4*f4
     &  + 4.D0*PC_q**2*ssm*svp*F21*F43r*u1*u3*c3*f4
     &  + 4.D0*PC_q**2*ssm*svp*F21*F42r*u1*u2*c2*f4
     &  + 4.D0*PC_q**2*ssm*svp*F21*F41r*u1**2*c1*f4
     &  + 4.D0*PC_q**2*ssm*svp*F21*F40r*u0*u1*c0*f4
     &  + 4.D0*PC_q**2*ssm*svp*F20*F45r*u0*u5*c5*f4
     &  + 4.D0*PC_q**2*ssm*svp*F20*F44r*u0*u4*c4*f4
     &  + 4.D0*PC_q**2*ssm*svp*F20*F43r*u0*u3*c3*f4
     &  + 4.D0*PC_q**2*ssm*svp*F20*F42r*u0*u2*c2*f4
     &  + 4.D0*PC_q**2*ssm*svp*F20*F41r*u0*u1*c1*f4
     &  + 4.D0*PC_q**2*ssm*svp*F20*F40r*u0**2*c0*f4
     &  - 4.D0*PC_q**2*ssm*svp*F15*F35r*u5**2*c5*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F15*F34r*u4*u5*c4*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*svp*F15*F33r*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F15*F32r*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F15*F31r*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F15*F30r*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F14*F35r*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F14*F34r*u4**2*c4*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F14*F33r*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F14*F32r*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F14*F31r*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F14*F30r*u0*u4*c0*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F13*F35r*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F13*F34r*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F13*F33r*u3**2*c3*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F13*F32r*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F13*F31r*u1*u3*c1*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*svp*F13*F30r*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F12*F35r*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F12*F34r*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F12*F33r*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F12*F32r*u2**2*c2*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F12*F31r*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F12*F30r*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F11*F35r*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F11*F34r*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F11*F33r*u1*u3*c3*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F11*F32r*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F11*F31r*u1**2*c1*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F11*F30r*u0*u1*c0*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F10*F35r*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F10*F34r*u0*u4*c4*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*svp*F10*F33r*u0*u3*c3*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F10*F32r*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F10*F31r*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q**2*ssm*svp*F10*F30r*u0**2*c0*f3*w1
     &  + 4.D0*PC_q**2*ssp*svm*F45*F25r*u5**2*c5*f2
     &  + 4.D0*PC_q**2*ssp*svm*F45*F24r*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*ssp*svm*F45*F23r*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*ssp*svm*F45*F22r*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*ssp*svm*F45*F21r*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*ssp*svm*F45*F20r*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*ssp*svm*F44*F25r*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*svm*F44*F24r*u4**2*c4*f2
     &  + 4.D0*PC_q**2*ssp*svm*F44*F23r*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*ssp*svm*F44*F22r*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*ssp*svm*F44*F21r*u1*u4*c1*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*svm*F44*F20r*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*ssp*svm*F43*F25r*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*svm*F43*F24r*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*ssp*svm*F43*F23r*u3**2*c3*f2
     &  + 4.D0*PC_q**2*ssp*svm*F43*F22r*u2*u3*c2*f2
     &  + 4.D0*PC_q**2*ssp*svm*F43*F21r*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*ssp*svm*F43*F20r*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*ssp*svm*F42*F25r*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*svm*F42*F24r*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*ssp*svm*F42*F23r*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*ssp*svm*F42*F22r*u2**2*c2*f2
     &  + 4.D0*PC_q**2*ssp*svm*F42*F21r*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*ssp*svm*F42*F20r*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*ssp*svm*F41*F25r*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*svm*F41*F24r*u1*u4*c4*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*svm*F41*F23r*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*ssp*svm*F41*F22r*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*ssp*svm*F41*F21r*u1**2*c1*f2
     &  + 4.D0*PC_q**2*ssp*svm*F41*F20r*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*ssp*svm*F40*F25r*u0*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*svm*F40*F24r*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*ssp*svm*F40*F23r*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*ssp*svm*F40*F22r*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*ssp*svm*F40*F21r*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*ssp*svm*F40*F20r*u0**2*c0*f2
     &  - 4.D0*PC_q**2*ssp*svm*F35*F15r*u5**2*c5*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F35*F14r*u4*u5*c4*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F35*F13r*u3*u5*c3*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F35*F12r*u2*u5*c2*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F35*F11r*u1*u5*c1*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*svm*F35*F10r*u0*u5*c0*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F34*F15r*u4*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F34*F14r*u4**2*c4*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F34*F13r*u3*u4*c3*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F34*F12r*u2*u4*c2*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F34*F11r*u1*u4*c1*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F34*F10r*u0*u4*c0*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F33*F15r*u3*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F33*F14r*u3*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F33*F13r*u3**2*c3*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F33*F12r*u2*u3*c2*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F33*F11r*u1*u3*c1*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F33*F10r*u0*u3*c0*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F32*F15r*u2*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F32*F14r*u2*u4*c4*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*svm*F32*F13r*u2*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F32*F12r*u2**2*c2*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F32*F11r*u1*u2*c1*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F32*F10r*u0*u2*c0*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F31*F15r*u1*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F31*F14r*u1*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F31*F13r*u1*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F31*F12r*u1*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F31*F11r*u1**2*c1*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F31*F10r*u0*u1*c0*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F30*F15r*u0*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F30*F14r*u0*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F30*F13r*u0*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F30*F12r*u0*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*ssp*svm*F30*F11r*u0*u1*c1*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*svm*F30*F10r*u0**2*c0*f1*w2
     &  + 4.D0*PC_q**2*ssp*svm*F25*F45r*u5**2*c5*f4
     &  + 4.D0*PC_q**2*ssp*svm*F25*F44r*u4*u5*c4*f4
     &  + 4.D0*PC_q**2*ssp*svm*F25*F43r*u3*u5*c3*f4
     &  + 4.D0*PC_q**2*ssp*svm*F25*F42r*u2*u5*c2*f4
     &  + 4.D0*PC_q**2*ssp*svm*F25*F41r*u1*u5*c1*f4
     &  + 4.D0*PC_q**2*ssp*svm*F25*F40r*u0*u5*c0*f4
     &  + 4.D0*PC_q**2*ssp*svm*F24*F45r*u4*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*svm*F24*F44r*u4**2*c4*f4
     &  + 4.D0*PC_q**2*ssp*svm*F24*F43r*u3*u4*c3*f4
     &  + 4.D0*PC_q**2*ssp*svm*F24*F42r*u2*u4*c2*f4
     &  + 4.D0*PC_q**2*ssp*svm*F24*F41r*u1*u4*c1*f4
     &  + 4.D0*PC_q**2*ssp*svm*F24*F40r*u0*u4*c0*f4
     &  + 4.D0*PC_q**2*ssp*svm*F23*F45r*u3*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*svm*F23*F44r*u3*u4*c4*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*svm*F23*F43r*u3**2*c3*f4
     &  + 4.D0*PC_q**2*ssp*svm*F23*F42r*u2*u3*c2*f4
     &  + 4.D0*PC_q**2*ssp*svm*F23*F41r*u1*u3*c1*f4
     &  + 4.D0*PC_q**2*ssp*svm*F23*F40r*u0*u3*c0*f4
     &  + 4.D0*PC_q**2*ssp*svm*F22*F45r*u2*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*svm*F22*F44r*u2*u4*c4*f4
     &  + 4.D0*PC_q**2*ssp*svm*F22*F43r*u2*u3*c3*f4
     &  + 4.D0*PC_q**2*ssp*svm*F22*F42r*u2**2*c2*f4
     &  + 4.D0*PC_q**2*ssp*svm*F22*F41r*u1*u2*c1*f4
     &  + 4.D0*PC_q**2*ssp*svm*F22*F40r*u0*u2*c0*f4
     &  + 4.D0*PC_q**2*ssp*svm*F21*F45r*u1*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*svm*F21*F44r*u1*u4*c4*f4
     &  + 4.D0*PC_q**2*ssp*svm*F21*F43r*u1*u3*c3*f4
     &  + 4.D0*PC_q**2*ssp*svm*F21*F42r*u1*u2*c2*f4
     &  + 4.D0*PC_q**2*ssp*svm*F21*F41r*u1**2*c1*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*svm*F21*F40r*u0*u1*c0*f4
     &  + 4.D0*PC_q**2*ssp*svm*F20*F45r*u0*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*svm*F20*F44r*u0*u4*c4*f4
     &  + 4.D0*PC_q**2*ssp*svm*F20*F43r*u0*u3*c3*f4
     &  + 4.D0*PC_q**2*ssp*svm*F20*F42r*u0*u2*c2*f4
     &  + 4.D0*PC_q**2*ssp*svm*F20*F41r*u0*u1*c1*f4
     &  + 4.D0*PC_q**2*ssp*svm*F20*F40r*u0**2*c0*f4
     &  - 4.D0*PC_q**2*ssp*svm*F15*F35r*u5**2*c5*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F15*F34r*u4*u5*c4*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F15*F33r*u3*u5*c3*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F15*F32r*u2*u5*c2*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F15*F31r*u1*u5*c1*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F15*F30r*u0*u5*c0*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F14*F35r*u4*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F14*F34r*u4**2*c4*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*svm*F14*F33r*u3*u4*c3*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F14*F32r*u2*u4*c2*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F14*F31r*u1*u4*c1*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F14*F30r*u0*u4*c0*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F13*F35r*u3*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F13*F34r*u3*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F13*F33r*u3**2*c3*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F13*F32r*u2*u3*c2*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F13*F31r*u1*u3*c1*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F13*F30r*u0*u3*c0*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F12*F35r*u2*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F12*F34r*u2*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F12*F33r*u2*u3*c3*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F12*F32r*u2**2*c2*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F12*F31r*u1*u2*c1*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*svm*F12*F30r*u0*u2*c0*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F11*F35r*u1*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F11*F34r*u1*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F11*F33r*u1*u3*c3*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F11*F32r*u1*u2*c2*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F11*F31r*u1**2*c1*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F11*F30r*u0*u1*c0*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F10*F35r*u0*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F10*F34r*u0*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F10*F33r*u0*u3*c3*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F10*F32r*u0*u2*c2*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F10*F31r*u0*u1*c1*f3*w2
     &  - 4.D0*PC_q**2*ssp*svm*F10*F30r*u0**2*c0*f3*w2
     &  + 4.D0*PC_q**2*ssp*ssm*F45*F45r*u5**2*c5*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F45*F44r*u4*u5*c4*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F45*F43r*u3*u5*c3*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F45*F42r*u2*u5*c2*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F45*F41r*u1*u5*c1*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F45*F40r*u0*u5*c0*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F44*F45r*u4*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F44*F44r*u4**2*c4*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F44*F43r*u3*u4*c3*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F44*F42r*u2*u4*c2*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F44*F41r*u1*u4*c1*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F44*F40r*u0*u4*c0*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F43*F45r*u3*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F43*F44r*u3*u4*c4*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F43*F43r*u3**2*c3*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F43*F42r*u2*u3*c2*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F43*F41r*u1*u3*c1*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F43*F40r*u0*u3*c0*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F42*F45r*u2*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F42*F44r*u2*u4*c4*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F42*F43r*u2*u3*c3*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F42*F42r*u2**2*c2*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F42*F41r*u1*u2*c1*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F42*F40r*u0*u2*c0*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F41*F45r*u1*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F41*F44r*u1*u4*c4*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F41*F43r*u1*u3*c3*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F41*F42r*u1*u2*c2*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F41*F41r*u1**2*c1*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F41*F40r*u0*u1*c0*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F40*F45r*u0*u5*c5*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F40*F44r*u0*u4*c4*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F40*F43r*u0*u3*c3*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F40*F42r*u0*u2*c2*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F40*F41r*u0*u1*c1*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F40*F40r*u0**2*c0*f4
     &  + 4.D0*PC_q**2*ssp*ssm*F35*F25r*u5**2*c5*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F35*F24r*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F35*F23r*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F35*F22r*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F35*F21r*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F35*F20r*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F34*F25r*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F34*F24r*u4**2*c4*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F34*F23r*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F34*F22r*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F34*F21r*u1*u4*c1*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F34*F20r*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F33*F25r*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F33*F24r*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F33*F23r*u3**2*c3*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F33*F22r*u2*u3*c2*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F33*F21r*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F33*F20r*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F32*F25r*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F32*F24r*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F32*F23r*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F32*F22r*u2**2*c2*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F32*F21r*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F32*F20r*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F31*F25r*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F31*F24r*u1*u4*c4*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F31*F23r*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F31*F22r*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F31*F21r*u1**2*c1*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F31*F20r*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F30*F25r*u0*u5*c5*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F30*F24r*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F30*F23r*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F30*F22r*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F30*F21r*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F30*F20r*u0**2*c0*f2
     &  + 4.D0*PC_q**2*ssp*ssm*F25*F35r*u5**2*c5*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F25*F34r*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F25*F33r*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F25*F32r*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F25*F31r*u1*u5*c1*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F25*F30r*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F24*F35r*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F24*F34r*u4**2*c4*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F24*F33r*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F24*F32r*u2*u4*c2*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F24*F31r*u1*u4*c1*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F24*F30r*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F23*F35r*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F23*F34r*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F23*F33r*u3**2*c3*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F23*F32r*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F23*F31r*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F23*F30r*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F22*F35r*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F22*F34r*u2*u4*c4*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F22*F33r*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F22*F32r*u2**2*c2*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F22*F31r*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F22*F30r*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F21*F35r*u1*u5*c5*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F21*F34r*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F21*F33r*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F21*F32r*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F21*F31r*u1**2*c1*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F21*F30r*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F20*F35r*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F20*F34r*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F20*F33r*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F20*F32r*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*ssp*ssm*F20*F31r*u0*u1*c1*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*ssm*F20*F30r*u0**2*c0*f3
     &  - 4.D0*PC_q**2*qq*svp*svm*F45*F45r*u5**2*c5*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F45*F44r*u4*u5*c4*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F45*F43r*u3*u5*c3*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F45*F42r*u2*u5*c2*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F45*F41r*u1*u5*c1*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F45*F40r*u0*u5*c0*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F44*F45r*u4*u5*c5*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F44*F44r*u4**2*c4*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F44*F43r*u3*u4*c3*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F44*F42r*u2*u4*c2*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F44*F41r*u1*u4*c1*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F44*F40r*u0*u4*c0*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F43*F45r*u3*u5*c5*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F43*F44r*u3*u4*c4*f4
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svp*svm*F43*F43r*u3**2*c3*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F43*F42r*u2*u3*c2*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F43*F41r*u1*u3*c1*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F43*F40r*u0*u3*c0*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F42*F45r*u2*u5*c5*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F42*F44r*u2*u4*c4*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F42*F43r*u2*u3*c3*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F42*F42r*u2**2*c2*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F42*F41r*u1*u2*c1*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F42*F40r*u0*u2*c0*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F41*F45r*u1*u5*c5*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F41*F44r*u1*u4*c4*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F41*F43r*u1*u3*c3*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F41*F42r*u1*u2*c2*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F41*F41r*u1**2*c1*f4
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svp*svm*F41*F40r*u0*u1*c0*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F40*F45r*u0*u5*c5*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F40*F44r*u0*u4*c4*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F40*F43r*u0*u3*c3*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F40*F42r*u0*u2*c2*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F40*F41r*u0*u1*c1*f4
     &  - 4.D0*PC_q**2*qq*svp*svm*F40*F40r*u0**2*c0*f4
     &  + 4.D0*PC_q**2*qq*svp*svm*F35*F25r*u5**2*c5*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F35*F24r*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F35*F23r*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F35*F22r*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F35*F21r*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F35*F20r*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F34*F25r*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F34*F24r*u4**2*c4*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*svm*F34*F23r*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F34*F22r*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F34*F21r*u1*u4*c1*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F34*F20r*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F33*F25r*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F33*F24r*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F33*F23r*u3**2*c3*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F33*F22r*u2*u3*c2*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F33*F21r*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F33*F20r*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F32*F25r*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F32*F24r*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F32*F23r*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F32*F22r*u2**2*c2*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F32*F21r*u1*u2*c1*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*svm*F32*F20r*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F31*F25r*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F31*F24r*u1*u4*c4*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F31*F23r*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F31*F22r*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F31*F21r*u1**2*c1*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F31*F20r*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F30*F25r*u0*u5*c5*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F30*F24r*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F30*F23r*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F30*F22r*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F30*F21r*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F30*F20r*u0**2*c0*f2
     &  + 4.D0*PC_q**2*qq*svp*svm*F25*F35r*u5**2*c5*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F25*F34r*u4*u5*c4*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*svm*F25*F33r*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F25*F32r*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F25*F31r*u1*u5*c1*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F25*F30r*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F24*F35r*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F24*F34r*u4**2*c4*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F24*F33r*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F24*F32r*u2*u4*c2*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F24*F31r*u1*u4*c1*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F24*F30r*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F23*F35r*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F23*F34r*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F23*F33r*u3**2*c3*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F23*F32r*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F23*F31r*u1*u3*c1*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*svm*F23*F30r*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F22*F35r*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F22*F34r*u2*u4*c4*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F22*F33r*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F22*F32r*u2**2*c2*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F22*F31r*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F22*F30r*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F21*F35r*u1*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F21*F34r*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F21*F33r*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F21*F32r*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F21*F31r*u1**2*c1*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F21*F30r*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F20*F35r*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F20*F34r*u0*u4*c4*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*svm*F20*F33r*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F20*F32r*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F20*F31r*u0*u1*c1*f3
     &  + 4.D0*PC_q**2*qq*svp*svm*F20*F30r*u0**2*c0*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F35*F35r*u5**2*c5*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F35*F34r*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F35*F33r*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F35*F32r*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F35*F31r*u1*u5*c1*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F35*F30r*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F34*F35r*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F34*F34r*u4**2*c4*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F34*F33r*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F34*F32r*u2*u4*c2*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F34*F31r*u1*u4*c1*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssp*ssm*F34*F30r*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F33*F35r*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F33*F34r*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F33*F33r*u3**2*c3*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F33*F32r*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F33*F31r*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F33*F30r*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F32*F35r*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F32*F34r*u2*u4*c4*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F32*F33r*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F32*F32r*u2**2*c2*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F32*F31r*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F32*F30r*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F31*F35r*u1*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F31*F34r*u1*u4*c4*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssp*ssm*F31*F33r*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F31*F32r*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F31*F31r*u1**2*c1*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F31*F30r*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F30*F35r*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F30*F34r*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F30*F33r*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F30*F32r*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F30*F31r*u0*u1*c1*f3
     &  + 4.D0*PC_q**2*qq*ssp*ssm*F30*F30r*u0**2*c0*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F35*F35r*u5**2*c5*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F35*F34r*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F35*F33r*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F35*F32r*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F35*F31r*u1*u5*c1*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svp*svm*F35*F30r*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F34*F35r*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F34*F34r*u4**2*c4*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F34*F33r*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F34*F32r*u2*u4*c2*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F34*F31r*u1*u4*c1*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F34*F30r*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F33*F35r*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F33*F34r*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F33*F33r*u3**2*c3*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F33*F32r*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F33*F31r*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F33*F30r*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F32*F35r*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F32*F34r*u2*u4*c4*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svp*svm*F32*F33r*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F32*F32r*u2**2*c2*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F32*F31r*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F32*F30r*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F31*F35r*u1*u5*c5*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F31*F34r*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F31*F33r*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F31*F32r*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F31*F31r*u1**2*c1*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F31*F30r*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F30*F35r*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F30*F34r*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F30*F33r*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F30*F32r*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*qq**2*svp*svm*F30*F31r*u0*u1*c1*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svp*svm*F30*F30r*u0**2*c0*f3
     &  + 4.D0*PC_q**2*PC2*svp*svm*F45*F45r*u5**2*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F45*F44r*u4*u5*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F45*F43r*u3*u5*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F45*F42r*u2*u5*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F45*F41r*u1*u5*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F45*F40r*u0*u5*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F44*F45r*u4*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F44*F44r*u4**2*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F44*F43r*u3*u4*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F44*F42r*u2*u4*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F44*F41r*u1*u4*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F44*F40r*u0*u4*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F43*F45r*u3*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F43*F44r*u3*u4*c4*f4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svp*svm*F43*F43r*u3**2*c3*f4*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F43*F42r*u2*u3*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F43*F41r*u1*u3*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F43*F40r*u0*u3*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F42*F45r*u2*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F42*F44r*u2*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F42*F43r*u2*u3*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F42*F42r*u2**2*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F42*F41r*u1*u2*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F42*F40r*u0*u2*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F41*F45r*u1*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F41*F44r*u1*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F41*F43r*u1*u3*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F41*F42r*u1*u2*c2*f4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svp*svm*F41*F41r*u1**2*c1*f4*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F41*F40r*u0*u1*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F40*F45r*u0*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F40*F44r*u0*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F40*F43r*u0*u3*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F40*F42r*u0*u2*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F40*F41r*u0*u1*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*svm*F40*F40r*u0**2*c0*f4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F35*F25r*u5**2*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F35*F24r*u4*u5*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F35*F23r*u3*u5*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F35*F22r*u2*u5*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F35*F21r*u1*u5*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F35*F20r*u0*u5*c0*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*svm*F34*F25r*u4*u5*c5*f2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F34*F24r*u4**2*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F34*F23r*u3*u4*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F34*F22r*u2*u4*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F34*F21r*u1*u4*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F34*F20r*u0*u4*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F33*F25r*u3*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F33*F24r*u3*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F33*F23r*u3**2*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F33*F22r*u2*u3*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F33*F21r*u1*u3*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F33*F20r*u0*u3*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F32*F25r*u2*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F32*F24r*u2*u4*c4*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*svm*F32*F23r*u2*u3*c3*f2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F32*F22r*u2**2*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F32*F21r*u1*u2*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F32*F20r*u0*u2*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F31*F25r*u1*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F31*F24r*u1*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F31*F23r*u1*u3*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F31*F22r*u1*u2*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F31*F21r*u1**2*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F31*F20r*u0*u1*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F30*F25r*u0*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F30*F24r*u0*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F30*F23r*u0*u3*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F30*F22r*u0*u2*c2*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*svm*F30*F21r*u0*u1*c1*f2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F30*F20r*u0**2*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F25*F35r*u5**2*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F25*F34r*u4*u5*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F25*F33r*u3*u5*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F25*F32r*u2*u5*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F25*F31r*u1*u5*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F25*F30r*u0*u5*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F24*F35r*u4*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F24*F34r*u4**2*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F24*F33r*u3*u4*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F24*F32r*u2*u4*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F24*F31r*u1*u4*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F24*F30r*u0*u4*c0*f3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*svm*F23*F35r*u3*u5*c5*f3*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F23*F34r*u3*u4*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F23*F33r*u3**2*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F23*F32r*u2*u3*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F23*F31r*u1*u3*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F23*F30r*u0*u3*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F22*F35r*u2*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F22*F34r*u2*u4*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F22*F33r*u2*u3*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F22*F32r*u2**2*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F22*F31r*u1*u2*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F22*F30r*u0*u2*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F21*F35r*u1*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F21*F34r*u1*u4*c4*f3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*svm*F21*F33r*u1*u3*c3*f3*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F21*F32r*u1*u2*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F21*F31r*u1**2*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F21*F30r*u0*u1*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F20*F35r*u0*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F20*F34r*u0*u4*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F20*F33r*u0*u3*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F20*F32r*u0*u2*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F20*F31r*u0*u1*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*svm*F20*F30r*u0**2*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F35*F35r*u5**2*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F35*F34r*u4*u5*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F35*F33r*u3*u5*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F35*F32r*u2*u5*c2*f3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svp*svm*F35*F31r*u1*u5*c1*
     & f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F35*F30r*u0*u5*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F34*F35r*u4*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F34*F34r*u4**2*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F34*F33r*u3*u4*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F34*F32r*u2*u4*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F34*F31r*u1*u4*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F34*F30r*u0*u4*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F33*F35r*u3*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F33*F34r*u3*u4*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F33*F33r*u3**2*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F33*F32r*u2*u3*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F33*F31r*u1*u3*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F33*F30r*u0*u3*c0*f3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svp*svm*F32*F35r*u2*u5*c5*
     & f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F32*F34r*u2*u4*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F32*F33r*u2*u3*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F32*F32r*u2**2*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F32*F31r*u1*u2*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F32*F30r*u0*u2*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F31*F35r*u1*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F31*F34r*u1*u4*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F31*F33r*u1*u3*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F31*F32r*u1*u2*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F31*F31r*u1**2*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F31*F30r*u0*u1*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F30*F35r*u0*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F30*F34r*u0*u4*c4*f3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svp*svm*F30*F33r*u0*u3*c3*
     & f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F30*F32r*u0*u2*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F30*F31r*u0*u1*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*svm*F30*F30r*u0**2*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F15i*u5**2*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F15i*u5**2*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F14i*u4*u5*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F14i*u4*u5*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F13i*u3*u5*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F13i*u3*u5*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F12i*u2*u5*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F12i*u2*u5*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F11i*u1*u5*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F11i*u1*u5*c1*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F45*F10i*u0*u5*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F45*F10i*u0*u5*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F15i*u4*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F15i*u4*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F14i*u4**2*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F14i*u4**2*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F13i*u3*u4*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F13i*u3*u4*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F12i*u2*u4*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F12i*u2*u4*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F11i*u1*u4*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F11i*u1*u4*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F10i*u0*u4*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F44*F10i*u0*u4*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F15i*u3*u5*c5*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F43*F15i*u3*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F14i*u3*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F14i*u3*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F13i*u3**2*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F13i*u3**2*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F12i*u2*u3*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F12i*u2*u3*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F11i*u1*u3*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F11i*u1*u3*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F10i*u0*u3*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F43*F10i*u0*u3*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F15i*u2*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F15i*u2*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F14i*u2*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F14i*u2*u4*c4*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F42*F13i*u2*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F13i*u2*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F12i*u2**2*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F12i*u2**2*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F11i*u1*u2*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F11i*u1*u2*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F10i*u0*u2*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F42*F10i*u0*u2*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F15i*u1*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F15i*u1*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F14i*u1*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F14i*u1*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F13i*u1*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F13i*u1*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F12i*u1*u2*c2*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F41*F12i*u1*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F11i*u1**2*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F11i*u1**2*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F10i*u0*u1*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F41*F10i*u0*u1*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F15i*u0*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F15i*u0*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F14i*u0*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F14i*u0*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F13i*u0*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F13i*u0*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F12i*u0*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F12i*u0*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F11i*u0*u1*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F11i*u0*u1*c1*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F40*F10i*u0**2*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F40*F10i*u0**2*c0*f1*w1
     &  + 8.D0*PC_q**2*i_*svp*svm*F25*F25i*u5**2*c5*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F25*F24i*u4*u5*c4*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F25*F23i*u3*u5*c3*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F25*F22i*u2*u5*c2*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F25*F21i*u1*u5*c1*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F25*F20i*u0*u5*c0*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F24*F25i*u4*u5*c5*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F24*F24i*u4**2*c4*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F24*F23i*u3*u4*c3*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F24*F22i*u2*u4*c2*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F24*F21i*u1*u4*c1*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F24*F20i*u0*u4*c0*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F23*F25i*u3*u5*c5*f2
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svp*svm*F23*F24i*u3*u4*c4*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F23*F23i*u3**2*c3*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F23*F22i*u2*u3*c2*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F23*F21i*u1*u3*c1*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F23*F20i*u0*u3*c0*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F22*F25i*u2*u5*c5*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F22*F24i*u2*u4*c4*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F22*F23i*u2*u3*c3*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F22*F22i*u2**2*c2*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F22*F21i*u1*u2*c1*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F22*F20i*u0*u2*c0*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F21*F25i*u1*u5*c5*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F21*F24i*u1*u4*c4*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F21*F23i*u1*u3*c3*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F21*F22i*u1*u2*c2*f2
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svp*svm*F21*F21i*u1**2*c1*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F21*F20i*u0*u1*c0*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F20*F25i*u0*u5*c5*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F20*F24i*u0*u4*c4*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F20*F23i*u0*u3*c3*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F20*F22i*u0*u2*c2*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F20*F21i*u0*u1*c1*f2
     &  + 8.D0*PC_q**2*i_*svp*svm*F20*F20i*u0**2*c0*f2
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F45i*u5**2*c5*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F45i*u5**2*c5*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F44i*u4*u5*c4*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F44i*u4*u5*c4*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F43i*u3*u5*c3*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F43i*u3*u5*c3*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F42i*u2*u5*c2*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F15*F42i*u2*u5*c2*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F41i*u1*u5*c1*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F41i*u1*u5*c1*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F40i*u0*u5*c0*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F15*F40i*u0*u5*c0*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F45i*u4*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F45i*u4*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F44i*u4**2*c4*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F44i*u4**2*c4*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F43i*u3*u4*c3*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F43i*u3*u4*c3*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F42i*u2*u4*c2*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F42i*u2*u4*c2*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F41i*u1*u4*c1*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F41i*u1*u4*c1*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F14*F40i*u0*u4*c0*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F14*F40i*u0*u4*c0*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F45i*u3*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F45i*u3*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F44i*u3*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F44i*u3*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F43i*u3**2*c3*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F43i*u3**2*c3*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F42i*u2*u3*c2*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F42i*u2*u3*c2*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F41i*u1*u3*c1*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F41i*u1*u3*c1*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F40i*u0*u3*c0*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F13*F40i*u0*u3*c0*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F45i*u2*u5*c5*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F12*F45i*u2*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F44i*u2*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F44i*u2*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F43i*u2*u3*c3*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F43i*u2*u3*c3*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F42i*u2**2*c2*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F42i*u2**2*c2*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F41i*u1*u2*c1*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F41i*u1*u2*c1*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F40i*u0*u2*c0*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F12*F40i*u0*u2*c0*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F45i*u1*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F45i*u1*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F44i*u1*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F44i*u1*u4*c4*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F11*F43i*u1*u3*c3*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F43i*u1*u3*c3*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F42i*u1*u2*c2*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F42i*u1*u2*c2*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F41i*u1**2*c1*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F41i*u1**2*c1*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F40i*u0*u1*c0*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F11*F40i*u0*u1*c0*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F45i*u0*u5*c5*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F45i*u0*u5*c5*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F44i*u0*u4*c4*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F44i*u0*u4*c4*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F43i*u0*u3*c3*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F43i*u0*u3*c3*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F42i*u0*u2*c2*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*svm*F10*F42i*u0*u2*c2*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F41i*u0*u1*c1*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F41i*u0*u1*c1*f4*w1
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F40i*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**2*i_*svp*svm*F10*F40i*u0**2*c0*f4*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F45*F25i*u5**2*c5*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F45*F24i*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F45*F23i*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F45*F22i*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F45*F21i*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F45*F20i*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F44*F25i*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F44*F24i*u4**2*c4*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F44*F23i*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F44*F22i*u2*u4*c2*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*svp*F44*F21i*u1*u4*c1*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F44*F20i*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F43*F25i*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F43*F24i*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F43*F23i*u3**2*c3*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F43*F22i*u2*u3*c2*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F43*F21i*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F43*F20i*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F42*F25i*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F42*F24i*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F42*F23i*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F42*F22i*u2**2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F42*F21i*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F42*F20i*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F41*F25i*u1*u5*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*svp*F41*F24i*u1*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F41*F23i*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F41*F22i*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F41*F21i*u1**2*c1*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F41*F20i*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F40*F25i*u0*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F40*F24i*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F40*F23i*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F40*F22i*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F40*F21i*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*i_*ssm*svp*F40*F20i*u0**2*c0*f2
     &  - 4.D0*PC_q**2*i_*ssm*svp*F35*F15i*u5**2*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F35*F14i*u4*u5*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F35*F13i*u3*u5*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F35*F12i*u2*u5*c2*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*svp*F35*F11i*u1*u5*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F35*F10i*u0*u5*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F34*F15i*u4*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F34*F14i*u4**2*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F34*F13i*u3*u4*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F34*F12i*u2*u4*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F34*F11i*u1*u4*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F34*F10i*u0*u4*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F33*F15i*u3*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F33*F14i*u3*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F33*F13i*u3**2*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F33*F12i*u2*u3*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F33*F11i*u1*u3*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F33*F10i*u0*u3*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F32*F15i*u2*u5*c5*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*svp*F32*F14i*u2*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F32*F13i*u2*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F32*F12i*u2**2*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F32*F11i*u1*u2*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F32*F10i*u0*u2*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F31*F15i*u1*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F31*F14i*u1*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F31*F13i*u1*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F31*F12i*u1*u2*c2*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F31*F11i*u1**2*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F31*F10i*u0*u1*c0*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F30*F15i*u0*u5*c5*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F30*F14i*u0*u4*c4*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F30*F13i*u0*u3*c3*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F30*F12i*u0*u2*c2*f1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*svp*F30*F11i*u0*u1*c1*f1*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F30*F10i*u0**2*c0*f1*w1
     &  + 4.D0*PC_q**2*i_*ssm*svp*F25*F45i*u5**2*c5*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F25*F44i*u4*u5*c4*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F25*F43i*u3*u5*c3*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F25*F42i*u2*u5*c2*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F25*F41i*u1*u5*c1*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F25*F40i*u0*u5*c0*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F24*F45i*u4*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F24*F44i*u4**2*c4*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F24*F43i*u3*u4*c3*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F24*F42i*u2*u4*c2*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F24*F41i*u1*u4*c1*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F24*F40i*u0*u4*c0*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F23*F45i*u3*u5*c5*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*svp*F23*F44i*u3*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F23*F43i*u3**2*c3*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F23*F42i*u2*u3*c2*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F23*F41i*u1*u3*c1*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F23*F40i*u0*u3*c0*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F22*F45i*u2*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F22*F44i*u2*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F22*F43i*u2*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F22*F42i*u2**2*c2*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F22*F41i*u1*u2*c1*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F22*F40i*u0*u2*c0*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F21*F45i*u1*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F21*F44i*u1*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F21*F43i*u1*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F21*F42i*u1*u2*c2*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*svp*F21*F41i*u1**2*c1*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F21*F40i*u0*u1*c0*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F20*F45i*u0*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F20*F44i*u0*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F20*F43i*u0*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F20*F42i*u0*u2*c2*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F20*F41i*u0*u1*c1*f4
     &  + 4.D0*PC_q**2*i_*ssm*svp*F20*F40i*u0**2*c0*f4
     &  - 4.D0*PC_q**2*i_*ssm*svp*F15*F35i*u5**2*c5*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F15*F34i*u4*u5*c4*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F15*F33i*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F15*F32i*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F15*F31i*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F15*F30i*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F14*F35i*u4*u5*c5*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*svp*F14*F34i*u4**2*c4*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F14*F33i*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F14*F32i*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F14*F31i*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F14*F30i*u0*u4*c0*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F13*F35i*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F13*F34i*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F13*F33i*u3**2*c3*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F13*F32i*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F13*F31i*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F13*F30i*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F12*F35i*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F12*F34i*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F12*F33i*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F12*F32i*u2**2*c2*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*svp*F12*F31i*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F12*F30i*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F11*F35i*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F11*F34i*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F11*F33i*u1*u3*c3*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F11*F32i*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F11*F31i*u1**2*c1*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F11*F30i*u0*u1*c0*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F10*F35i*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F10*F34i*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F10*F33i*u0*u3*c3*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F10*F32i*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F10*F31i*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q**2*i_*ssm*svp*F10*F30i*u0**2*c0*f3*w1
     &  + 4.D0*PC_q**2*i_*ssp*svm*F45*F25i*u5**2*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*svm*F45*F24i*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F45*F23i*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F45*F22i*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F45*F21i*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F45*F20i*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F44*F25i*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F44*F24i*u4**2*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F44*F23i*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F44*F22i*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F44*F21i*u1*u4*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F44*F20i*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F43*F25i*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F43*F24i*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F43*F23i*u3**2*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F43*F22i*u2*u3*c2*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*svm*F43*F21i*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F43*F20i*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F42*F25i*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F42*F24i*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F42*F23i*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F42*F22i*u2**2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F42*F21i*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F42*F20i*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F41*F25i*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F41*F24i*u1*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F41*F23i*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F41*F22i*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F41*F21i*u1**2*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F41*F20i*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F40*F25i*u0*u5*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*svm*F40*F24i*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F40*F23i*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F40*F22i*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F40*F21i*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F40*F20i*u0**2*c0*f2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F35*F15i*u5**2*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F35*F14i*u4*u5*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F35*F13i*u3*u5*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F35*F12i*u2*u5*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F35*F11i*u1*u5*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F35*F10i*u0*u5*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F34*F15i*u4*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F34*F14i*u4**2*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F34*F13i*u3*u4*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F34*F12i*u2*u4*c2*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*svm*F34*F11i*u1*u4*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F34*F10i*u0*u4*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F33*F15i*u3*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F33*F14i*u3*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F33*F13i*u3**2*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F33*F12i*u2*u3*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F33*F11i*u1*u3*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F33*F10i*u0*u3*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F32*F15i*u2*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F32*F14i*u2*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F32*F13i*u2*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F32*F12i*u2**2*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F32*F11i*u1*u2*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F32*F10i*u0*u2*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F31*F15i*u1*u5*c5*f1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*svm*F31*F14i*u1*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F31*F13i*u1*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F31*F12i*u1*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F31*F11i*u1**2*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F31*F10i*u0*u1*c0*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F30*F15i*u0*u5*c5*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F30*F14i*u0*u4*c4*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F30*F13i*u0*u3*c3*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F30*F12i*u0*u2*c2*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F30*F11i*u0*u1*c1*f1*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F30*F10i*u0**2*c0*f1*w2
     &  + 4.D0*PC_q**2*i_*ssp*svm*F25*F45i*u5**2*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F25*F44i*u4*u5*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F25*F43i*u3*u5*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F25*F42i*u2*u5*c2*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*svm*F25*F41i*u1*u5*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F25*F40i*u0*u5*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F24*F45i*u4*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F24*F44i*u4**2*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F24*F43i*u3*u4*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F24*F42i*u2*u4*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F24*F41i*u1*u4*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F24*F40i*u0*u4*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F23*F45i*u3*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F23*F44i*u3*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F23*F43i*u3**2*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F23*F42i*u2*u3*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F23*F41i*u1*u3*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F23*F40i*u0*u3*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F22*F45i*u2*u5*c5*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*svm*F22*F44i*u2*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F22*F43i*u2*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F22*F42i*u2**2*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F22*F41i*u1*u2*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F22*F40i*u0*u2*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F21*F45i*u1*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F21*F44i*u1*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F21*F43i*u1*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F21*F42i*u1*u2*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F21*F41i*u1**2*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F21*F40i*u0*u1*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F20*F45i*u0*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F20*F44i*u0*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F20*F43i*u0*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F20*F42i*u0*u2*c2*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*svm*F20*F41i*u0*u1*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*svm*F20*F40i*u0**2*c0*f4
     &  - 4.D0*PC_q**2*i_*ssp*svm*F15*F35i*u5**2*c5*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F15*F34i*u4*u5*c4*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F15*F33i*u3*u5*c3*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F15*F32i*u2*u5*c2*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F15*F31i*u1*u5*c1*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F15*F30i*u0*u5*c0*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F14*F35i*u4*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F14*F34i*u4**2*c4*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F14*F33i*u3*u4*c3*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F14*F32i*u2*u4*c2*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F14*F31i*u1*u4*c1*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F14*F30i*u0*u4*c0*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F13*F35i*u3*u5*c5*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*svm*F13*F34i*u3*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F13*F33i*u3**2*c3*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F13*F32i*u2*u3*c2*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F13*F31i*u1*u3*c1*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F13*F30i*u0*u3*c0*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F12*F35i*u2*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F12*F34i*u2*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F12*F33i*u2*u3*c3*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F12*F32i*u2**2*c2*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F12*F31i*u1*u2*c1*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F12*F30i*u0*u2*c0*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F11*F35i*u1*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F11*F34i*u1*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F11*F33i*u1*u3*c3*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F11*F32i*u1*u2*c2*f3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*svm*F11*F31i*u1**2*c1*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F11*F30i*u0*u1*c0*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F10*F35i*u0*u5*c5*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F10*F34i*u0*u4*c4*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F10*F33i*u0*u3*c3*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F10*F32i*u0*u2*c2*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F10*F31i*u0*u1*c1*f3*w2
     &  - 4.D0*PC_q**2*i_*ssp*svm*F10*F30i*u0**2*c0*f3*w2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F45*F45i*u5**2*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F45*F44i*u4*u5*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F45*F43i*u3*u5*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F45*F42i*u2*u5*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F45*F41i*u1*u5*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F45*F40i*u0*u5*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F44*F45i*u4*u5*c5*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F44*F44i*u4**2*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F44*F43i*u3*u4*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F44*F42i*u2*u4*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F44*F41i*u1*u4*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F44*F40i*u0*u4*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F43*F45i*u3*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F43*F44i*u3*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F43*F43i*u3**2*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F43*F42i*u2*u3*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F43*F41i*u1*u3*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F43*F40i*u0*u3*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F42*F45i*u2*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F42*F44i*u2*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F42*F43i*u2*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F42*F42i*u2**2*c2*f4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F42*F41i*u1*u2*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F42*F40i*u0*u2*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F41*F45i*u1*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F41*F44i*u1*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F41*F43i*u1*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F41*F42i*u1*u2*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F41*F41i*u1**2*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F41*F40i*u0*u1*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F40*F45i*u0*u5*c5*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F40*F44i*u0*u4*c4*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F40*F43i*u0*u3*c3*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F40*F42i*u0*u2*c2*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F40*F41i*u0*u1*c1*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F40*F40i*u0**2*c0*f4
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F35*F25i*u5**2*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F35*F24i*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F35*F23i*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F35*F22i*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F35*F21i*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F35*F20i*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F34*F25i*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F34*F24i*u4**2*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F34*F23i*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F34*F22i*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F34*F21i*u1*u4*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F34*F20i*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F33*F25i*u3*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F33*F24i*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F33*F23i*u3**2*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F33*F22i*u2*u3*c2*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F33*F21i*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F33*F20i*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F32*F25i*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F32*F24i*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F32*F23i*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F32*F22i*u2**2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F32*F21i*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F32*F20i*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F31*F25i*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F31*F24i*u1*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F31*F23i*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F31*F22i*u1*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F31*F21i*u1**2*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F31*F20i*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F30*F25i*u0*u5*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F30*F24i*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F30*F23i*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F30*F22i*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F30*F21i*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F30*F20i*u0**2*c0*f2
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F25*F35i*u5**2*c5*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F25*F34i*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F25*F33i*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F25*F32i*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F25*F31i*u1*u5*c1*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F25*F30i*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F24*F35i*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F24*F34i*u4**2*c4*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F24*F33i*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F24*F32i*u2*u4*c2*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F24*F31i*u1*u4*c1*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F24*F30i*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F23*F35i*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F23*F34i*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F23*F33i*u3**2*c3*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F23*F32i*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F23*F31i*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F23*F30i*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F22*F35i*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F22*F34i*u2*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F22*F33i*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F22*F32i*u2**2*c2*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F22*F31i*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F22*F30i*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F21*F35i*u1*u5*c5*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*ssm*F21*F34i*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F21*F33i*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F21*F32i*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F21*F31i*u1**2*c1*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F21*F30i*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F20*F35i*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F20*F34i*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F20*F33i*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F20*F32i*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F20*F31i*u0*u1*c1*f3
     &  + 4.D0*PC_q**2*i_*ssp*ssm*F20*F30i*u0**2*c0*f3
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F45*F45i*u5**2*c5*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F45*F44i*u4*u5*c4*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F45*F43i*u3*u5*c3*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F45*F42i*u2*u5*c2*f4
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svp*svm*F45*F41i*u1*u5*c1*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F45*F40i*u0*u5*c0*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F44*F45i*u4*u5*c5*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F44*F44i*u4**2*c4*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F44*F43i*u3*u4*c3*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F44*F42i*u2*u4*c2*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F44*F41i*u1*u4*c1*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F44*F40i*u0*u4*c0*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F43*F45i*u3*u5*c5*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F43*F44i*u3*u4*c4*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F43*F43i*u3**2*c3*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F43*F42i*u2*u3*c2*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F43*F41i*u1*u3*c1*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F43*F40i*u0*u3*c0*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F42*F45i*u2*u5*c5*f4
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svp*svm*F42*F44i*u2*u4*c4*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F42*F43i*u2*u3*c3*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F42*F42i*u2**2*c2*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F42*F41i*u1*u2*c1*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F42*F40i*u0*u2*c0*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F41*F45i*u1*u5*c5*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F41*F44i*u1*u4*c4*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F41*F43i*u1*u3*c3*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F41*F42i*u1*u2*c2*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F41*F41i*u1**2*c1*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F41*F40i*u0*u1*c0*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F40*F45i*u0*u5*c5*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F40*F44i*u0*u4*c4*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F40*F43i*u0*u3*c3*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F40*F42i*u0*u2*c2*f4
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svp*svm*F40*F41i*u0*u1*c1*f4
     &  - 4.D0*PC_q**2*i_*qq*svp*svm*F40*F40i*u0**2*c0*f4
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F35*F25i*u5**2*c5*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F35*F24i*u4*u5*c4*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F35*F23i*u3*u5*c3*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F35*F22i*u2*u5*c2*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F35*F21i*u1*u5*c1*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F35*F20i*u0*u5*c0*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F34*F25i*u4*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F34*F24i*u4**2*c4*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F34*F23i*u3*u4*c3*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F34*F22i*u2*u4*c2*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F34*F21i*u1*u4*c1*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F34*F20i*u0*u4*c0*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F33*F25i*u3*u5*c5*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*svm*F33*F24i*u3*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F33*F23i*u3**2*c3*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F33*F22i*u2*u3*c2*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F33*F21i*u1*u3*c1*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F33*F20i*u0*u3*c0*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F32*F25i*u2*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F32*F24i*u2*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F32*F23i*u2*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F32*F22i*u2**2*c2*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F32*F21i*u1*u2*c1*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F32*F20i*u0*u2*c0*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F31*F25i*u1*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F31*F24i*u1*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F31*F23i*u1*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F31*F22i*u1*u2*c2*f2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*svm*F31*F21i*u1**2*c1*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F31*F20i*u0*u1*c0*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F30*F25i*u0*u5*c5*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F30*F24i*u0*u4*c4*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F30*F23i*u0*u3*c3*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F30*F22i*u0*u2*c2*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F30*F21i*u0*u1*c1*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F30*F20i*u0**2*c0*f2
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F25*F35i*u5**2*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F25*F34i*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F25*F33i*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F25*F32i*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F25*F31i*u1*u5*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F25*F30i*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F24*F35i*u4*u5*c5*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*svm*F24*F34i*u4**2*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F24*F33i*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F24*F32i*u2*u4*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F24*F31i*u1*u4*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F24*F30i*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F23*F35i*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F23*F34i*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F23*F33i*u3**2*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F23*F32i*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F23*F31i*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F23*F30i*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F22*F35i*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F22*F34i*u2*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F22*F33i*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F22*F32i*u2**2*c2*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*svm*F22*F31i*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F22*F30i*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F21*F35i*u1*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F21*F34i*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F21*F33i*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F21*F32i*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F21*F31i*u1**2*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F21*F30i*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F20*F35i*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F20*F34i*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F20*F33i*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F20*F32i*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F20*F31i*u0*u1*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*svp*svm*F20*F30i*u0**2*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F35*F35i*u5**2*c5*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssp*ssm*F35*F34i*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F35*F33i*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F35*F32i*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F35*F31i*u1*u5*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F35*F30i*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F34*F35i*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F34*F34i*u4**2*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F34*F33i*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F34*F32i*u2*u4*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F34*F31i*u1*u4*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F34*F30i*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F33*F35i*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F33*F34i*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F33*F33i*u3**2*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F33*F32i*u2*u3*c2*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssp*ssm*F33*F31i*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F33*F30i*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F32*F35i*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F32*F34i*u2*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F32*F33i*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F32*F32i*u2**2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F32*F31i*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F32*F30i*u0*u2*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F31*F35i*u1*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F31*F34i*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F31*F33i*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F31*F32i*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F31*F31i*u1**2*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F31*F30i*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F30*F35i*u0*u5*c5*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssp*ssm*F30*F34i*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F30*F33i*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F30*F32i*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F30*F31i*u0*u1*c1*f3
     &  + 4.D0*PC_q**2*i_*qq*ssp*ssm*F30*F30i*u0**2*c0*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F35*F35i*u5**2*c5*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F35*F34i*u4*u5*c4*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F35*F33i*u3*u5*c3*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F35*F32i*u2*u5*c2*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F35*F31i*u1*u5*c1*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F35*F30i*u0*u5*c0*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F34*F35i*u4*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F34*F34i*u4**2*c4*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F34*F33i*u3*u4*c3*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F34*F32i*u2*u4*c2*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svp*svm*F34*F31i*u1*u4*c1
     & *f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F34*F30i*u0*u4*c0*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F33*F35i*u3*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F33*F34i*u3*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F33*F33i*u3**2*c3*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F33*F32i*u2*u3*c2*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F33*F31i*u1*u3*c1*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F33*F30i*u0*u3*c0*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F32*F35i*u2*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F32*F34i*u2*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F32*F33i*u2*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F32*F32i*u2**2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F32*F31i*u1*u2*c1*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F32*F30i*u0*u2*c0*f3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svp*svm*F31*F35i*u1*u5*c5
     & *f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F31*F34i*u1*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F31*F33i*u1*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F31*F32i*u1*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F31*F31i*u1**2*c1*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F31*F30i*u0*u1*c0*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F30*F35i*u0*u5*c5*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F30*F34i*u0*u4*c4*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F30*F33i*u0*u3*c3*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F30*F32i*u0*u2*c2*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F30*F31i*u0*u1*c1*f3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*svm*F30*F30i*u0**2*c0*f3
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F45*F45i*u5**2*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F45*F44i*u4*u5*c4*f4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svp*svm*F45*F43i*u3*u5*c3*
     & f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F45*F42i*u2*u5*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F45*F41i*u1*u5*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F45*F40i*u0*u5*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F44*F45i*u4*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F44*F44i*u4**2*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F44*F43i*u3*u4*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F44*F42i*u2*u4*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F44*F41i*u1*u4*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F44*F40i*u0*u4*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F43*F45i*u3*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F43*F44i*u3*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F43*F43i*u3**2*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F43*F42i*u2*u3*c2*f4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svp*svm*F43*F41i*u1*u3*c1*
     & f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F43*F40i*u0*u3*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F42*F45i*u2*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F42*F44i*u2*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F42*F43i*u2*u3*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F42*F42i*u2**2*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F42*F41i*u1*u2*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F42*F40i*u0*u2*c0*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F41*F45i*u1*u5*c5*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F41*F44i*u1*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F41*F43i*u1*u3*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F41*F42i*u1*u2*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F41*F41i*u1**2*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F41*F40i*u0*u1*c0*f4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svp*svm*F40*F45i*u0*u5*c5*
     & f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F40*F44i*u0*u4*c4*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F40*F43i*u0*u3*c3*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F40*F42i*u0*u2*c2*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F40*F41i*u0*u1*c1*f4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*svm*F40*F40i*u0**2*c0*f4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F35*F25i*u5**2*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F35*F24i*u4*u5*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F35*F23i*u3*u5*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F35*F22i*u2*u5*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F35*F21i*u1*u5*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F35*F20i*u0*u5*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F34*F25i*u4*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F34*F24i*u4**2*c4*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*svm*F34*F23i*u3*u4*c3*
     & f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F34*F22i*u2*u4*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F34*F21i*u1*u4*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F34*F20i*u0*u4*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F33*F25i*u3*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F33*F24i*u3*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F33*F23i*u3**2*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F33*F22i*u2*u3*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F33*F21i*u1*u3*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F33*F20i*u0*u3*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F32*F25i*u2*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F32*F24i*u2*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F32*F23i*u2*u3*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F32*F22i*u2**2*c2*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*svm*F32*F21i*u1*u2*c1*
     & f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F32*F20i*u0*u2*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F31*F25i*u1*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F31*F24i*u1*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F31*F23i*u1*u3*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F31*F22i*u1*u2*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F31*F21i*u1**2*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F31*F20i*u0*u1*c0*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F30*F25i*u0*u5*c5*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F30*F24i*u0*u4*c4*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F30*F23i*u0*u3*c3*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F30*F22i*u0*u2*c2*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F30*F21i*u0*u1*c1*f2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F30*F20i*u0**2*c0*f2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*svm*F25*F35i*u5**2*c5*
     & f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F25*F34i*u4*u5*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F25*F33i*u3*u5*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F25*F32i*u2*u5*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F25*F31i*u1*u5*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F25*F30i*u0*u5*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F24*F35i*u4*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F24*F34i*u4**2*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F24*F33i*u3*u4*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F24*F32i*u2*u4*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F24*F31i*u1*u4*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F24*F30i*u0*u4*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F23*F35i*u3*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F23*F34i*u3*u4*c4*f3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*svm*F23*F33i*u3**2*c3*
     & f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F23*F32i*u2*u3*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F23*F31i*u1*u3*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F23*F30i*u0*u3*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F22*F35i*u2*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F22*F34i*u2*u4*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F22*F33i*u2*u3*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F22*F32i*u2**2*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F22*F31i*u1*u2*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F22*F30i*u0*u2*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F21*F35i*u1*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F21*F34i*u1*u4*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F21*F33i*u1*u3*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F21*F32i*u1*u2*c2*f3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*svm*F21*F31i*u1**2*c1*
     & f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F21*F30i*u0*u1*c0*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F20*F35i*u0*u5*c5*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F20*F34i*u0*u4*c4*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F20*F33i*u0*u3*c3*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F20*F32i*u0*u2*c2*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F20*F31i*u0*u1*c1*f3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*svm*F20*F30i*u0**2*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F35*F35i*u5**2*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F35*F34i*u4*u5*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F35*F33i*u3*u5*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F35*F32i*u2*u5*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F35*F31i*u1*u5*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F35*F30i*u0*u5*c0*f3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F34*F35i*u4*u5*
     & c5*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F34*F34i*u4**2*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F34*F33i*u3*u4*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F34*F32i*u2*u4*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F34*F31i*u1*u4*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F34*F30i*u0*u4*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F33*F35i*u3*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F33*F34i*u3*u4*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F33*F33i*u3**2*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F33*F32i*u2*u3*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F33*F31i*u1*u3*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F33*F30i*u0*u3*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F32*F35i*u2*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F32*F34i*u2*u4*c4*f3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F32*F33i*u2*u3*
     & c3*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F32*F32i*u2**2*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F32*F31i*u1*u2*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F32*F30i*u0*u2*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F31*F35i*u1*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F31*F34i*u1*u4*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F31*F33i*u1*u3*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F31*F32i*u1*u2*c2*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F31*F31i*u1**2*c1*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F31*F30i*u0*u1*c0*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F30*F35i*u0*u5*c5*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F30*F34i*u0*u4*c4*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F30*F33i*u0*u3*c3*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F30*F32i*u0*u2*c2*f3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F30*F31i*u0*u1*
     & c1*f3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*svm*F30*F30i*u0**2*c0*f3*w1*w2
     &  + 4.D0*PC_q**3*svp*svm*F45*F45r*u5**2*c5*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F45*F45r*u5**2*c5*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F45*F44r*u4*u5*c4*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F45*F44r*u4*u5*c4*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F45*F43r*u3*u5*c3*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F45*F43r*u3*u5*c3*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F45*F42r*u2*u5*c2*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F45*F42r*u2*u5*c2*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F45*F41r*u1*u5*c1*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F45*F41r*u1*u5*c1*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F45*F40r*u0*u5*c0*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F45*F40r*u0*u5*c0*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svp*svm*F44*F45r*u4*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F44*F45r*u4*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F44*F44r*u4**2*c4*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F44*F44r*u4**2*c4*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F44*F43r*u3*u4*c3*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F44*F43r*u3*u4*c3*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F44*F42r*u2*u4*c2*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F44*F42r*u2*u4*c2*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F44*F41r*u1*u4*c1*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F44*F41r*u1*u4*c1*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F44*F40r*u0*u4*c0*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F44*F40r*u0*u4*c0*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F43*F45r*u3*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F43*F45r*u3*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F43*F44r*u3*u4*c4*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*svm*F43*F44r*u3*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F43*F43r*u3**2*c3*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F43*F43r*u3**2*c3*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F43*F42r*u2*u3*c2*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F43*F42r*u2*u3*c2*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F43*F41r*u1*u3*c1*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F43*F41r*u1*u3*c1*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F43*F40r*u0*u3*c0*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F43*F40r*u0*u3*c0*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F42*F45r*u2*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F42*F45r*u2*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F42*F44r*u2*u4*c4*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F42*F44r*u2*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F42*F43r*u2*u3*c3*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F42*F43r*u2*u3*c3*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svp*svm*F42*F42r*u2**2*c2*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F42*F42r*u2**2*c2*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F42*F41r*u1*u2*c1*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F42*F41r*u1*u2*c1*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F42*F40r*u0*u2*c0*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F42*F40r*u0*u2*c0*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F41*F45r*u1*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F41*F45r*u1*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F41*F44r*u1*u4*c4*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F41*F44r*u1*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F41*F43r*u1*u3*c3*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F41*F43r*u1*u3*c3*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F41*F42r*u1*u2*c2*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F41*F42r*u1*u2*c2*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F41*F41r*u1**2*c1*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*svm*F41*F41r*u1**2*c1*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F41*F40r*u0*u1*c0*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F41*F40r*u0*u1*c0*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F40*F45r*u0*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F40*F45r*u0*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F40*F44r*u0*u4*c4*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F40*F44r*u0*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F40*F43r*u0*u3*c3*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F40*F43r*u0*u3*c3*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F40*F42r*u0*u2*c2*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F40*F42r*u0*u2*c2*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F40*F41r*u0*u1*c1*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F40*F41r*u0*u1*c1*f4*w1
     &  + 4.D0*PC_q**3*svp*svm*F40*F40r*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**3*svp*svm*F40*F40r*u0**2*c0*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*svp*F45*F35r*u5**2*c5*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F45*F34r*u4*u5*c4*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F45*F33r*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F45*F32r*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F45*F31r*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F45*F30r*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F44*F35r*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F44*F34r*u4**2*c4*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F44*F33r*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F44*F32r*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F44*F31r*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F44*F30r*u0*u4*c0*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F43*F35r*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F43*F34r*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F43*F33r*u3**2*c3*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*svp*F43*F32r*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F43*F31r*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F43*F30r*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F42*F35r*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F42*F34r*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F42*F33r*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F42*F32r*u2**2*c2*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F42*F31r*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F42*F30r*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F41*F35r*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F41*F34r*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F41*F33r*u1*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F41*F32r*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F41*F31r*u1**2*c1*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F41*F30r*u0*u1*c0*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*svp*F40*F35r*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F40*F34r*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F40*F33r*u0*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F40*F32r*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F40*F31r*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F40*F30r*u0**2*c0*f3*w1
     &  - 4.D0*PC_q**3*ssm*svp*F35*F45r*u5**2*c5*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F35*F44r*u4*u5*c4*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F35*F43r*u3*u5*c3*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F35*F42r*u2*u5*c2*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F35*F41r*u1*u5*c1*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F35*F40r*u0*u5*c0*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F34*F45r*u4*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F34*F44r*u4**2*c4*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F34*F43r*u3*u4*c3*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*svp*F34*F42r*u2*u4*c2*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F34*F41r*u1*u4*c1*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F34*F40r*u0*u4*c0*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F33*F45r*u3*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F33*F44r*u3*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F33*F43r*u3**2*c3*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F33*F42r*u2*u3*c2*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F33*F41r*u1*u3*c1*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F33*F40r*u0*u3*c0*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F32*F45r*u2*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F32*F44r*u2*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F32*F43r*u2*u3*c3*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F32*F42r*u2**2*c2*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F32*F41r*u1*u2*c1*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F32*F40r*u0*u2*c0*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*svp*F31*F45r*u1*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F31*F44r*u1*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F31*F43r*u1*u3*c3*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F31*F42r*u1*u2*c2*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F31*F41r*u1**2*c1*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F31*F40r*u0*u1*c0*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F30*F45r*u0*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F30*F44r*u0*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F30*F43r*u0*u3*c3*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F30*F42r*u0*u2*c2*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F30*F41r*u0*u1*c1*f4*w1
     &  - 4.D0*PC_q**3*ssm*svp*F30*F40r*u0**2*c0*f4*w1
     &  + 4.D0*PC_q**3*ssp*svm*F45*F35r*u5**2*c5*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F45*F34r*u4*u5*c4*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F45*F33r*u3*u5*c3*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*svm*F45*F32r*u2*u5*c2*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F45*F31r*u1*u5*c1*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F45*F30r*u0*u5*c0*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F44*F35r*u4*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F44*F34r*u4**2*c4*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F44*F33r*u3*u4*c3*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F44*F32r*u2*u4*c2*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F44*F31r*u1*u4*c1*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F44*F30r*u0*u4*c0*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F43*F35r*u3*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F43*F34r*u3*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F43*F33r*u3**2*c3*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F43*F32r*u2*u3*c2*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F43*F31r*u1*u3*c1*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F43*F30r*u0*u3*c0*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*svm*F42*F35r*u2*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F42*F34r*u2*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F42*F33r*u2*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F42*F32r*u2**2*c2*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F42*F31r*u1*u2*c1*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F42*F30r*u0*u2*c0*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F41*F35r*u1*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F41*F34r*u1*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F41*F33r*u1*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F41*F32r*u1*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F41*F31r*u1**2*c1*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F41*F30r*u0*u1*c0*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F40*F35r*u0*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F40*F34r*u0*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F40*F33r*u0*u3*c3*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*svm*F40*F32r*u0*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F40*F31r*u0*u1*c1*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F40*F30r*u0**2*c0*f3*w2
     &  + 4.D0*PC_q**3*ssp*svm*F35*F45r*u5**2*c5*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F35*F44r*u4*u5*c4*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F35*F43r*u3*u5*c3*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F35*F42r*u2*u5*c2*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F35*F41r*u1*u5*c1*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F35*F40r*u0*u5*c0*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F34*F45r*u4*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F34*F44r*u4**2*c4*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F34*F43r*u3*u4*c3*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F34*F42r*u2*u4*c2*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F34*F41r*u1*u4*c1*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F34*F40r*u0*u4*c0*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*svm*F33*F45r*u3*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F33*F44r*u3*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F33*F43r*u3**2*c3*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F33*F42r*u2*u3*c2*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F33*F41r*u1*u3*c1*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F33*F40r*u0*u3*c0*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F32*F45r*u2*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F32*F44r*u2*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F32*F43r*u2*u3*c3*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F32*F42r*u2**2*c2*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F32*F41r*u1*u2*c1*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F32*F40r*u0*u2*c0*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F31*F45r*u1*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F31*F44r*u1*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F31*F43r*u1*u3*c3*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*svm*F31*F42r*u1*u2*c2*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F31*F41r*u1**2*c1*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F31*F40r*u0*u1*c0*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F30*F45r*u0*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F30*F44r*u0*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F30*F43r*u0*u3*c3*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F30*F42r*u0*u2*c2*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F30*F41r*u0*u1*c1*f4*w2
     &  + 4.D0*PC_q**3*ssp*svm*F30*F40r*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**3*qq*svp*svm*F35*F35r*u5**2*c5*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F35*F35r*u5**2*c5*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F35*F34r*u4*u5*c4*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F35*F34r*u4*u5*c4*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F35*F33r*u3*u5*c3*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F35*F33r*u3*u5*c3*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*qq*svp*svm*F35*F32r*u2*u5*c2*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F35*F32r*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F35*F31r*u1*u5*c1*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F35*F31r*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F35*F30r*u0*u5*c0*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F35*F30r*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F34*F35r*u4*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F34*F35r*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F34*F34r*u4**2*c4*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F34*F34r*u4**2*c4*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F34*F33r*u3*u4*c3*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F34*F33r*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F34*F32r*u2*u4*c2*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F34*F32r*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F34*F31r*u1*u4*c1*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*svm*F34*F31r*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F34*F30r*u0*u4*c0*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F34*F30r*u0*u4*c0*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F33*F35r*u3*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F33*F35r*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F33*F34r*u3*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F33*F34r*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F33*F33r*u3**2*c3*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F33*F33r*u3**2*c3*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F33*F32r*u2*u3*c2*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F33*F32r*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F33*F31r*u1*u3*c1*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F33*F31r*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F33*F30r*u0*u3*c0*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F33*F30r*u0*u3*c0*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*qq*svp*svm*F32*F35r*u2*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F32*F35r*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F32*F34r*u2*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F32*F34r*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F32*F33r*u2*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F32*F33r*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F32*F32r*u2**2*c2*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F32*F32r*u2**2*c2*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F32*F31r*u1*u2*c1*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F32*F31r*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F32*F30r*u0*u2*c0*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F32*F30r*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F31*F35r*u1*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F31*F35r*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F31*F34r*u1*u4*c4*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*svm*F31*F34r*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F31*F33r*u1*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F31*F33r*u1*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F31*F32r*u1*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F31*F32r*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F31*F31r*u1**2*c1*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F31*F31r*u1**2*c1*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F31*F30r*u0*u1*c0*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F31*F30r*u0*u1*c0*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F30*F35r*u0*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F30*F35r*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F30*F34r*u0*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F30*F34r*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F30*F33r*u0*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F30*F33r*u0*u3*c3*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*qq*svp*svm*F30*F32r*u0*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F30*F32r*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F30*F31r*u0*u1*c1*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F30*F31r*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q**3*qq*svp*svm*F30*F30r*u0**2*c0*f3*w2
     &  + 4.D0*PC_q**3*qq*svp*svm*F30*F30r*u0**2*c0*f3*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F45*F45i*u5**2*c5*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F45*F45i*u5**2*c5*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F45*F44i*u4*u5*c4*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F45*F44i*u4*u5*c4*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F45*F43i*u3*u5*c3*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F45*F43i*u3*u5*c3*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F45*F42i*u2*u5*c2*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F45*F42i*u2*u5*c2*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F45*F41i*u1*u5*c1*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*svm*F45*F41i*u1*u5*c1*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F45*F40i*u0*u5*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F45*F40i*u0*u5*c0*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F44*F45i*u4*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F44*F45i*u4*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F44*F44i*u4**2*c4*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F44*F44i*u4**2*c4*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F44*F43i*u3*u4*c3*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F44*F43i*u3*u4*c3*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F44*F42i*u2*u4*c2*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F44*F42i*u2*u4*c2*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F44*F41i*u1*u4*c1*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F44*F41i*u1*u4*c1*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F44*F40i*u0*u4*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F44*F40i*u0*u4*c0*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svp*svm*F43*F45i*u3*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F43*F45i*u3*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F43*F44i*u3*u4*c4*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F43*F44i*u3*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F43*F43i*u3**2*c3*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F43*F43i*u3**2*c3*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F43*F42i*u2*u3*c2*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F43*F42i*u2*u3*c2*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F43*F41i*u1*u3*c1*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F43*F41i*u1*u3*c1*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F43*F40i*u0*u3*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F43*F40i*u0*u3*c0*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F42*F45i*u2*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F42*F45i*u2*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F42*F44i*u2*u4*c4*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*svm*F42*F44i*u2*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F42*F43i*u2*u3*c3*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F42*F43i*u2*u3*c3*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F42*F42i*u2**2*c2*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F42*F42i*u2**2*c2*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F42*F41i*u1*u2*c1*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F42*F41i*u1*u2*c1*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F42*F40i*u0*u2*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F42*F40i*u0*u2*c0*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F41*F45i*u1*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F41*F45i*u1*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F41*F44i*u1*u4*c4*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F41*F44i*u1*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F41*F43i*u1*u3*c3*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F41*F43i*u1*u3*c3*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svp*svm*F41*F42i*u1*u2*c2*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F41*F42i*u1*u2*c2*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F41*F41i*u1**2*c1*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F41*F41i*u1**2*c1*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F41*F40i*u0*u1*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F41*F40i*u0*u1*c0*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F40*F45i*u0*u5*c5*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F40*F45i*u0*u5*c5*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F40*F44i*u0*u4*c4*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F40*F44i*u0*u4*c4*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F40*F43i*u0*u3*c3*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F40*F43i*u0*u3*c3*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F40*F42i*u0*u2*c2*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F40*F42i*u0*u2*c2*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F40*F41i*u0*u1*c1*f4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*svm*F40*F41i*u0*u1*c1*f4*w1
     &  + 4.D0*PC_q**3*i_*svp*svm*F40*F40i*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*svp*svm*F40*F40i*u0**2*c0*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F45*F35i*u5**2*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F45*F34i*u4*u5*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F45*F33i*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F45*F32i*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F45*F31i*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F45*F30i*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F44*F35i*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F44*F34i*u4**2*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F44*F33i*u3*u4*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F44*F32i*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F44*F31i*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F44*F30i*u0*u4*c0*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*svp*F43*F35i*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F43*F34i*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F43*F33i*u3**2*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F43*F32i*u2*u3*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F43*F31i*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F43*F30i*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F42*F35i*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F42*F34i*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F42*F33i*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F42*F32i*u2**2*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F42*F31i*u1*u2*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F42*F30i*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F41*F35i*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F41*F34i*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F41*F33i*u1*u3*c3*f3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*svp*F41*F32i*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F41*F31i*u1**2*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F41*F30i*u0*u1*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F40*F35i*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F40*F34i*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F40*F33i*u0*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F40*F32i*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F40*F31i*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F40*F30i*u0**2*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F35*F45i*u5**2*c5*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F35*F44i*u4*u5*c4*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F35*F43i*u3*u5*c3*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F35*F42i*u2*u5*c2*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F35*F41i*u1*u5*c1*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F35*F40i*u0*u5*c0*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*svp*F34*F45i*u4*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F34*F44i*u4**2*c4*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F34*F43i*u3*u4*c3*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F34*F42i*u2*u4*c2*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F34*F41i*u1*u4*c1*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F34*F40i*u0*u4*c0*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F33*F45i*u3*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F33*F44i*u3*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F33*F43i*u3**2*c3*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F33*F42i*u2*u3*c2*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F33*F41i*u1*u3*c1*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F33*F40i*u0*u3*c0*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F32*F45i*u2*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F32*F44i*u2*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F32*F43i*u2*u3*c3*f4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*svp*F32*F42i*u2**2*c2*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F32*F41i*u1*u2*c1*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F32*F40i*u0*u2*c0*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F31*F45i*u1*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F31*F44i*u1*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F31*F43i*u1*u3*c3*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F31*F42i*u1*u2*c2*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F31*F41i*u1**2*c1*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F31*F40i*u0*u1*c0*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F30*F45i*u0*u5*c5*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F30*F44i*u0*u4*c4*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F30*F43i*u0*u3*c3*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F30*F42i*u0*u2*c2*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F30*F41i*u0*u1*c1*f4*w1
     &  - 4.D0*PC_q**3*i_*ssm*svp*F30*F40i*u0**2*c0*f4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*svm*F45*F35i*u5**2*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F45*F34i*u4*u5*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F45*F33i*u3*u5*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F45*F32i*u2*u5*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F45*F31i*u1*u5*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F45*F30i*u0*u5*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F44*F35i*u4*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F44*F34i*u4**2*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F44*F33i*u3*u4*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F44*F32i*u2*u4*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F44*F31i*u1*u4*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F44*F30i*u0*u4*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F43*F35i*u3*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F43*F34i*u3*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F43*F33i*u3**2*c3*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*svm*F43*F32i*u2*u3*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F43*F31i*u1*u3*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F43*F30i*u0*u3*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F42*F35i*u2*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F42*F34i*u2*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F42*F33i*u2*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F42*F32i*u2**2*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F42*F31i*u1*u2*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F42*F30i*u0*u2*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F41*F35i*u1*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F41*F34i*u1*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F41*F33i*u1*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F41*F32i*u1*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F41*F31i*u1**2*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F41*F30i*u0*u1*c0*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*svm*F40*F35i*u0*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F40*F34i*u0*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F40*F33i*u0*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F40*F32i*u0*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F40*F31i*u0*u1*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F40*F30i*u0**2*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F35*F45i*u5**2*c5*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F35*F44i*u4*u5*c4*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F35*F43i*u3*u5*c3*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F35*F42i*u2*u5*c2*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F35*F41i*u1*u5*c1*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F35*F40i*u0*u5*c0*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F34*F45i*u4*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F34*F44i*u4**2*c4*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F34*F43i*u3*u4*c3*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*svm*F34*F42i*u2*u4*c2*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F34*F41i*u1*u4*c1*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F34*F40i*u0*u4*c0*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F33*F45i*u3*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F33*F44i*u3*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F33*F43i*u3**2*c3*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F33*F42i*u2*u3*c2*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F33*F41i*u1*u3*c1*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F33*F40i*u0*u3*c0*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F32*F45i*u2*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F32*F44i*u2*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F32*F43i*u2*u3*c3*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F32*F42i*u2**2*c2*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F32*F41i*u1*u2*c1*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F32*F40i*u0*u2*c0*f4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*svm*F31*F45i*u1*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F31*F44i*u1*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F31*F43i*u1*u3*c3*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F31*F42i*u1*u2*c2*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F31*F41i*u1**2*c1*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F31*F40i*u0*u1*c0*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F30*F45i*u0*u5*c5*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F30*F44i*u0*u4*c4*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F30*F43i*u0*u3*c3*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F30*F42i*u0*u2*c2*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F30*F41i*u0*u1*c1*f4*w2
     &  + 4.D0*PC_q**3*i_*ssp*svm*F30*F40i*u0**2*c0*f4*w2
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F35*F35i*u5**2*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F35*F35i*u5**2*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F35*F34i*u4*u5*c4*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*svm*F35*F34i*u4*u5*c4*f3
     & *w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F35*F33i*u3*u5*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F35*F33i*u3*u5*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F35*F32i*u2*u5*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F35*F32i*u2*u5*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F35*F31i*u1*u5*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F35*F31i*u1*u5*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F35*F30i*u0*u5*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F35*F30i*u0*u5*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F34*F35i*u4*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F34*F35i*u4*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F34*F34i*u4**2*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F34*F34i*u4**2*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F34*F33i*u3*u4*c3*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*svm*F34*F33i*u3*u4*c3*f3
     & *w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F34*F32i*u2*u4*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F34*F32i*u2*u4*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F34*F31i*u1*u4*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F34*F31i*u1*u4*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F34*F30i*u0*u4*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F34*F30i*u0*u4*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F33*F35i*u3*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F33*F35i*u3*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F33*F34i*u3*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F33*F34i*u3*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F33*F33i*u3**2*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F33*F33i*u3**2*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F33*F32i*u2*u3*c2*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*svm*F33*F32i*u2*u3*c2*f3
     & *w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F33*F31i*u1*u3*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F33*F31i*u1*u3*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F33*F30i*u0*u3*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F33*F30i*u0*u3*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F32*F35i*u2*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F32*F35i*u2*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F32*F34i*u2*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F32*F34i*u2*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F32*F33i*u2*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F32*F33i*u2*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F32*F32i*u2**2*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F32*F32i*u2**2*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F32*F31i*u1*u2*c1*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*svm*F32*F31i*u1*u2*c1*f3
     & *w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F32*F30i*u0*u2*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F32*F30i*u0*u2*c0*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F31*F35i*u1*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F31*F35i*u1*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F31*F34i*u1*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F31*F34i*u1*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F31*F33i*u1*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F31*F33i*u1*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F31*F32i*u1*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F31*F32i*u1*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F31*F31i*u1**2*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F31*F31i*u1**2*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F31*F30i*u0*u1*c0*f3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*svm*F31*F30i*u0*u1*c0*f3
     & *w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F30*F35i*u0*u5*c5*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F30*F35i*u0*u5*c5*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F30*F34i*u0*u4*c4*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F30*F34i*u0*u4*c4*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F30*F33i*u0*u3*c3*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F30*F33i*u0*u3*c3*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F30*F32i*u0*u2*c2*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F30*F32i*u0*u2*c2*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F30*F31i*u0*u1*c1*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F30*F31i*u0*u1*c1*f3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*svm*F30*F30i*u0**2*c0*f3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*svm*F30*F30i*u0**2*c0*f3*w1
     &  - 8.D0*PC_q**4*svp*svm*F35*F35r*u5**2*c5*f3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svp*svm*F35*F34r*u4*u5*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F35*F33r*u3*u5*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F35*F32r*u2*u5*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F35*F31r*u1*u5*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F35*F30r*u0*u5*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F34*F35r*u4*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F34*F34r*u4**2*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F34*F33r*u3*u4*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F34*F32r*u2*u4*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F34*F31r*u1*u4*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F34*F30r*u0*u4*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F33*F35r*u3*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F33*F34r*u3*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F33*F33r*u3**2*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F33*F32r*u2*u3*c2*f3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svp*svm*F33*F31r*u1*u3*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F33*F30r*u0*u3*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F32*F35r*u2*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F32*F34r*u2*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F32*F33r*u2*u3*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F32*F32r*u2**2*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F32*F31r*u1*u2*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F32*F30r*u0*u2*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F31*F35r*u1*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F31*F34r*u1*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F31*F33r*u1*u3*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F31*F32r*u1*u2*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F31*F31r*u1**2*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F31*F30r*u0*u1*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F30*F35r*u0*u5*c5*f3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svp*svm*F30*F34r*u0*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F30*F33r*u0*u3*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F30*F32r*u0*u2*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F30*F31r*u0*u1*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*svp*svm*F30*F30r*u0**2*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F35*F35i*u5**2*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F35*F34i*u4*u5*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F35*F33i*u3*u5*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F35*F32i*u2*u5*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F35*F31i*u1*u5*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F35*F30i*u0*u5*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F34*F35i*u4*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F34*F34i*u4**2*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F34*F33i*u3*u4*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F34*F32i*u2*u4*c2*f3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svp*svm*F34*F31i*u1*u4*c1*f3*w1
     & *w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F34*F30i*u0*u4*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F33*F35i*u3*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F33*F34i*u3*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F33*F33i*u3**2*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F33*F32i*u2*u3*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F33*F31i*u1*u3*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F33*F30i*u0*u3*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F32*F35i*u2*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F32*F34i*u2*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F32*F33i*u2*u3*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F32*F32i*u2**2*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F32*F31i*u1*u2*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F32*F30i*u0*u2*c0*f3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svp*svm*F31*F35i*u1*u5*c5*f3*w1
     & *w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F31*F34i*u1*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F31*F33i*u1*u3*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F31*F32i*u1*u2*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F31*F31i*u1**2*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F31*F30i*u0*u1*c0*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F30*F35i*u0*u5*c5*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F30*F34i*u0*u4*c4*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F30*F33i*u0*u3*c3*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F30*F32i*u0*u2*c2*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F30*F31i*u0*u1*c1*f3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*svm*F30*F30i*u0**2*c0*f3*w1*w2
     &

   
       if(print.eqv..true.)then
        write(2,*)"traza-variables_____________________________________"       
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
       write(2,*)"f1,f2,f3,f4",f1,f2,f3,f4
       write(2,*)"c0,c1,c2,c3,c4,c5",c0,c1,c2,c3,c4,c5
       write(2,*)"z1",z1
       write(2,*)"uh(1),uh(2),uh(3),uh(4)",uh(1),uh(2),uh(3),uh(4)
       write(2,*)"pcv(i)"
       write(2,*)"pcv(1),pcv(2),pcv(3),pcv(4)",
     . pcv(1),pcv(2),pcv(3),pcv(4)
       write(2,*)"pc2,pc_q,traza1",pc2,pc_q,traza1
       write(2,*)"ssp,ssm,svp,svm",ssp,ssm,svp,svm
       write(2,*)"dssp,dssm,dsvp,dsvm",dssp,dssm,dsvp,dsvm
       write(2,*)"pc_uh,q_uh",pc_uh,q_uh
       print = .false.
       else
       endif


       return
       end
