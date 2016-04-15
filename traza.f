       subroutine traza(print,w1,z1,qv,pcv,ssp,ssm,svp,
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

      traza1 = - 4.D0*ssm*dssp*F15*F15r*c5*f1*u5**2
     &  - 4.D0*ssm*dssp*F15*F14r*c4*f1*u4*u5
     &  - 4.D0*ssm*dssp*F15*F13r*c3*f1*u3*u5
     &  - 4.D0*ssm*dssp*F15*F12r*c2*f1*u2*u5
     &  - 4.D0*ssm*dssp*F15*F11r*c1*f1*u1*u5
     &  - 4.D0*ssm*dssp*F15*F10r*c0*f1*u0*u5
     &  - 4.D0*ssm*dssp*F14*F15r*c5*f1*u4*u5
     &  - 4.D0*ssm*dssp*F14*F14r*c4*f1*u4**2
     &  - 4.D0*ssm*dssp*F14*F13r*c3*f1*u3*u4
     &  - 4.D0*ssm*dssp*F14*F12r*c2*f1*u2*u4
     &  - 4.D0*ssm*dssp*F14*F11r*c1*f1*u1*u4
     &  - 4.D0*ssm*dssp*F14*F10r*c0*f1*u0*u4
     &  - 4.D0*ssm*dssp*F13*F15r*c5*f1*u3*u5
     &  - 4.D0*ssm*dssp*F13*F14r*c4*f1*u3*u4
     &  - 4.D0*ssm*dssp*F13*F13r*c3*f1*u3**2
     &
      traza1 = traza1 - 4.D0*ssm*dssp*F13*F12r*c2*f1*u2*u3
     &  - 4.D0*ssm*dssp*F13*F11r*c1*f1*u1*u3
     &  - 4.D0*ssm*dssp*F13*F10r*c0*f1*u0*u3
     &  - 4.D0*ssm*dssp*F12*F15r*c5*f1*u2*u5
     &  - 4.D0*ssm*dssp*F12*F14r*c4*f1*u2*u4
     &  - 4.D0*ssm*dssp*F12*F13r*c3*f1*u2*u3
     &  - 4.D0*ssm*dssp*F12*F12r*c2*f1*u2**2
     &  - 4.D0*ssm*dssp*F12*F11r*c1*f1*u1*u2
     &  - 4.D0*ssm*dssp*F12*F10r*c0*f1*u0*u2
     &  - 4.D0*ssm*dssp*F11*F15r*c5*f1*u1*u5
     &  - 4.D0*ssm*dssp*F11*F14r*c4*f1*u1*u4
     &  - 4.D0*ssm*dssp*F11*F13r*c3*f1*u1*u3
     &  - 4.D0*ssm*dssp*F11*F12r*c2*f1*u1*u2
     &  - 4.D0*ssm*dssp*F11*F11r*c1*f1*u1**2
     &  - 4.D0*ssm*dssp*F11*F10r*c0*f1*u0*u1
     &
      traza1 = traza1 - 4.D0*ssm*dssp*F10*F15r*c5*f1*u0*u5
     &  - 4.D0*ssm*dssp*F10*F14r*c4*f1*u0*u4
     &  - 4.D0*ssm*dssp*F10*F13r*c3*f1*u0*u3
     &  - 4.D0*ssm*dssp*F10*F12r*c2*f1*u0*u2
     &  - 4.D0*ssm*dssp*F10*F11r*c1*f1*u0*u1
     &  - 4.D0*ssm*dssp*F10*F10r*c0*f1*u0**2
     &  - 4.D0*ssp*dssm*F15*F15r*c5*f1*u5**2
     &  - 4.D0*ssp*dssm*F15*F14r*c4*f1*u4*u5
     &  - 4.D0*ssp*dssm*F15*F13r*c3*f1*u3*u5
     &  - 4.D0*ssp*dssm*F15*F12r*c2*f1*u2*u5
     &  - 4.D0*ssp*dssm*F15*F11r*c1*f1*u1*u5
     &  - 4.D0*ssp*dssm*F15*F10r*c0*f1*u0*u5
     &  - 4.D0*ssp*dssm*F14*F15r*c5*f1*u4*u5
     &  - 4.D0*ssp*dssm*F14*F14r*c4*f1*u4**2
     &  - 4.D0*ssp*dssm*F14*F13r*c3*f1*u3*u4
     &
      traza1 = traza1 - 4.D0*ssp*dssm*F14*F12r*c2*f1*u2*u4
     &  - 4.D0*ssp*dssm*F14*F11r*c1*f1*u1*u4
     &  - 4.D0*ssp*dssm*F14*F10r*c0*f1*u0*u4
     &  - 4.D0*ssp*dssm*F13*F15r*c5*f1*u3*u5
     &  - 4.D0*ssp*dssm*F13*F14r*c4*f1*u3*u4
     &  - 4.D0*ssp*dssm*F13*F13r*c3*f1*u3**2
     &  - 4.D0*ssp*dssm*F13*F12r*c2*f1*u2*u3
     &  - 4.D0*ssp*dssm*F13*F11r*c1*f1*u1*u3
     &  - 4.D0*ssp*dssm*F13*F10r*c0*f1*u0*u3
     &  - 4.D0*ssp*dssm*F12*F15r*c5*f1*u2*u5
     &  - 4.D0*ssp*dssm*F12*F14r*c4*f1*u2*u4
     &  - 4.D0*ssp*dssm*F12*F13r*c3*f1*u2*u3
     &  - 4.D0*ssp*dssm*F12*F12r*c2*f1*u2**2
     &  - 4.D0*ssp*dssm*F12*F11r*c1*f1*u1*u2
     &  - 4.D0*ssp*dssm*F12*F10r*c0*f1*u0*u2
     &
      traza1 = traza1 - 4.D0*ssp*dssm*F11*F15r*c5*f1*u1*u5
     &  - 4.D0*ssp*dssm*F11*F14r*c4*f1*u1*u4
     &  - 4.D0*ssp*dssm*F11*F13r*c3*f1*u1*u3
     &  - 4.D0*ssp*dssm*F11*F12r*c2*f1*u1*u2
     &  - 4.D0*ssp*dssm*F11*F11r*c1*f1*u1**2
     &  - 4.D0*ssp*dssm*F11*F10r*c0*f1*u0*u1
     &  - 4.D0*ssp*dssm*F10*F15r*c5*f1*u0*u5
     &  - 4.D0*ssp*dssm*F10*F14r*c4*f1*u0*u4
     &  - 4.D0*ssp*dssm*F10*F13r*c3*f1*u0*u3
     &  - 4.D0*ssp*dssm*F10*F12r*c2*f1*u0*u2
     &  - 4.D0*ssp*dssm*F10*F11r*c1*f1*u0*u1
     &  - 4.D0*ssp*dssm*F10*F10r*c0*f1*u0**2
     &  - 4.D0*qq*svm*dsvp*F15*F15r*c5*f1*u5**2
     &  - 4.D0*qq*svm*dsvp*F15*F14r*c4*f1*u4*u5
     &  - 4.D0*qq*svm*dsvp*F15*F13r*c3*f1*u3*u5
     &
      traza1 = traza1 - 4.D0*qq*svm*dsvp*F15*F12r*c2*f1*u2*u5
     &  - 4.D0*qq*svm*dsvp*F15*F11r*c1*f1*u1*u5
     &  - 4.D0*qq*svm*dsvp*F15*F10r*c0*f1*u0*u5
     &  - 4.D0*qq*svm*dsvp*F14*F15r*c5*f1*u4*u5
     &  - 4.D0*qq*svm*dsvp*F14*F14r*c4*f1*u4**2
     &  - 4.D0*qq*svm*dsvp*F14*F13r*c3*f1*u3*u4
     &  - 4.D0*qq*svm*dsvp*F14*F12r*c2*f1*u2*u4
     &  - 4.D0*qq*svm*dsvp*F14*F11r*c1*f1*u1*u4
     &  - 4.D0*qq*svm*dsvp*F14*F10r*c0*f1*u0*u4
     &  - 4.D0*qq*svm*dsvp*F13*F15r*c5*f1*u3*u5
     &  - 4.D0*qq*svm*dsvp*F13*F14r*c4*f1*u3*u4
     &  - 4.D0*qq*svm*dsvp*F13*F13r*c3*f1*u3**2
     &  - 4.D0*qq*svm*dsvp*F13*F12r*c2*f1*u2*u3
     &  - 4.D0*qq*svm*dsvp*F13*F11r*c1*f1*u1*u3
     &  - 4.D0*qq*svm*dsvp*F13*F10r*c0*f1*u0*u3
     &
      traza1 = traza1 - 4.D0*qq*svm*dsvp*F12*F15r*c5*f1*u2*u5
     &  - 4.D0*qq*svm*dsvp*F12*F14r*c4*f1*u2*u4
     &  - 4.D0*qq*svm*dsvp*F12*F13r*c3*f1*u2*u3
     &  - 4.D0*qq*svm*dsvp*F12*F12r*c2*f1*u2**2
     &  - 4.D0*qq*svm*dsvp*F12*F11r*c1*f1*u1*u2
     &  - 4.D0*qq*svm*dsvp*F12*F10r*c0*f1*u0*u2
     &  - 4.D0*qq*svm*dsvp*F11*F15r*c5*f1*u1*u5
     &  - 4.D0*qq*svm*dsvp*F11*F14r*c4*f1*u1*u4
     &  - 4.D0*qq*svm*dsvp*F11*F13r*c3*f1*u1*u3
     &  - 4.D0*qq*svm*dsvp*F11*F12r*c2*f1*u1*u2
     &  - 4.D0*qq*svm*dsvp*F11*F11r*c1*f1*u1**2
     &  - 4.D0*qq*svm*dsvp*F11*F10r*c0*f1*u0*u1
     &  - 4.D0*qq*svm*dsvp*F10*F15r*c5*f1*u0*u5
     &  - 4.D0*qq*svm*dsvp*F10*F14r*c4*f1*u0*u4
     &  - 4.D0*qq*svm*dsvp*F10*F13r*c3*f1*u0*u3
     &
      traza1 = traza1 - 4.D0*qq*svm*dsvp*F10*F12r*c2*f1*u0*u2
     &  - 4.D0*qq*svm*dsvp*F10*F11r*c1*f1*u0*u1
     &  - 4.D0*qq*svm*dsvp*F10*F10r*c0*f1*u0**2
     &  - 4.D0*qq*svp*dsvm*F15*F15r*c5*f1*u5**2
     &  - 4.D0*qq*svp*dsvm*F15*F14r*c4*f1*u4*u5
     &  - 4.D0*qq*svp*dsvm*F15*F13r*c3*f1*u3*u5
     &  - 4.D0*qq*svp*dsvm*F15*F12r*c2*f1*u2*u5
     &  - 4.D0*qq*svp*dsvm*F15*F11r*c1*f1*u1*u5
     &  - 4.D0*qq*svp*dsvm*F15*F10r*c0*f1*u0*u5
     &  - 4.D0*qq*svp*dsvm*F14*F15r*c5*f1*u4*u5
     &  - 4.D0*qq*svp*dsvm*F14*F14r*c4*f1*u4**2
     &  - 4.D0*qq*svp*dsvm*F14*F13r*c3*f1*u3*u4
     &  - 4.D0*qq*svp*dsvm*F14*F12r*c2*f1*u2*u4
     &  - 4.D0*qq*svp*dsvm*F14*F11r*c1*f1*u1*u4
     &  - 4.D0*qq*svp*dsvm*F14*F10r*c0*f1*u0*u4
     &
      traza1 = traza1 - 4.D0*qq*svp*dsvm*F13*F15r*c5*f1*u3*u5
     &  - 4.D0*qq*svp*dsvm*F13*F14r*c4*f1*u3*u4
     &  - 4.D0*qq*svp*dsvm*F13*F13r*c3*f1*u3**2
     &  - 4.D0*qq*svp*dsvm*F13*F12r*c2*f1*u2*u3
     &  - 4.D0*qq*svp*dsvm*F13*F11r*c1*f1*u1*u3
     &  - 4.D0*qq*svp*dsvm*F13*F10r*c0*f1*u0*u3
     &  - 4.D0*qq*svp*dsvm*F12*F15r*c5*f1*u2*u5
     &  - 4.D0*qq*svp*dsvm*F12*F14r*c4*f1*u2*u4
     &  - 4.D0*qq*svp*dsvm*F12*F13r*c3*f1*u2*u3
     &  - 4.D0*qq*svp*dsvm*F12*F12r*c2*f1*u2**2
     &  - 4.D0*qq*svp*dsvm*F12*F11r*c1*f1*u1*u2
     &  - 4.D0*qq*svp*dsvm*F12*F10r*c0*f1*u0*u2
     &  - 4.D0*qq*svp*dsvm*F11*F15r*c5*f1*u1*u5
     &  - 4.D0*qq*svp*dsvm*F11*F14r*c4*f1*u1*u4
     &  - 4.D0*qq*svp*dsvm*F11*F13r*c3*f1*u1*u3
     &
      traza1 = traza1 - 4.D0*qq*svp*dsvm*F11*F12r*c2*f1*u1*u2
     &  - 4.D0*qq*svp*dsvm*F11*F11r*c1*f1*u1**2
     &  - 4.D0*qq*svp*dsvm*F11*F10r*c0*f1*u0*u1
     &  - 4.D0*qq*svp*dsvm*F10*F15r*c5*f1*u0*u5
     &  - 4.D0*qq*svp*dsvm*F10*F14r*c4*f1*u0*u4
     &  - 4.D0*qq*svp*dsvm*F10*F13r*c3*f1*u0*u3
     &  - 4.D0*qq*svp*dsvm*F10*F12r*c2*f1*u0*u2
     &  - 4.D0*qq*svp*dsvm*F10*F11r*c1*f1*u0*u1
     &  - 4.D0*qq*svp*dsvm*F10*F10r*c0*f1*u0**2
     &  + 4.D0*PC2*svm*dsvp*F15*F15r*c5*f1*u5**2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F15*F14r*c4*f1*u4*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F15*F13r*c3*f1*u3*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F15*F12r*c2*f1*u2*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F15*F11r*c1*f1*u1*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F15*F10r*c0*f1*u0*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svm*dsvp*F14*F15r*c5*f1*u4*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F14*F14r*c4*f1*u4**2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F14*F13r*c3*f1*u3*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F14*F12r*c2*f1*u2*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F14*F11r*c1*f1*u1*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F14*F10r*c0*f1*u0*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F13*F15r*c5*f1*u3*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F13*F14r*c4*f1*u3*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F13*F13r*c3*f1*u3**2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F13*F12r*c2*f1*u2*u3*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F13*F11r*c1*f1*u1*u3*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F13*F10r*c0*f1*u0*u3*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F12*F15r*c5*f1*u2*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F12*F14r*c4*f1*u2*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F12*F13r*c3*f1*u2*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svm*dsvp*F12*F12r*c2*f1*u2**2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F12*F11r*c1*f1*u1*u2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F12*F10r*c0*f1*u0*u2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F11*F15r*c5*f1*u1*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F11*F14r*c4*f1*u1*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F11*F13r*c3*f1*u1*u3*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F11*F12r*c2*f1*u1*u2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F11*F11r*c1*f1*u1**2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F11*F10r*c0*f1*u0*u1*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F10*F15r*c5*f1*u0*u5*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F10*F14r*c4*f1*u0*u4*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F10*F13r*c3*f1*u0*u3*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F10*F12r*c2*f1*u0*u2*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F10*F11r*c1*f1*u0*u1*w1*w2
     &  + 4.D0*PC2*svm*dsvp*F10*F10r*c0*f1*u0**2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2*svm*dssp*F25*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC2*svm*dssp*F25*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC2*svm*dssp*F25*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC2*svm*dssp*F25*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC2*svm*dssp*F25*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC2*svm*dssp*F25*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC2*svm*dssp*F24*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC2*svm*dssp*F24*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC2*svm*dssp*F24*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC2*svm*dssp*F24*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC2*svm*dssp*F24*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC2*svm*dssp*F24*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC2*svm*dssp*F23*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC2*svm*dssp*F23*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC2*svm*dssp*F23*F13r*c3*f1*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC2*svm*dssp*F23*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC2*svm*dssp*F23*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC2*svm*dssp*F23*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC2*svm*dssp*F22*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC2*svm*dssp*F22*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC2*svm*dssp*F22*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC2*svm*dssp*F22*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC2*svm*dssp*F22*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC2*svm*dssp*F22*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC2*svm*dssp*F21*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC2*svm*dssp*F21*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC2*svm*dssp*F21*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC2*svm*dssp*F21*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC2*svm*dssp*F21*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC2*svm*dssp*F21*F10r*c0*f1*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC2*svm*dssp*F20*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC2*svm*dssp*F20*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC2*svm*dssp*F20*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC2*svm*dssp*F20*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC2*svm*dssp*F20*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC2*svm*dssp*F20*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC2*svm*dssp*F15*F25r*c5*f2*u5**2*w2
     &  - 4.D0*PC2*svm*dssp*F15*F24r*c4*f2*u4*u5*w2
     &  - 4.D0*PC2*svm*dssp*F15*F23r*c3*f2*u3*u5*w2
     &  - 4.D0*PC2*svm*dssp*F15*F22r*c2*f2*u2*u5*w2
     &  - 4.D0*PC2*svm*dssp*F15*F21r*c1*f2*u1*u5*w2
     &  - 4.D0*PC2*svm*dssp*F15*F20r*c0*f2*u0*u5*w2
     &  - 4.D0*PC2*svm*dssp*F14*F25r*c5*f2*u4*u5*w2
     &  - 4.D0*PC2*svm*dssp*F14*F24r*c4*f2*u4**2*w2
     &  - 4.D0*PC2*svm*dssp*F14*F23r*c3*f2*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*PC2*svm*dssp*F14*F22r*c2*f2*u2*u4*w2
     &  - 4.D0*PC2*svm*dssp*F14*F21r*c1*f2*u1*u4*w2
     &  - 4.D0*PC2*svm*dssp*F14*F20r*c0*f2*u0*u4*w2
     &  - 4.D0*PC2*svm*dssp*F13*F25r*c5*f2*u3*u5*w2
     &  - 4.D0*PC2*svm*dssp*F13*F24r*c4*f2*u3*u4*w2
     &  - 4.D0*PC2*svm*dssp*F13*F23r*c3*f2*u3**2*w2
     &  - 4.D0*PC2*svm*dssp*F13*F22r*c2*f2*u2*u3*w2
     &  - 4.D0*PC2*svm*dssp*F13*F21r*c1*f2*u1*u3*w2
     &  - 4.D0*PC2*svm*dssp*F13*F20r*c0*f2*u0*u3*w2
     &  - 4.D0*PC2*svm*dssp*F12*F25r*c5*f2*u2*u5*w2
     &  - 4.D0*PC2*svm*dssp*F12*F24r*c4*f2*u2*u4*w2
     &  - 4.D0*PC2*svm*dssp*F12*F23r*c3*f2*u2*u3*w2
     &  - 4.D0*PC2*svm*dssp*F12*F22r*c2*f2*u2**2*w2
     &  - 4.D0*PC2*svm*dssp*F12*F21r*c1*f2*u1*u2*w2
     &  - 4.D0*PC2*svm*dssp*F12*F20r*c0*f2*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*PC2*svm*dssp*F11*F25r*c5*f2*u1*u5*w2
     &  - 4.D0*PC2*svm*dssp*F11*F24r*c4*f2*u1*u4*w2
     &  - 4.D0*PC2*svm*dssp*F11*F23r*c3*f2*u1*u3*w2
     &  - 4.D0*PC2*svm*dssp*F11*F22r*c2*f2*u1*u2*w2
     &  - 4.D0*PC2*svm*dssp*F11*F21r*c1*f2*u1**2*w2
     &  - 4.D0*PC2*svm*dssp*F11*F20r*c0*f2*u0*u1*w2
     &  - 4.D0*PC2*svm*dssp*F10*F25r*c5*f2*u0*u5*w2
     &  - 4.D0*PC2*svm*dssp*F10*F24r*c4*f2*u0*u4*w2
     &  - 4.D0*PC2*svm*dssp*F10*F23r*c3*f2*u0*u3*w2
     &  - 4.D0*PC2*svm*dssp*F10*F22r*c2*f2*u0*u2*w2
     &  - 4.D0*PC2*svm*dssp*F10*F21r*c1*f2*u0*u1*w2
     &  - 4.D0*PC2*svm*dssp*F10*F20r*c0*f2*u0**2*w2
     &  + 4.D0*PC2*svp*dsvm*F15*F15r*c5*f1*u5**2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F15*F14r*c4*f1*u4*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F15*F13r*c3*f1*u3*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svp*dsvm*F15*F12r*c2*f1*u2*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F15*F11r*c1*f1*u1*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F15*F10r*c0*f1*u0*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F14*F15r*c5*f1*u4*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F14*F14r*c4*f1*u4**2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F14*F13r*c3*f1*u3*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F14*F12r*c2*f1*u2*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F14*F11r*c1*f1*u1*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F14*F10r*c0*f1*u0*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F13*F15r*c5*f1*u3*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F13*F14r*c4*f1*u3*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F13*F13r*c3*f1*u3**2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F13*F12r*c2*f1*u2*u3*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F13*F11r*c1*f1*u1*u3*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F13*F10r*c0*f1*u0*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svp*dsvm*F12*F15r*c5*f1*u2*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F12*F14r*c4*f1*u2*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F12*F13r*c3*f1*u2*u3*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F12*F12r*c2*f1*u2**2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F12*F11r*c1*f1*u1*u2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F12*F10r*c0*f1*u0*u2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F11*F15r*c5*f1*u1*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F11*F14r*c4*f1*u1*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F11*F13r*c3*f1*u1*u3*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F11*F12r*c2*f1*u1*u2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F11*F11r*c1*f1*u1**2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F11*F10r*c0*f1*u0*u1*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F10*F15r*c5*f1*u0*u5*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F10*F14r*c4*f1*u0*u4*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F10*F13r*c3*f1*u0*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC2*svp*dsvm*F10*F12r*c2*f1*u0*u2*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F10*F11r*c1*f1*u0*u1*w1*w2
     &  + 4.D0*PC2*svp*dsvm*F10*F10r*c0*f1*u0**2*w1*w2
     &  - 4.D0*PC2*svp*dssm*F25*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC2*svp*dssm*F25*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC2*svp*dssm*F25*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC2*svp*dssm*F25*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC2*svp*dssm*F25*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC2*svp*dssm*F25*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC2*svp*dssm*F24*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC2*svp*dssm*F24*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC2*svp*dssm*F24*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC2*svp*dssm*F24*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC2*svp*dssm*F24*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC2*svp*dssm*F24*F10r*c0*f1*u0*u4*w1
     &
      traza1 = traza1 - 4.D0*PC2*svp*dssm*F23*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC2*svp*dssm*F23*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC2*svp*dssm*F23*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC2*svp*dssm*F23*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC2*svp*dssm*F23*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC2*svp*dssm*F23*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC2*svp*dssm*F22*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC2*svp*dssm*F22*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC2*svp*dssm*F22*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC2*svp*dssm*F22*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC2*svp*dssm*F22*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC2*svp*dssm*F22*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC2*svp*dssm*F21*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC2*svp*dssm*F21*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC2*svp*dssm*F21*F13r*c3*f1*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC2*svp*dssm*F21*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC2*svp*dssm*F21*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC2*svp*dssm*F21*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC2*svp*dssm*F20*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC2*svp*dssm*F20*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC2*svp*dssm*F20*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC2*svp*dssm*F20*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC2*svp*dssm*F20*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC2*svp*dssm*F20*F10r*c0*f1*u0**2*w1
     &  - 4.D0*PC2*svp*dssm*F15*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC2*svp*dssm*F15*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC2*svp*dssm*F15*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC2*svp*dssm*F15*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC2*svp*dssm*F15*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC2*svp*dssm*F15*F20r*c0*f2*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC2*svp*dssm*F14*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*PC2*svp*dssm*F14*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC2*svp*dssm*F14*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC2*svp*dssm*F14*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC2*svp*dssm*F14*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC2*svp*dssm*F14*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC2*svp*dssm*F13*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC2*svp*dssm*F13*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC2*svp*dssm*F13*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC2*svp*dssm*F13*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC2*svp*dssm*F13*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC2*svp*dssm*F13*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC2*svp*dssm*F12*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*PC2*svp*dssm*F12*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC2*svp*dssm*F12*F23r*c3*f2*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC2*svp*dssm*F12*F22r*c2*f2*u2**2*w1
     &  - 4.D0*PC2*svp*dssm*F12*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC2*svp*dssm*F12*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC2*svp*dssm*F11*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC2*svp*dssm*F11*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC2*svp*dssm*F11*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC2*svp*dssm*F11*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC2*svp*dssm*F11*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC2*svp*dssm*F11*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC2*svp*dssm*F10*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC2*svp*dssm*F10*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC2*svp*dssm*F10*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC2*svp*dssm*F10*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC2*svp*dssm*F10*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC2*svp*dssm*F10*F20r*c0*f2*u0**2*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*dsvp*F25*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F25*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F25*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F25*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F25*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F25*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F24*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F24*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F24*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F24*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F24*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F24*F10r*c0*f1*u0*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F23*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F23*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F23*F13r*c3*f1*u3**2*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*dsvp*F23*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F23*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F23*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F22*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F22*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F22*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F22*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F22*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F22*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F21*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F21*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F21*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F21*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F21*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F21*F10r*c0*f1*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*dsvp*F20*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F20*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F20*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F20*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F20*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC2*ssm*dsvp*F20*F10r*c0*f1*u0**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F15*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F15*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F15*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F15*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F15*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F15*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F14*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F14*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F14*F23r*c3*f2*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*dsvp*F14*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F14*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F14*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F13*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F13*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F13*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F13*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F13*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F13*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F12*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F12*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F12*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F12*F22r*c2*f2*u2**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F12*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F12*F20r*c0*f2*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC2*ssm*dsvp*F11*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F11*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F11*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F11*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F11*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC2*ssm*dsvp*F11*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC2*ssm*dsvp*F10*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC2*ssm*dsvp*F10*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC2*ssm*dsvp*F10*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC2*ssm*dsvp*F10*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC2*ssm*dsvp*F10*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC2*ssm*dsvp*F10*F20r*c0*f2*u0**2*w1
     &  + 4.D0*PC2*ssm*dssp*F25*F25r*c5*f2*u5**2
     &  + 4.D0*PC2*ssm*dssp*F25*F24r*c4*f2*u4*u5
     &  + 4.D0*PC2*ssm*dssp*F25*F23r*c3*f2*u3*u5
     &
      traza1 = traza1 + 4.D0*PC2*ssm*dssp*F25*F22r*c2*f2*u2*u5
     &  + 4.D0*PC2*ssm*dssp*F25*F21r*c1*f2*u1*u5
     &  + 4.D0*PC2*ssm*dssp*F25*F20r*c0*f2*u0*u5
     &  + 4.D0*PC2*ssm*dssp*F24*F25r*c5*f2*u4*u5
     &  + 4.D0*PC2*ssm*dssp*F24*F24r*c4*f2*u4**2
     &  + 4.D0*PC2*ssm*dssp*F24*F23r*c3*f2*u3*u4
     &  + 4.D0*PC2*ssm*dssp*F24*F22r*c2*f2*u2*u4
     &  + 4.D0*PC2*ssm*dssp*F24*F21r*c1*f2*u1*u4
     &  + 4.D0*PC2*ssm*dssp*F24*F20r*c0*f2*u0*u4
     &  + 4.D0*PC2*ssm*dssp*F23*F25r*c5*f2*u3*u5
     &  + 4.D0*PC2*ssm*dssp*F23*F24r*c4*f2*u3*u4
     &  + 4.D0*PC2*ssm*dssp*F23*F23r*c3*f2*u3**2
     &  + 4.D0*PC2*ssm*dssp*F23*F22r*c2*f2*u2*u3
     &  + 4.D0*PC2*ssm*dssp*F23*F21r*c1*f2*u1*u3
     &  + 4.D0*PC2*ssm*dssp*F23*F20r*c0*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC2*ssm*dssp*F22*F25r*c5*f2*u2*u5
     &  + 4.D0*PC2*ssm*dssp*F22*F24r*c4*f2*u2*u4
     &  + 4.D0*PC2*ssm*dssp*F22*F23r*c3*f2*u2*u3
     &  + 4.D0*PC2*ssm*dssp*F22*F22r*c2*f2*u2**2
     &  + 4.D0*PC2*ssm*dssp*F22*F21r*c1*f2*u1*u2
     &  + 4.D0*PC2*ssm*dssp*F22*F20r*c0*f2*u0*u2
     &  + 4.D0*PC2*ssm*dssp*F21*F25r*c5*f2*u1*u5
     &  + 4.D0*PC2*ssm*dssp*F21*F24r*c4*f2*u1*u4
     &  + 4.D0*PC2*ssm*dssp*F21*F23r*c3*f2*u1*u3
     &  + 4.D0*PC2*ssm*dssp*F21*F22r*c2*f2*u1*u2
     &  + 4.D0*PC2*ssm*dssp*F21*F21r*c1*f2*u1**2
     &  + 4.D0*PC2*ssm*dssp*F21*F20r*c0*f2*u0*u1
     &  + 4.D0*PC2*ssm*dssp*F20*F25r*c5*f2*u0*u5
     &  + 4.D0*PC2*ssm*dssp*F20*F24r*c4*f2*u0*u4
     &  + 4.D0*PC2*ssm*dssp*F20*F23r*c3*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC2*ssm*dssp*F20*F22r*c2*f2*u0*u2
     &  + 4.D0*PC2*ssm*dssp*F20*F21r*c1*f2*u0*u1
     &  + 4.D0*PC2*ssm*dssp*F20*F20r*c0*f2*u0**2
     &  - 4.D0*PC2*ssp*dsvm*F25*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F25*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F25*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F25*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F25*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F25*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F24*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F24*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F24*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F24*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F24*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F24*F10r*c0*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*dsvm*F23*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F23*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F23*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F23*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F23*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F23*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F22*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F22*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F22*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F22*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F22*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F22*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F21*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F21*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F21*F13r*c3*f1*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*dsvm*F21*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F21*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F21*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC2*ssp*dsvm*F20*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F20*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F20*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F20*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F20*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC2*ssp*dsvm*F20*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F15*F25r*c5*f2*u5**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F15*F24r*c4*f2*u4*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F15*F23r*c3*f2*u3*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F15*F22r*c2*f2*u2*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F15*F21r*c1*f2*u1*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F15*F20r*c0*f2*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*dsvm*F14*F25r*c5*f2*u4*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F14*F24r*c4*f2*u4**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F14*F23r*c3*f2*u3*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F14*F22r*c2*f2*u2*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F14*F21r*c1*f2*u1*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F14*F20r*c0*f2*u0*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F13*F25r*c5*f2*u3*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F13*F24r*c4*f2*u3*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F13*F23r*c3*f2*u3**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F13*F22r*c2*f2*u2*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F13*F21r*c1*f2*u1*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F13*F20r*c0*f2*u0*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F12*F25r*c5*f2*u2*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F12*F24r*c4*f2*u2*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F12*F23r*c3*f2*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC2*ssp*dsvm*F12*F22r*c2*f2*u2**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F12*F21r*c1*f2*u1*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F12*F20r*c0*f2*u0*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F11*F25r*c5*f2*u1*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F11*F24r*c4*f2*u1*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F11*F23r*c3*f2*u1*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F11*F22r*c2*f2*u1*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F11*F21r*c1*f2*u1**2*w2
     &  - 4.D0*PC2*ssp*dsvm*F11*F20r*c0*f2*u0*u1*w2
     &  - 4.D0*PC2*ssp*dsvm*F10*F25r*c5*f2*u0*u5*w2
     &  - 4.D0*PC2*ssp*dsvm*F10*F24r*c4*f2*u0*u4*w2
     &  - 4.D0*PC2*ssp*dsvm*F10*F23r*c3*f2*u0*u3*w2
     &  - 4.D0*PC2*ssp*dsvm*F10*F22r*c2*f2*u0*u2*w2
     &  - 4.D0*PC2*ssp*dsvm*F10*F21r*c1*f2*u0*u1*w2
     &  - 4.D0*PC2*ssp*dsvm*F10*F20r*c0*f2*u0**2*w2
     &
      traza1 = traza1 + 4.D0*PC2*ssp*dssm*F25*F25r*c5*f2*u5**2
     &  + 4.D0*PC2*ssp*dssm*F25*F24r*c4*f2*u4*u5
     &  + 4.D0*PC2*ssp*dssm*F25*F23r*c3*f2*u3*u5
     &  + 4.D0*PC2*ssp*dssm*F25*F22r*c2*f2*u2*u5
     &  + 4.D0*PC2*ssp*dssm*F25*F21r*c1*f2*u1*u5
     &  + 4.D0*PC2*ssp*dssm*F25*F20r*c0*f2*u0*u5
     &  + 4.D0*PC2*ssp*dssm*F24*F25r*c5*f2*u4*u5
     &  + 4.D0*PC2*ssp*dssm*F24*F24r*c4*f2*u4**2
     &  + 4.D0*PC2*ssp*dssm*F24*F23r*c3*f2*u3*u4
     &  + 4.D0*PC2*ssp*dssm*F24*F22r*c2*f2*u2*u4
     &  + 4.D0*PC2*ssp*dssm*F24*F21r*c1*f2*u1*u4
     &  + 4.D0*PC2*ssp*dssm*F24*F20r*c0*f2*u0*u4
     &  + 4.D0*PC2*ssp*dssm*F23*F25r*c5*f2*u3*u5
     &  + 4.D0*PC2*ssp*dssm*F23*F24r*c4*f2*u3*u4
     &  + 4.D0*PC2*ssp*dssm*F23*F23r*c3*f2*u3**2
     &
      traza1 = traza1 + 4.D0*PC2*ssp*dssm*F23*F22r*c2*f2*u2*u3
     &  + 4.D0*PC2*ssp*dssm*F23*F21r*c1*f2*u1*u3
     &  + 4.D0*PC2*ssp*dssm*F23*F20r*c0*f2*u0*u3
     &  + 4.D0*PC2*ssp*dssm*F22*F25r*c5*f2*u2*u5
     &  + 4.D0*PC2*ssp*dssm*F22*F24r*c4*f2*u2*u4
     &  + 4.D0*PC2*ssp*dssm*F22*F23r*c3*f2*u2*u3
     &  + 4.D0*PC2*ssp*dssm*F22*F22r*c2*f2*u2**2
     &  + 4.D0*PC2*ssp*dssm*F22*F21r*c1*f2*u1*u2
     &  + 4.D0*PC2*ssp*dssm*F22*F20r*c0*f2*u0*u2
     &  + 4.D0*PC2*ssp*dssm*F21*F25r*c5*f2*u1*u5
     &  + 4.D0*PC2*ssp*dssm*F21*F24r*c4*f2*u1*u4
     &  + 4.D0*PC2*ssp*dssm*F21*F23r*c3*f2*u1*u3
     &  + 4.D0*PC2*ssp*dssm*F21*F22r*c2*f2*u1*u2
     &  + 4.D0*PC2*ssp*dssm*F21*F21r*c1*f2*u1**2
     &  + 4.D0*PC2*ssp*dssm*F21*F20r*c0*f2*u0*u1
     &
      traza1 = traza1 + 4.D0*PC2*ssp*dssm*F20*F25r*c5*f2*u0*u5
     &  + 4.D0*PC2*ssp*dssm*F20*F24r*c4*f2*u0*u4
     &  + 4.D0*PC2*ssp*dssm*F20*F23r*c3*f2*u0*u3
     &  + 4.D0*PC2*ssp*dssm*F20*F22r*c2*f2*u0*u2
     &  + 4.D0*PC2*ssp*dssm*F20*F21r*c1*f2*u0*u1
     &  + 4.D0*PC2*ssp*dssm*F20*F20r*c0*f2*u0**2
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F15r*c5*f1*u5**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F15r*c5*f1*u5**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F14r*c4*f1*u4*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F14r*c4*f1*u4*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F13r*c3*f1*u3*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F13r*c3*f1*u3*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F12r*c2*f1*u2*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F12r*c2*f1*u2*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F11r*c1*f1*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F45*F11r*c1*f1*u1*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F10r*c0*f1*u0*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F45*F10r*c0*f1*u0*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F15r*c5*f1*u4*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F15r*c5*f1*u4*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F14r*c4*f1*u4**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F14r*c4*f1*u4**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F13r*c3*f1*u3*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F13r*c3*f1*u3*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F12r*c2*f1*u2*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F12r*c2*f1*u2*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F11r*c1*f1*u1*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F11r*c1*f1*u1*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F10r*c0*f1*u0*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F44*F10r*c0*f1*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F43*F15r*c5*f1*u3*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F15r*c5*f1*u3*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F14r*c4*f1*u3*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F14r*c4*f1*u3*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F13r*c3*f1*u3**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F13r*c3*f1*u3**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F12r*c2*f1*u2*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F12r*c2*f1*u2*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F11r*c1*f1*u1*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F11r*c1*f1*u1*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F10r*c0*f1*u0*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F43*F10r*c0*f1*u0*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F15r*c5*f1*u2*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F15r*c5*f1*u2*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F14r*c4*f1*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F42*F14r*c4*f1*u2*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F13r*c3*f1*u2*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F13r*c3*f1*u2*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F12r*c2*f1*u2**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F12r*c2*f1*u2**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F11r*c1*f1*u1*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F11r*c1*f1*u1*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F10r*c0*f1*u0*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F42*F10r*c0*f1*u0*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F15r*c5*f1*u1*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F15r*c5*f1*u1*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F14r*c4*f1*u1*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F14r*c4*f1*u1*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F13r*c3*f1*u1*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F13r*c3*f1*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F41*F12r*c2*f1*u1*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F12r*c2*f1*u1*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F11r*c1*f1*u1**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F11r*c1*f1*u1**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F10r*c0*f1*u0*u1*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F41*F10r*c0*f1*u0*u1*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F15r*c5*f1*u0*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F15r*c5*f1*u0*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F14r*c4*f1*u0*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F14r*c4*f1*u0*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F13r*c3*f1*u0*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F13r*c3*f1*u0*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F12r*c2*f1*u0*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F12r*c2*f1*u0*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F11r*c1*f1*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F40*F11r*c1*f1*u0*u1*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F10r*c0*f1*u0**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F40*F10r*c0*f1*u0**2*w1
     &  - 4.D0*PC2*qq*svm*dsvp*F25*F25r*c5*f2*u5**2
     &  - 4.D0*PC2*qq*svm*dsvp*F25*F24r*c4*f2*u4*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F25*F23r*c3*f2*u3*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F25*F22r*c2*f2*u2*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F25*F21r*c1*f2*u1*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F25*F20r*c0*f2*u0*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F24*F25r*c5*f2*u4*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F24*F24r*c4*f2*u4**2
     &  - 4.D0*PC2*qq*svm*dsvp*F24*F23r*c3*f2*u3*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F24*F22r*c2*f2*u2*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F24*F21r*c1*f2*u1*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F24*F20r*c0*f2*u0*u4
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dsvp*F23*F25r*c5*f2*u3*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F23*F24r*c4*f2*u3*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F23*F23r*c3*f2*u3**2
     &  - 4.D0*PC2*qq*svm*dsvp*F23*F22r*c2*f2*u2*u3
     &  - 4.D0*PC2*qq*svm*dsvp*F23*F21r*c1*f2*u1*u3
     &  - 4.D0*PC2*qq*svm*dsvp*F23*F20r*c0*f2*u0*u3
     &  - 4.D0*PC2*qq*svm*dsvp*F22*F25r*c5*f2*u2*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F22*F24r*c4*f2*u2*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F22*F23r*c3*f2*u2*u3
     &  - 4.D0*PC2*qq*svm*dsvp*F22*F22r*c2*f2*u2**2
     &  - 4.D0*PC2*qq*svm*dsvp*F22*F21r*c1*f2*u1*u2
     &  - 4.D0*PC2*qq*svm*dsvp*F22*F20r*c0*f2*u0*u2
     &  - 4.D0*PC2*qq*svm*dsvp*F21*F25r*c5*f2*u1*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F21*F24r*c4*f2*u1*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F21*F23r*c3*f2*u1*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dsvp*F21*F22r*c2*f2*u1*u2
     &  - 4.D0*PC2*qq*svm*dsvp*F21*F21r*c1*f2*u1**2
     &  - 4.D0*PC2*qq*svm*dsvp*F21*F20r*c0*f2*u0*u1
     &  - 4.D0*PC2*qq*svm*dsvp*F20*F25r*c5*f2*u0*u5
     &  - 4.D0*PC2*qq*svm*dsvp*F20*F24r*c4*f2*u0*u4
     &  - 4.D0*PC2*qq*svm*dsvp*F20*F23r*c3*f2*u0*u3
     &  - 4.D0*PC2*qq*svm*dsvp*F20*F22r*c2*f2*u0*u2
     &  - 4.D0*PC2*qq*svm*dsvp*F20*F21r*c1*f2*u0*u1
     &  - 4.D0*PC2*qq*svm*dsvp*F20*F20r*c0*f2*u0**2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F43r*c3*f4*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F15*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F15*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F41r*c1*f4*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F14*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F14*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F13*F40r*c0*f4*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F12*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F12*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F44r*c4*f4*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F11*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F11*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F43r*c3*f4*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svm*dsvp*F10*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC2*qq*svm*dsvp*F10*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC2*qq*svm*dssp*F45*F25r*c5*f2*u5**2
     &  - 4.D0*PC2*qq*svm*dssp*F45*F24r*c4*f2*u4*u5
     &  - 4.D0*PC2*qq*svm*dssp*F45*F23r*c3*f2*u3*u5
     &  - 4.D0*PC2*qq*svm*dssp*F45*F22r*c2*f2*u2*u5
     &  - 4.D0*PC2*qq*svm*dssp*F45*F21r*c1*f2*u1*u5
     &  - 4.D0*PC2*qq*svm*dssp*F45*F20r*c0*f2*u0*u5
     &  - 4.D0*PC2*qq*svm*dssp*F44*F25r*c5*f2*u4*u5
     &  - 4.D0*PC2*qq*svm*dssp*F44*F24r*c4*f2*u4**2
     &  - 4.D0*PC2*qq*svm*dssp*F44*F23r*c3*f2*u3*u4
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dssp*F44*F22r*c2*f2*u2*u4
     &  - 4.D0*PC2*qq*svm*dssp*F44*F21r*c1*f2*u1*u4
     &  - 4.D0*PC2*qq*svm*dssp*F44*F20r*c0*f2*u0*u4
     &  - 4.D0*PC2*qq*svm*dssp*F43*F25r*c5*f2*u3*u5
     &  - 4.D0*PC2*qq*svm*dssp*F43*F24r*c4*f2*u3*u4
     &  - 4.D0*PC2*qq*svm*dssp*F43*F23r*c3*f2*u3**2
     &  - 4.D0*PC2*qq*svm*dssp*F43*F22r*c2*f2*u2*u3
     &  - 4.D0*PC2*qq*svm*dssp*F43*F21r*c1*f2*u1*u3
     &  - 4.D0*PC2*qq*svm*dssp*F43*F20r*c0*f2*u0*u3
     &  - 4.D0*PC2*qq*svm*dssp*F42*F25r*c5*f2*u2*u5
     &  - 4.D0*PC2*qq*svm*dssp*F42*F24r*c4*f2*u2*u4
     &  - 4.D0*PC2*qq*svm*dssp*F42*F23r*c3*f2*u2*u3
     &  - 4.D0*PC2*qq*svm*dssp*F42*F22r*c2*f2*u2**2
     &  - 4.D0*PC2*qq*svm*dssp*F42*F21r*c1*f2*u1*u2
     &  - 4.D0*PC2*qq*svm*dssp*F42*F20r*c0*f2*u0*u2
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dssp*F41*F25r*c5*f2*u1*u5
     &  - 4.D0*PC2*qq*svm*dssp*F41*F24r*c4*f2*u1*u4
     &  - 4.D0*PC2*qq*svm*dssp*F41*F23r*c3*f2*u1*u3
     &  - 4.D0*PC2*qq*svm*dssp*F41*F22r*c2*f2*u1*u2
     &  - 4.D0*PC2*qq*svm*dssp*F41*F21r*c1*f2*u1**2
     &  - 4.D0*PC2*qq*svm*dssp*F41*F20r*c0*f2*u0*u1
     &  - 4.D0*PC2*qq*svm*dssp*F40*F25r*c5*f2*u0*u5
     &  - 4.D0*PC2*qq*svm*dssp*F40*F24r*c4*f2*u0*u4
     &  - 4.D0*PC2*qq*svm*dssp*F40*F23r*c3*f2*u0*u3
     &  - 4.D0*PC2*qq*svm*dssp*F40*F22r*c2*f2*u0*u2
     &  - 4.D0*PC2*qq*svm*dssp*F40*F21r*c1*f2*u0*u1
     &  - 4.D0*PC2*qq*svm*dssp*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC2*qq*svm*dssp*F25*F45r*c5*f4*u5**2
     &  - 4.D0*PC2*qq*svm*dssp*F25*F44r*c4*f4*u4*u5
     &  - 4.D0*PC2*qq*svm*dssp*F25*F43r*c3*f4*u3*u5
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dssp*F25*F42r*c2*f4*u2*u5
     &  - 4.D0*PC2*qq*svm*dssp*F25*F41r*c1*f4*u1*u5
     &  - 4.D0*PC2*qq*svm*dssp*F25*F40r*c0*f4*u0*u5
     &  - 4.D0*PC2*qq*svm*dssp*F24*F45r*c5*f4*u4*u5
     &  - 4.D0*PC2*qq*svm*dssp*F24*F44r*c4*f4*u4**2
     &  - 4.D0*PC2*qq*svm*dssp*F24*F43r*c3*f4*u3*u4
     &  - 4.D0*PC2*qq*svm*dssp*F24*F42r*c2*f4*u2*u4
     &  - 4.D0*PC2*qq*svm*dssp*F24*F41r*c1*f4*u1*u4
     &  - 4.D0*PC2*qq*svm*dssp*F24*F40r*c0*f4*u0*u4
     &  - 4.D0*PC2*qq*svm*dssp*F23*F45r*c5*f4*u3*u5
     &  - 4.D0*PC2*qq*svm*dssp*F23*F44r*c4*f4*u3*u4
     &  - 4.D0*PC2*qq*svm*dssp*F23*F43r*c3*f4*u3**2
     &  - 4.D0*PC2*qq*svm*dssp*F23*F42r*c2*f4*u2*u3
     &  - 4.D0*PC2*qq*svm*dssp*F23*F41r*c1*f4*u1*u3
     &  - 4.D0*PC2*qq*svm*dssp*F23*F40r*c0*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dssp*F22*F45r*c5*f4*u2*u5
     &  - 4.D0*PC2*qq*svm*dssp*F22*F44r*c4*f4*u2*u4
     &  - 4.D0*PC2*qq*svm*dssp*F22*F43r*c3*f4*u2*u3
     &  - 4.D0*PC2*qq*svm*dssp*F22*F42r*c2*f4*u2**2
     &  - 4.D0*PC2*qq*svm*dssp*F22*F41r*c1*f4*u1*u2
     &  - 4.D0*PC2*qq*svm*dssp*F22*F40r*c0*f4*u0*u2
     &  - 4.D0*PC2*qq*svm*dssp*F21*F45r*c5*f4*u1*u5
     &  - 4.D0*PC2*qq*svm*dssp*F21*F44r*c4*f4*u1*u4
     &  - 4.D0*PC2*qq*svm*dssp*F21*F43r*c3*f4*u1*u3
     &  - 4.D0*PC2*qq*svm*dssp*F21*F42r*c2*f4*u1*u2
     &  - 4.D0*PC2*qq*svm*dssp*F21*F41r*c1*f4*u1**2
     &  - 4.D0*PC2*qq*svm*dssp*F21*F40r*c0*f4*u0*u1
     &  - 4.D0*PC2*qq*svm*dssp*F20*F45r*c5*f4*u0*u5
     &  - 4.D0*PC2*qq*svm*dssp*F20*F44r*c4*f4*u0*u4
     &  - 4.D0*PC2*qq*svm*dssp*F20*F43r*c3*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*svm*dssp*F20*F42r*c2*f4*u0*u2
     &  - 4.D0*PC2*qq*svm*dssp*F20*F41r*c1*f4*u0*u1
     &  - 4.D0*PC2*qq*svm*dssp*F20*F40r*c0*f4*u0**2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F15r*c5*f1*u5**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F15r*c5*f1*u5**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F14r*c4*f1*u4*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F14r*c4*f1*u4*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F13r*c3*f1*u3*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F13r*c3*f1*u3*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F12r*c2*f1*u2*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F12r*c2*f1*u2*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F11r*c1*f1*u1*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F11r*c1*f1*u1*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F10r*c0*f1*u0*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F45*F10r*c0*f1*u0*u5*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F44*F15r*c5*f1*u4*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F15r*c5*f1*u4*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F14r*c4*f1*u4**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F14r*c4*f1*u4**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F13r*c3*f1*u3*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F13r*c3*f1*u3*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F12r*c2*f1*u2*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F12r*c2*f1*u2*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F11r*c1*f1*u1*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F11r*c1*f1*u1*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F10r*c0*f1*u0*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F44*F10r*c0*f1*u0*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F15r*c5*f1*u3*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F15r*c5*f1*u3*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F14r*c4*f1*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F43*F14r*c4*f1*u3*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F13r*c3*f1*u3**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F13r*c3*f1*u3**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F12r*c2*f1*u2*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F12r*c2*f1*u2*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F11r*c1*f1*u1*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F11r*c1*f1*u1*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F10r*c0*f1*u0*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F43*F10r*c0*f1*u0*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F15r*c5*f1*u2*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F15r*c5*f1*u2*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F14r*c4*f1*u2*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F14r*c4*f1*u2*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F13r*c3*f1*u2*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F13r*c3*f1*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F42*F12r*c2*f1*u2**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F12r*c2*f1*u2**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F11r*c1*f1*u1*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F11r*c1*f1*u1*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F10r*c0*f1*u0*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F42*F10r*c0*f1*u0*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F15r*c5*f1*u1*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F15r*c5*f1*u1*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F14r*c4*f1*u1*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F14r*c4*f1*u1*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F13r*c3*f1*u1*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F13r*c3*f1*u1*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F12r*c2*f1*u1*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F12r*c2*f1*u1*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F11r*c1*f1*u1**2*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F41*F11r*c1*f1*u1**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F10r*c0*f1*u0*u1*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F41*F10r*c0*f1*u0*u1*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F15r*c5*f1*u0*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F15r*c5*f1*u0*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F14r*c4*f1*u0*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F14r*c4*f1*u0*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F13r*c3*f1*u0*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F13r*c3*f1*u0*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F12r*c2*f1*u0*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F12r*c2*f1*u0*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F11r*c1*f1*u0*u1*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F11r*c1*f1*u0*u1*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F10r*c0*f1*u0**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F40*F10r*c0*f1*u0**2*w1
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dsvm*F25*F25r*c5*f2*u5**2
     &  - 4.D0*PC2*qq*svp*dsvm*F25*F24r*c4*f2*u4*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F25*F23r*c3*f2*u3*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F25*F22r*c2*f2*u2*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F25*F21r*c1*f2*u1*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F25*F20r*c0*f2*u0*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F24*F25r*c5*f2*u4*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F24*F24r*c4*f2*u4**2
     &  - 4.D0*PC2*qq*svp*dsvm*F24*F23r*c3*f2*u3*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F24*F22r*c2*f2*u2*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F24*F21r*c1*f2*u1*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F24*F20r*c0*f2*u0*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F23*F25r*c5*f2*u3*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F23*F24r*c4*f2*u3*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F23*F23r*c3*f2*u3**2
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dsvm*F23*F22r*c2*f2*u2*u3
     &  - 4.D0*PC2*qq*svp*dsvm*F23*F21r*c1*f2*u1*u3
     &  - 4.D0*PC2*qq*svp*dsvm*F23*F20r*c0*f2*u0*u3
     &  - 4.D0*PC2*qq*svp*dsvm*F22*F25r*c5*f2*u2*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F22*F24r*c4*f2*u2*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F22*F23r*c3*f2*u2*u3
     &  - 4.D0*PC2*qq*svp*dsvm*F22*F22r*c2*f2*u2**2
     &  - 4.D0*PC2*qq*svp*dsvm*F22*F21r*c1*f2*u1*u2
     &  - 4.D0*PC2*qq*svp*dsvm*F22*F20r*c0*f2*u0*u2
     &  - 4.D0*PC2*qq*svp*dsvm*F21*F25r*c5*f2*u1*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F21*F24r*c4*f2*u1*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F21*F23r*c3*f2*u1*u3
     &  - 4.D0*PC2*qq*svp*dsvm*F21*F22r*c2*f2*u1*u2
     &  - 4.D0*PC2*qq*svp*dsvm*F21*F21r*c1*f2*u1**2
     &  - 4.D0*PC2*qq*svp*dsvm*F21*F20r*c0*f2*u0*u1
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dsvm*F20*F25r*c5*f2*u0*u5
     &  - 4.D0*PC2*qq*svp*dsvm*F20*F24r*c4*f2*u0*u4
     &  - 4.D0*PC2*qq*svp*dsvm*F20*F23r*c3*f2*u0*u3
     &  - 4.D0*PC2*qq*svp*dsvm*F20*F22r*c2*f2*u0*u2
     &  - 4.D0*PC2*qq*svp*dsvm*F20*F21r*c1*f2*u0*u1
     &  - 4.D0*PC2*qq*svp*dsvm*F20*F20r*c0*f2*u0**2
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F41r*c1*f4*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F15*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F15*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F14*F40r*c0*f4*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F13*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F13*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F44r*c4*f4*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F12*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F12*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F43r*c3*f4*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F11*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F11*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F41r*c1*f4*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC2*qq*svp*dsvm*F10*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC2*qq*svp*dsvm*F10*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC2*qq*svp*dssm*F45*F25r*c5*f2*u5**2
     &  - 4.D0*PC2*qq*svp*dssm*F45*F24r*c4*f2*u4*u5
     &  - 4.D0*PC2*qq*svp*dssm*F45*F23r*c3*f2*u3*u5
     &  - 4.D0*PC2*qq*svp*dssm*F45*F22r*c2*f2*u2*u5
     &  - 4.D0*PC2*qq*svp*dssm*F45*F21r*c1*f2*u1*u5
     &  - 4.D0*PC2*qq*svp*dssm*F45*F20r*c0*f2*u0*u5
     &  - 4.D0*PC2*qq*svp*dssm*F44*F25r*c5*f2*u4*u5
     &  - 4.D0*PC2*qq*svp*dssm*F44*F24r*c4*f2*u4**2
     &  - 4.D0*PC2*qq*svp*dssm*F44*F23r*c3*f2*u3*u4
     &  - 4.D0*PC2*qq*svp*dssm*F44*F22r*c2*f2*u2*u4
     &  - 4.D0*PC2*qq*svp*dssm*F44*F21r*c1*f2*u1*u4
     &  - 4.D0*PC2*qq*svp*dssm*F44*F20r*c0*f2*u0*u4
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dssm*F43*F25r*c5*f2*u3*u5
     &  - 4.D0*PC2*qq*svp*dssm*F43*F24r*c4*f2*u3*u4
     &  - 4.D0*PC2*qq*svp*dssm*F43*F23r*c3*f2*u3**2
     &  - 4.D0*PC2*qq*svp*dssm*F43*F22r*c2*f2*u2*u3
     &  - 4.D0*PC2*qq*svp*dssm*F43*F21r*c1*f2*u1*u3
     &  - 4.D0*PC2*qq*svp*dssm*F43*F20r*c0*f2*u0*u3
     &  - 4.D0*PC2*qq*svp*dssm*F42*F25r*c5*f2*u2*u5
     &  - 4.D0*PC2*qq*svp*dssm*F42*F24r*c4*f2*u2*u4
     &  - 4.D0*PC2*qq*svp*dssm*F42*F23r*c3*f2*u2*u3
     &  - 4.D0*PC2*qq*svp*dssm*F42*F22r*c2*f2*u2**2
     &  - 4.D0*PC2*qq*svp*dssm*F42*F21r*c1*f2*u1*u2
     &  - 4.D0*PC2*qq*svp*dssm*F42*F20r*c0*f2*u0*u2
     &  - 4.D0*PC2*qq*svp*dssm*F41*F25r*c5*f2*u1*u5
     &  - 4.D0*PC2*qq*svp*dssm*F41*F24r*c4*f2*u1*u4
     &  - 4.D0*PC2*qq*svp*dssm*F41*F23r*c3*f2*u1*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dssm*F41*F22r*c2*f2*u1*u2
     &  - 4.D0*PC2*qq*svp*dssm*F41*F21r*c1*f2*u1**2
     &  - 4.D0*PC2*qq*svp*dssm*F41*F20r*c0*f2*u0*u1
     &  - 4.D0*PC2*qq*svp*dssm*F40*F25r*c5*f2*u0*u5
     &  - 4.D0*PC2*qq*svp*dssm*F40*F24r*c4*f2*u0*u4
     &  - 4.D0*PC2*qq*svp*dssm*F40*F23r*c3*f2*u0*u3
     &  - 4.D0*PC2*qq*svp*dssm*F40*F22r*c2*f2*u0*u2
     &  - 4.D0*PC2*qq*svp*dssm*F40*F21r*c1*f2*u0*u1
     &  - 4.D0*PC2*qq*svp*dssm*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC2*qq*svp*dssm*F25*F45r*c5*f4*u5**2
     &  - 4.D0*PC2*qq*svp*dssm*F25*F44r*c4*f4*u4*u5
     &  - 4.D0*PC2*qq*svp*dssm*F25*F43r*c3*f4*u3*u5
     &  - 4.D0*PC2*qq*svp*dssm*F25*F42r*c2*f4*u2*u5
     &  - 4.D0*PC2*qq*svp*dssm*F25*F41r*c1*f4*u1*u5
     &  - 4.D0*PC2*qq*svp*dssm*F25*F40r*c0*f4*u0*u5
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dssm*F24*F45r*c5*f4*u4*u5
     &  - 4.D0*PC2*qq*svp*dssm*F24*F44r*c4*f4*u4**2
     &  - 4.D0*PC2*qq*svp*dssm*F24*F43r*c3*f4*u3*u4
     &  - 4.D0*PC2*qq*svp*dssm*F24*F42r*c2*f4*u2*u4
     &  - 4.D0*PC2*qq*svp*dssm*F24*F41r*c1*f4*u1*u4
     &  - 4.D0*PC2*qq*svp*dssm*F24*F40r*c0*f4*u0*u4
     &  - 4.D0*PC2*qq*svp*dssm*F23*F45r*c5*f4*u3*u5
     &  - 4.D0*PC2*qq*svp*dssm*F23*F44r*c4*f4*u3*u4
     &  - 4.D0*PC2*qq*svp*dssm*F23*F43r*c3*f4*u3**2
     &  - 4.D0*PC2*qq*svp*dssm*F23*F42r*c2*f4*u2*u3
     &  - 4.D0*PC2*qq*svp*dssm*F23*F41r*c1*f4*u1*u3
     &  - 4.D0*PC2*qq*svp*dssm*F23*F40r*c0*f4*u0*u3
     &  - 4.D0*PC2*qq*svp*dssm*F22*F45r*c5*f4*u2*u5
     &  - 4.D0*PC2*qq*svp*dssm*F22*F44r*c4*f4*u2*u4
     &  - 4.D0*PC2*qq*svp*dssm*F22*F43r*c3*f4*u2*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*svp*dssm*F22*F42r*c2*f4*u2**2
     &  - 4.D0*PC2*qq*svp*dssm*F22*F41r*c1*f4*u1*u2
     &  - 4.D0*PC2*qq*svp*dssm*F22*F40r*c0*f4*u0*u2
     &  - 4.D0*PC2*qq*svp*dssm*F21*F45r*c5*f4*u1*u5
     &  - 4.D0*PC2*qq*svp*dssm*F21*F44r*c4*f4*u1*u4
     &  - 4.D0*PC2*qq*svp*dssm*F21*F43r*c3*f4*u1*u3
     &  - 4.D0*PC2*qq*svp*dssm*F21*F42r*c2*f4*u1*u2
     &  - 4.D0*PC2*qq*svp*dssm*F21*F41r*c1*f4*u1**2
     &  - 4.D0*PC2*qq*svp*dssm*F21*F40r*c0*f4*u0*u1
     &  - 4.D0*PC2*qq*svp*dssm*F20*F45r*c5*f4*u0*u5
     &  - 4.D0*PC2*qq*svp*dssm*F20*F44r*c4*f4*u0*u4
     &  - 4.D0*PC2*qq*svp*dssm*F20*F43r*c3*f4*u0*u3
     &  - 4.D0*PC2*qq*svp*dssm*F20*F42r*c2*f4*u0*u2
     &  - 4.D0*PC2*qq*svp*dssm*F20*F41r*c1*f4*u0*u1
     &  - 4.D0*PC2*qq*svp*dssm*F20*F40r*c0*f4*u0**2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dsvp*F45*F25r*c5*f2*u5**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F45*F24r*c4*f2*u4*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F45*F23r*c3*f2*u3*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F45*F22r*c2*f2*u2*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F45*F21r*c1*f2*u1*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F45*F20r*c0*f2*u0*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F44*F25r*c5*f2*u4*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F44*F24r*c4*f2*u4**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F44*F23r*c3*f2*u3*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F44*F22r*c2*f2*u2*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F44*F21r*c1*f2*u1*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F44*F20r*c0*f2*u0*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F43*F25r*c5*f2*u3*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F43*F24r*c4*f2*u3*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F43*F23r*c3*f2*u3**2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dsvp*F43*F22r*c2*f2*u2*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F43*F21r*c1*f2*u1*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F43*F20r*c0*f2*u0*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F42*F25r*c5*f2*u2*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F42*F24r*c4*f2*u2*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F42*F23r*c3*f2*u2*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F42*F22r*c2*f2*u2**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F42*F21r*c1*f2*u1*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F42*F20r*c0*f2*u0*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F41*F25r*c5*f2*u1*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F41*F24r*c4*f2*u1*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F41*F23r*c3*f2*u1*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F41*F22r*c2*f2*u1*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F41*F21r*c1*f2*u1**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F41*F20r*c0*f2*u0*u1
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dsvp*F40*F25r*c5*f2*u0*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F40*F24r*c4*f2*u0*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F40*F23r*c3*f2*u0*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F40*F22r*c2*f2*u0*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F40*F21r*c1*f2*u0*u1
     &  - 4.D0*PC2*qq*ssm*dsvp*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F25*F45r*c5*f4*u5**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F25*F44r*c4*f4*u4*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F25*F43r*c3*f4*u3*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F25*F42r*c2*f4*u2*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F25*F41r*c1*f4*u1*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F25*F40r*c0*f4*u0*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F24*F45r*c5*f4*u4*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F24*F44r*c4*f4*u4**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F24*F43r*c3*f4*u3*u4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dsvp*F24*F42r*c2*f4*u2*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F24*F41r*c1*f4*u1*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F24*F40r*c0*f4*u0*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F23*F45r*c5*f4*u3*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F23*F44r*c4*f4*u3*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F23*F43r*c3*f4*u3**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F23*F42r*c2*f4*u2*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F23*F41r*c1*f4*u1*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F23*F40r*c0*f4*u0*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F22*F45r*c5*f4*u2*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F22*F44r*c4*f4*u2*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F22*F43r*c3*f4*u2*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F22*F42r*c2*f4*u2**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F22*F41r*c1*f4*u1*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F22*F40r*c0*f4*u0*u2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dsvp*F21*F45r*c5*f4*u1*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F21*F44r*c4*f4*u1*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F21*F43r*c3*f4*u1*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F21*F42r*c2*f4*u1*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F21*F41r*c1*f4*u1**2
     &  - 4.D0*PC2*qq*ssm*dsvp*F21*F40r*c0*f4*u0*u1
     &  - 4.D0*PC2*qq*ssm*dsvp*F20*F45r*c5*f4*u0*u5
     &  - 4.D0*PC2*qq*ssm*dsvp*F20*F44r*c4*f4*u0*u4
     &  - 4.D0*PC2*qq*ssm*dsvp*F20*F43r*c3*f4*u0*u3
     &  - 4.D0*PC2*qq*ssm*dsvp*F20*F42r*c2*f4*u0*u2
     &  - 4.D0*PC2*qq*ssm*dsvp*F20*F41r*c1*f4*u0*u1
     &  - 4.D0*PC2*qq*ssm*dsvp*F20*F40r*c0*f4*u0**2
     &  - 4.D0*PC2*qq*ssm*dssp*F45*F45r*c5*f4*u5**2
     &  - 4.D0*PC2*qq*ssm*dssp*F45*F44r*c4*f4*u4*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F45*F43r*c3*f4*u3*u5
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dssp*F45*F42r*c2*f4*u2*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F45*F41r*c1*f4*u1*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F45*F40r*c0*f4*u0*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F44*F45r*c5*f4*u4*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F44*F44r*c4*f4*u4**2
     &  - 4.D0*PC2*qq*ssm*dssp*F44*F43r*c3*f4*u3*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F44*F42r*c2*f4*u2*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F44*F41r*c1*f4*u1*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F44*F40r*c0*f4*u0*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F43*F45r*c5*f4*u3*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F43*F44r*c4*f4*u3*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F43*F43r*c3*f4*u3**2
     &  - 4.D0*PC2*qq*ssm*dssp*F43*F42r*c2*f4*u2*u3
     &  - 4.D0*PC2*qq*ssm*dssp*F43*F41r*c1*f4*u1*u3
     &  - 4.D0*PC2*qq*ssm*dssp*F43*F40r*c0*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dssp*F42*F45r*c5*f4*u2*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F42*F44r*c4*f4*u2*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F42*F43r*c3*f4*u2*u3
     &  - 4.D0*PC2*qq*ssm*dssp*F42*F42r*c2*f4*u2**2
     &  - 4.D0*PC2*qq*ssm*dssp*F42*F41r*c1*f4*u1*u2
     &  - 4.D0*PC2*qq*ssm*dssp*F42*F40r*c0*f4*u0*u2
     &  - 4.D0*PC2*qq*ssm*dssp*F41*F45r*c5*f4*u1*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F41*F44r*c4*f4*u1*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F41*F43r*c3*f4*u1*u3
     &  - 4.D0*PC2*qq*ssm*dssp*F41*F42r*c2*f4*u1*u2
     &  - 4.D0*PC2*qq*ssm*dssp*F41*F41r*c1*f4*u1**2
     &  - 4.D0*PC2*qq*ssm*dssp*F41*F40r*c0*f4*u0*u1
     &  - 4.D0*PC2*qq*ssm*dssp*F40*F45r*c5*f4*u0*u5
     &  - 4.D0*PC2*qq*ssm*dssp*F40*F44r*c4*f4*u0*u4
     &  - 4.D0*PC2*qq*ssm*dssp*F40*F43r*c3*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssm*dssp*F40*F42r*c2*f4*u0*u2
     &  - 4.D0*PC2*qq*ssm*dssp*F40*F41r*c1*f4*u0*u1
     &  - 4.D0*PC2*qq*ssm*dssp*F40*F40r*c0*f4*u0**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F45*F25r*c5*f2*u5**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F45*F24r*c4*f2*u4*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F45*F23r*c3*f2*u3*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F45*F22r*c2*f2*u2*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F45*F21r*c1*f2*u1*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F45*F20r*c0*f2*u0*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F44*F25r*c5*f2*u4*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F44*F24r*c4*f2*u4**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F44*F23r*c3*f2*u3*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F44*F22r*c2*f2*u2*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F44*F21r*c1*f2*u1*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F44*F20r*c0*f2*u0*u4
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dsvm*F43*F25r*c5*f2*u3*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F43*F24r*c4*f2*u3*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F43*F23r*c3*f2*u3**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F43*F22r*c2*f2*u2*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F43*F21r*c1*f2*u1*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F43*F20r*c0*f2*u0*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F42*F25r*c5*f2*u2*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F42*F24r*c4*f2*u2*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F42*F23r*c3*f2*u2*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F42*F22r*c2*f2*u2**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F42*F21r*c1*f2*u1*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F42*F20r*c0*f2*u0*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F41*F25r*c5*f2*u1*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F41*F24r*c4*f2*u1*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F41*F23r*c3*f2*u1*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dsvm*F41*F22r*c2*f2*u1*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F41*F21r*c1*f2*u1**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F41*F20r*c0*f2*u0*u1
     &  - 4.D0*PC2*qq*ssp*dsvm*F40*F25r*c5*f2*u0*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F40*F24r*c4*f2*u0*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F40*F23r*c3*f2*u0*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F40*F22r*c2*f2*u0*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F40*F21r*c1*f2*u0*u1
     &  - 4.D0*PC2*qq*ssp*dsvm*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F25*F45r*c5*f4*u5**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F25*F44r*c4*f4*u4*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F25*F43r*c3*f4*u3*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F25*F42r*c2*f4*u2*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F25*F41r*c1*f4*u1*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F25*F40r*c0*f4*u0*u5
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dsvm*F24*F45r*c5*f4*u4*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F24*F44r*c4*f4*u4**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F24*F43r*c3*f4*u3*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F24*F42r*c2*f4*u2*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F24*F41r*c1*f4*u1*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F24*F40r*c0*f4*u0*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F23*F45r*c5*f4*u3*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F23*F44r*c4*f4*u3*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F23*F43r*c3*f4*u3**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F23*F42r*c2*f4*u2*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F23*F41r*c1*f4*u1*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F23*F40r*c0*f4*u0*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F22*F45r*c5*f4*u2*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F22*F44r*c4*f4*u2*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F22*F43r*c3*f4*u2*u3
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dsvm*F22*F42r*c2*f4*u2**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F22*F41r*c1*f4*u1*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F22*F40r*c0*f4*u0*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F21*F45r*c5*f4*u1*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F21*F44r*c4*f4*u1*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F21*F43r*c3*f4*u1*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F21*F42r*c2*f4*u1*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F21*F41r*c1*f4*u1**2
     &  - 4.D0*PC2*qq*ssp*dsvm*F21*F40r*c0*f4*u0*u1
     &  - 4.D0*PC2*qq*ssp*dsvm*F20*F45r*c5*f4*u0*u5
     &  - 4.D0*PC2*qq*ssp*dsvm*F20*F44r*c4*f4*u0*u4
     &  - 4.D0*PC2*qq*ssp*dsvm*F20*F43r*c3*f4*u0*u3
     &  - 4.D0*PC2*qq*ssp*dsvm*F20*F42r*c2*f4*u0*u2
     &  - 4.D0*PC2*qq*ssp*dsvm*F20*F41r*c1*f4*u0*u1
     &  - 4.D0*PC2*qq*ssp*dsvm*F20*F40r*c0*f4*u0**2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dssm*F45*F45r*c5*f4*u5**2
     &  - 4.D0*PC2*qq*ssp*dssm*F45*F44r*c4*f4*u4*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F45*F43r*c3*f4*u3*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F45*F42r*c2*f4*u2*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F45*F41r*c1*f4*u1*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F45*F40r*c0*f4*u0*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F44*F45r*c5*f4*u4*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F44*F44r*c4*f4*u4**2
     &  - 4.D0*PC2*qq*ssp*dssm*F44*F43r*c3*f4*u3*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F44*F42r*c2*f4*u2*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F44*F41r*c1*f4*u1*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F44*F40r*c0*f4*u0*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F43*F45r*c5*f4*u3*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F43*F44r*c4*f4*u3*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F43*F43r*c3*f4*u3**2
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dssm*F43*F42r*c2*f4*u2*u3
     &  - 4.D0*PC2*qq*ssp*dssm*F43*F41r*c1*f4*u1*u3
     &  - 4.D0*PC2*qq*ssp*dssm*F43*F40r*c0*f4*u0*u3
     &  - 4.D0*PC2*qq*ssp*dssm*F42*F45r*c5*f4*u2*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F42*F44r*c4*f4*u2*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F42*F43r*c3*f4*u2*u3
     &  - 4.D0*PC2*qq*ssp*dssm*F42*F42r*c2*f4*u2**2
     &  - 4.D0*PC2*qq*ssp*dssm*F42*F41r*c1*f4*u1*u2
     &  - 4.D0*PC2*qq*ssp*dssm*F42*F40r*c0*f4*u0*u2
     &  - 4.D0*PC2*qq*ssp*dssm*F41*F45r*c5*f4*u1*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F41*F44r*c4*f4*u1*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F41*F43r*c3*f4*u1*u3
     &  - 4.D0*PC2*qq*ssp*dssm*F41*F42r*c2*f4*u1*u2
     &  - 4.D0*PC2*qq*ssp*dssm*F41*F41r*c1*f4*u1**2
     &  - 4.D0*PC2*qq*ssp*dssm*F41*F40r*c0*f4*u0*u1
     &
      traza1 = traza1 - 4.D0*PC2*qq*ssp*dssm*F40*F45r*c5*f4*u0*u5
     &  - 4.D0*PC2*qq*ssp*dssm*F40*F44r*c4*f4*u0*u4
     &  - 4.D0*PC2*qq*ssp*dssm*F40*F43r*c3*f4*u0*u3
     &  - 4.D0*PC2*qq*ssp*dssm*F40*F42r*c2*f4*u0*u2
     &  - 4.D0*PC2*qq*ssp*dssm*F40*F41r*c1*f4*u0*u1
     &  - 4.D0*PC2*qq*ssp*dssm*F40*F40r*c0*f4*u0**2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F45*F45r*c5*f4*u5**2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F45*F44r*c4*f4*u4*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F45*F43r*c3*f4*u3*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F45*F42r*c2*f4*u2*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F45*F41r*c1*f4*u1*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F45*F40r*c0*f4*u0*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F44*F45r*c5*f4*u4*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F44*F44r*c4*f4*u4**2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F44*F43r*c3*f4*u3*u4
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svm*dsvp*F44*F42r*c2*f4*u2*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F44*F41r*c1*f4*u1*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F44*F40r*c0*f4*u0*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F43*F45r*c5*f4*u3*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F43*F44r*c4*f4*u3*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F43*F43r*c3*f4*u3**2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F43*F42r*c2*f4*u2*u3
     &  + 4.D0*PC2*qq**2*svm*dsvp*F43*F41r*c1*f4*u1*u3
     &  + 4.D0*PC2*qq**2*svm*dsvp*F43*F40r*c0*f4*u0*u3
     &  + 4.D0*PC2*qq**2*svm*dsvp*F42*F45r*c5*f4*u2*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F42*F44r*c4*f4*u2*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F42*F43r*c3*f4*u2*u3
     &  + 4.D0*PC2*qq**2*svm*dsvp*F42*F42r*c2*f4*u2**2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F42*F41r*c1*f4*u1*u2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F42*F40r*c0*f4*u0*u2
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svm*dsvp*F41*F45r*c5*f4*u1*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F41*F44r*c4*f4*u1*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F41*F43r*c3*f4*u1*u3
     &  + 4.D0*PC2*qq**2*svm*dsvp*F41*F42r*c2*f4*u1*u2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F41*F41r*c1*f4*u1**2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F41*F40r*c0*f4*u0*u1
     &  + 4.D0*PC2*qq**2*svm*dsvp*F40*F45r*c5*f4*u0*u5
     &  + 4.D0*PC2*qq**2*svm*dsvp*F40*F44r*c4*f4*u0*u4
     &  + 4.D0*PC2*qq**2*svm*dsvp*F40*F43r*c3*f4*u0*u3
     &  + 4.D0*PC2*qq**2*svm*dsvp*F40*F42r*c2*f4*u0*u2
     &  + 4.D0*PC2*qq**2*svm*dsvp*F40*F41r*c1*f4*u0*u1
     &  + 4.D0*PC2*qq**2*svm*dsvp*F40*F40r*c0*f4*u0**2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F45*F45r*c5*f4*u5**2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F45*F44r*c4*f4*u4*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F45*F43r*c3*f4*u3*u5
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svp*dsvm*F45*F42r*c2*f4*u2*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F45*F41r*c1*f4*u1*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F45*F40r*c0*f4*u0*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F44*F45r*c5*f4*u4*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F44*F44r*c4*f4*u4**2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F44*F43r*c3*f4*u3*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F44*F42r*c2*f4*u2*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F44*F41r*c1*f4*u1*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F44*F40r*c0*f4*u0*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F43*F45r*c5*f4*u3*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F43*F44r*c4*f4*u3*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F43*F43r*c3*f4*u3**2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F43*F42r*c2*f4*u2*u3
     &  + 4.D0*PC2*qq**2*svp*dsvm*F43*F41r*c1*f4*u1*u3
     &  + 4.D0*PC2*qq**2*svp*dsvm*F43*F40r*c0*f4*u0*u3
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svp*dsvm*F42*F45r*c5*f4*u2*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F42*F44r*c4*f4*u2*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F42*F43r*c3*f4*u2*u3
     &  + 4.D0*PC2*qq**2*svp*dsvm*F42*F42r*c2*f4*u2**2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F42*F41r*c1*f4*u1*u2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F42*F40r*c0*f4*u0*u2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F41*F45r*c5*f4*u1*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F41*F44r*c4*f4*u1*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F41*F43r*c3*f4*u1*u3
     &  + 4.D0*PC2*qq**2*svp*dsvm*F41*F42r*c2*f4*u1*u2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F41*F41r*c1*f4*u1**2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F41*F40r*c0*f4*u0*u1
     &  + 4.D0*PC2*qq**2*svp*dsvm*F40*F45r*c5*f4*u0*u5
     &  + 4.D0*PC2*qq**2*svp*dsvm*F40*F44r*c4*f4*u0*u4
     &  + 4.D0*PC2*qq**2*svp*dsvm*F40*F43r*c3*f4*u0*u3
     &
      traza1 = traza1 + 4.D0*PC2*qq**2*svp*dsvm*F40*F42r*c2*f4*u0*u2
     &  + 4.D0*PC2*qq**2*svp*dsvm*F40*F41r*c1*f4*u0*u1
     &  + 4.D0*PC2*qq**2*svp*dsvm*F40*F40r*c0*f4*u0**2
     &  - 4.D0*PC2**2*svm*dsvp*F25*F25r*c5*f2*u5**2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F25*F24r*c4*f2*u4*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F25*F23r*c3*f2*u3*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F25*F22r*c2*f2*u2*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F25*F21r*c1*f2*u1*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F25*F20r*c0*f2*u0*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F24*F25r*c5*f2*u4*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F24*F24r*c4*f2*u4**2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F24*F23r*c3*f2*u3*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F24*F22r*c2*f2*u2*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F24*F21r*c1*f2*u1*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F24*F20r*c0*f2*u0*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*svm*dsvp*F23*F25r*c5*f2*u3*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F23*F24r*c4*f2*u3*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F23*F23r*c3*f2*u3**2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F23*F22r*c2*f2*u2*u3*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F23*F21r*c1*f2*u1*u3*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F23*F20r*c0*f2*u0*u3*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F22*F25r*c5*f2*u2*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F22*F24r*c4*f2*u2*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F22*F23r*c3*f2*u2*u3*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F22*F22r*c2*f2*u2**2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F22*F21r*c1*f2*u1*u2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F22*F20r*c0*f2*u0*u2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F21*F25r*c5*f2*u1*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F21*F24r*c4*f2*u1*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F21*F23r*c3*f2*u1*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*svm*dsvp*F21*F22r*c2*f2*u1*u2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F21*F21r*c1*f2*u1**2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F21*F20r*c0*f2*u0*u1*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F20*F25r*c5*f2*u0*u5*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F20*F24r*c4*f2*u0*u4*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F20*F23r*c3*f2*u0*u3*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F20*F22r*c2*f2*u0*u2*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F20*F21r*c1*f2*u0*u1*w1*w2
     &  - 4.D0*PC2**2*svm*dsvp*F20*F20r*c0*f2*u0**2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F25*F25r*c5*f2*u5**2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F25*F24r*c4*f2*u4*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F25*F23r*c3*f2*u3*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F25*F22r*c2*f2*u2*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F25*F21r*c1*f2*u1*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F25*F20r*c0*f2*u0*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*svp*dsvm*F24*F25r*c5*f2*u4*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F24*F24r*c4*f2*u4**2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F24*F23r*c3*f2*u3*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F24*F22r*c2*f2*u2*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F24*F21r*c1*f2*u1*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F24*F20r*c0*f2*u0*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F23*F25r*c5*f2*u3*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F23*F24r*c4*f2*u3*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F23*F23r*c3*f2*u3**2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F23*F22r*c2*f2*u2*u3*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F23*F21r*c1*f2*u1*u3*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F23*F20r*c0*f2*u0*u3*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F22*F25r*c5*f2*u2*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F22*F24r*c4*f2*u2*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F22*F23r*c3*f2*u2*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*svp*dsvm*F22*F22r*c2*f2*u2**2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F22*F21r*c1*f2*u1*u2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F22*F20r*c0*f2*u0*u2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F21*F25r*c5*f2*u1*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F21*F24r*c4*f2*u1*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F21*F23r*c3*f2*u1*u3*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F21*F22r*c2*f2*u1*u2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F21*F21r*c1*f2*u1**2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F21*F20r*c0*f2*u0*u1*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F20*F25r*c5*f2*u0*u5*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F20*F24r*c4*f2*u0*u4*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F20*F23r*c3*f2*u0*u3*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F20*F22r*c2*f2*u0*u2*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F20*F21r*c1*f2*u0*u1*w1*w2
     &  - 4.D0*PC2**2*svp*dsvm*F20*F20r*c0*f2*u0**2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svm*dsvp*F45*F45r*c5*f4*u5**2*w1
     & *w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F45*F44r*c4*f4*u4*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F45*F43r*c3*f4*u3*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F45*F42r*c2*f4*u2*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F45*F41r*c1*f4*u1*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F45*F40r*c0*f4*u0*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F44*F45r*c5*f4*u4*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F44*F44r*c4*f4*u4**2*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F44*F43r*c3*f4*u3*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F44*F42r*c2*f4*u2*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F44*F41r*c1*f4*u1*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F44*F40r*c0*f4*u0*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F43*F45r*c5*f4*u3*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F43*F44r*c4*f4*u3*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svm*dsvp*F43*F43r*c3*f4*u3**2*w1
     & *w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F43*F42r*c2*f4*u2*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F43*F41r*c1*f4*u1*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F43*F40r*c0*f4*u0*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F42*F45r*c5*f4*u2*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F42*F44r*c4*f4*u2*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F42*F43r*c3*f4*u2*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F42*F42r*c2*f4*u2**2*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F42*F41r*c1*f4*u1*u2*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F42*F40r*c0*f4*u0*u2*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F41*F45r*c5*f4*u1*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F41*F44r*c4*f4*u1*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F41*F43r*c3*f4*u1*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F41*F42r*c2*f4*u1*u2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svm*dsvp*F41*F41r*c1*f4*u1**2*w1
     & *w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F41*F40r*c0*f4*u0*u1*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F40*F45r*c5*f4*u0*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F40*F44r*c4*f4*u0*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F40*F43r*c3*f4*u0*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F40*F42r*c2*f4*u0*u2*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F40*F41r*c1*f4*u0*u1*w1*w2
     &  - 4.D0*PC2**2*qq*svm*dsvp*F40*F40r*c0*f4*u0**2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F45*F45r*c5*f4*u5**2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F45*F44r*c4*f4*u4*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F45*F43r*c3*f4*u3*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F45*F42r*c2*f4*u2*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F45*F41r*c1*f4*u1*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F45*F40r*c0*f4*u0*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svp*dsvm*F44*F45r*c5*f4*u4*u5*w1
     & *w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F44*F44r*c4*f4*u4**2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F44*F43r*c3*f4*u3*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F44*F42r*c2*f4*u2*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F44*F41r*c1*f4*u1*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F44*F40r*c0*f4*u0*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F43*F45r*c5*f4*u3*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F43*F44r*c4*f4*u3*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F43*F43r*c3*f4*u3**2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F43*F42r*c2*f4*u2*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F43*F41r*c1*f4*u1*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F43*F40r*c0*f4*u0*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F42*F45r*c5*f4*u2*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F42*F44r*c4*f4*u2*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svp*dsvm*F42*F43r*c3*f4*u2*u3*w1
     & *w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F42*F42r*c2*f4*u2**2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F42*F41r*c1*f4*u1*u2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F42*F40r*c0*f4*u0*u2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F41*F45r*c5*f4*u1*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F41*F44r*c4*f4*u1*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F41*F43r*c3*f4*u1*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F41*F42r*c2*f4*u1*u2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F41*F41r*c1*f4*u1**2*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F41*F40r*c0*f4*u0*u1*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F40*F45r*c5*f4*u0*u5*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F40*F44r*c4*f4*u0*u4*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F40*F43r*c3*f4*u0*u3*w1*w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F40*F42r*c2*f4*u0*u2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC2**2*qq*svp*dsvm*F40*F41r*c1*f4*u0*u1*w1
     & *w2
     &  - 4.D0*PC2**2*qq*svp*dsvm*F40*F40r*c0*f4*u0**2*w1*w2
     &  - 4.D0*i_*ssm*dssp*F15*F15i*c5*f1*u5**2
     &  - 4.D0*i_*ssm*dssp*F15*F14i*c4*f1*u4*u5
     &  - 4.D0*i_*ssm*dssp*F15*F13i*c3*f1*u3*u5
     &  - 4.D0*i_*ssm*dssp*F15*F12i*c2*f1*u2*u5
     &  - 4.D0*i_*ssm*dssp*F15*F11i*c1*f1*u1*u5
     &  - 4.D0*i_*ssm*dssp*F15*F10i*c0*f1*u0*u5
     &  - 4.D0*i_*ssm*dssp*F14*F15i*c5*f1*u4*u5
     &  - 4.D0*i_*ssm*dssp*F14*F14i*c4*f1*u4**2
     &  - 4.D0*i_*ssm*dssp*F14*F13i*c3*f1*u3*u4
     &  - 4.D0*i_*ssm*dssp*F14*F12i*c2*f1*u2*u4
     &  - 4.D0*i_*ssm*dssp*F14*F11i*c1*f1*u1*u4
     &  - 4.D0*i_*ssm*dssp*F14*F10i*c0*f1*u0*u4
     &
      traza1 = traza1 - 4.D0*i_*ssm*dssp*F13*F15i*c5*f1*u3*u5
     &  - 4.D0*i_*ssm*dssp*F13*F14i*c4*f1*u3*u4
     &  - 4.D0*i_*ssm*dssp*F13*F13i*c3*f1*u3**2
     &  - 4.D0*i_*ssm*dssp*F13*F12i*c2*f1*u2*u3
     &  - 4.D0*i_*ssm*dssp*F13*F11i*c1*f1*u1*u3
     &  - 4.D0*i_*ssm*dssp*F13*F10i*c0*f1*u0*u3
     &  - 4.D0*i_*ssm*dssp*F12*F15i*c5*f1*u2*u5
     &  - 4.D0*i_*ssm*dssp*F12*F14i*c4*f1*u2*u4
     &  - 4.D0*i_*ssm*dssp*F12*F13i*c3*f1*u2*u3
     &  - 4.D0*i_*ssm*dssp*F12*F12i*c2*f1*u2**2
     &  - 4.D0*i_*ssm*dssp*F12*F11i*c1*f1*u1*u2
     &  - 4.D0*i_*ssm*dssp*F12*F10i*c0*f1*u0*u2
     &  - 4.D0*i_*ssm*dssp*F11*F15i*c5*f1*u1*u5
     &  - 4.D0*i_*ssm*dssp*F11*F14i*c4*f1*u1*u4
     &  - 4.D0*i_*ssm*dssp*F11*F13i*c3*f1*u1*u3
     &
      traza1 = traza1 - 4.D0*i_*ssm*dssp*F11*F12i*c2*f1*u1*u2
     &  - 4.D0*i_*ssm*dssp*F11*F11i*c1*f1*u1**2
     &  - 4.D0*i_*ssm*dssp*F11*F10i*c0*f1*u0*u1
     &  - 4.D0*i_*ssm*dssp*F10*F15i*c5*f1*u0*u5
     &  - 4.D0*i_*ssm*dssp*F10*F14i*c4*f1*u0*u4
     &  - 4.D0*i_*ssm*dssp*F10*F13i*c3*f1*u0*u3
     &  - 4.D0*i_*ssm*dssp*F10*F12i*c2*f1*u0*u2
     &  - 4.D0*i_*ssm*dssp*F10*F11i*c1*f1*u0*u1
     &  - 4.D0*i_*ssm*dssp*F10*F10i*c0*f1*u0**2
     &  - 4.D0*i_*ssp*dssm*F15*F15i*c5*f1*u5**2
     &  - 4.D0*i_*ssp*dssm*F15*F14i*c4*f1*u4*u5
     &  - 4.D0*i_*ssp*dssm*F15*F13i*c3*f1*u3*u5
     &  - 4.D0*i_*ssp*dssm*F15*F12i*c2*f1*u2*u5
     &  - 4.D0*i_*ssp*dssm*F15*F11i*c1*f1*u1*u5
     &  - 4.D0*i_*ssp*dssm*F15*F10i*c0*f1*u0*u5
     &
      traza1 = traza1 - 4.D0*i_*ssp*dssm*F14*F15i*c5*f1*u4*u5
     &  - 4.D0*i_*ssp*dssm*F14*F14i*c4*f1*u4**2
     &  - 4.D0*i_*ssp*dssm*F14*F13i*c3*f1*u3*u4
     &  - 4.D0*i_*ssp*dssm*F14*F12i*c2*f1*u2*u4
     &  - 4.D0*i_*ssp*dssm*F14*F11i*c1*f1*u1*u4
     &  - 4.D0*i_*ssp*dssm*F14*F10i*c0*f1*u0*u4
     &  - 4.D0*i_*ssp*dssm*F13*F15i*c5*f1*u3*u5
     &  - 4.D0*i_*ssp*dssm*F13*F14i*c4*f1*u3*u4
     &  - 4.D0*i_*ssp*dssm*F13*F13i*c3*f1*u3**2
     &  - 4.D0*i_*ssp*dssm*F13*F12i*c2*f1*u2*u3
     &  - 4.D0*i_*ssp*dssm*F13*F11i*c1*f1*u1*u3
     &  - 4.D0*i_*ssp*dssm*F13*F10i*c0*f1*u0*u3
     &  - 4.D0*i_*ssp*dssm*F12*F15i*c5*f1*u2*u5
     &  - 4.D0*i_*ssp*dssm*F12*F14i*c4*f1*u2*u4
     &  - 4.D0*i_*ssp*dssm*F12*F13i*c3*f1*u2*u3
     &
      traza1 = traza1 - 4.D0*i_*ssp*dssm*F12*F12i*c2*f1*u2**2
     &  - 4.D0*i_*ssp*dssm*F12*F11i*c1*f1*u1*u2
     &  - 4.D0*i_*ssp*dssm*F12*F10i*c0*f1*u0*u2
     &  - 4.D0*i_*ssp*dssm*F11*F15i*c5*f1*u1*u5
     &  - 4.D0*i_*ssp*dssm*F11*F14i*c4*f1*u1*u4
     &  - 4.D0*i_*ssp*dssm*F11*F13i*c3*f1*u1*u3
     &  - 4.D0*i_*ssp*dssm*F11*F12i*c2*f1*u1*u2
     &  - 4.D0*i_*ssp*dssm*F11*F11i*c1*f1*u1**2
     &  - 4.D0*i_*ssp*dssm*F11*F10i*c0*f1*u0*u1
     &  - 4.D0*i_*ssp*dssm*F10*F15i*c5*f1*u0*u5
     &  - 4.D0*i_*ssp*dssm*F10*F14i*c4*f1*u0*u4
     &  - 4.D0*i_*ssp*dssm*F10*F13i*c3*f1*u0*u3
     &  - 4.D0*i_*ssp*dssm*F10*F12i*c2*f1*u0*u2
     &  - 4.D0*i_*ssp*dssm*F10*F11i*c1*f1*u0*u1
     &  - 4.D0*i_*ssp*dssm*F10*F10i*c0*f1*u0**2
     &
      traza1 = traza1 - 4.D0*i_*qq*svm*dsvp*F15*F15i*c5*f1*u5**2
     &  - 4.D0*i_*qq*svm*dsvp*F15*F14i*c4*f1*u4*u5
     &  - 4.D0*i_*qq*svm*dsvp*F15*F13i*c3*f1*u3*u5
     &  - 4.D0*i_*qq*svm*dsvp*F15*F12i*c2*f1*u2*u5
     &  - 4.D0*i_*qq*svm*dsvp*F15*F11i*c1*f1*u1*u5
     &  - 4.D0*i_*qq*svm*dsvp*F15*F10i*c0*f1*u0*u5
     &  - 4.D0*i_*qq*svm*dsvp*F14*F15i*c5*f1*u4*u5
     &  - 4.D0*i_*qq*svm*dsvp*F14*F14i*c4*f1*u4**2
     &  - 4.D0*i_*qq*svm*dsvp*F14*F13i*c3*f1*u3*u4
     &  - 4.D0*i_*qq*svm*dsvp*F14*F12i*c2*f1*u2*u4
     &  - 4.D0*i_*qq*svm*dsvp*F14*F11i*c1*f1*u1*u4
     &  - 4.D0*i_*qq*svm*dsvp*F14*F10i*c0*f1*u0*u4
     &  - 4.D0*i_*qq*svm*dsvp*F13*F15i*c5*f1*u3*u5
     &  - 4.D0*i_*qq*svm*dsvp*F13*F14i*c4*f1*u3*u4
     &  - 4.D0*i_*qq*svm*dsvp*F13*F13i*c3*f1*u3**2
     &
      traza1 = traza1 - 4.D0*i_*qq*svm*dsvp*F13*F12i*c2*f1*u2*u3
     &  - 4.D0*i_*qq*svm*dsvp*F13*F11i*c1*f1*u1*u3
     &  - 4.D0*i_*qq*svm*dsvp*F13*F10i*c0*f1*u0*u3
     &  - 4.D0*i_*qq*svm*dsvp*F12*F15i*c5*f1*u2*u5
     &  - 4.D0*i_*qq*svm*dsvp*F12*F14i*c4*f1*u2*u4
     &  - 4.D0*i_*qq*svm*dsvp*F12*F13i*c3*f1*u2*u3
     &  - 4.D0*i_*qq*svm*dsvp*F12*F12i*c2*f1*u2**2
     &  - 4.D0*i_*qq*svm*dsvp*F12*F11i*c1*f1*u1*u2
     &  - 4.D0*i_*qq*svm*dsvp*F12*F10i*c0*f1*u0*u2
     &  - 4.D0*i_*qq*svm*dsvp*F11*F15i*c5*f1*u1*u5
     &  - 4.D0*i_*qq*svm*dsvp*F11*F14i*c4*f1*u1*u4
     &  - 4.D0*i_*qq*svm*dsvp*F11*F13i*c3*f1*u1*u3
     &  - 4.D0*i_*qq*svm*dsvp*F11*F12i*c2*f1*u1*u2
     &  - 4.D0*i_*qq*svm*dsvp*F11*F11i*c1*f1*u1**2
     &  - 4.D0*i_*qq*svm*dsvp*F11*F10i*c0*f1*u0*u1
     &
      traza1 = traza1 - 4.D0*i_*qq*svm*dsvp*F10*F15i*c5*f1*u0*u5
     &  - 4.D0*i_*qq*svm*dsvp*F10*F14i*c4*f1*u0*u4
     &  - 4.D0*i_*qq*svm*dsvp*F10*F13i*c3*f1*u0*u3
     &  - 4.D0*i_*qq*svm*dsvp*F10*F12i*c2*f1*u0*u2
     &  - 4.D0*i_*qq*svm*dsvp*F10*F11i*c1*f1*u0*u1
     &  - 4.D0*i_*qq*svm*dsvp*F10*F10i*c0*f1*u0**2
     &  - 4.D0*i_*qq*svp*dsvm*F15*F15i*c5*f1*u5**2
     &  - 4.D0*i_*qq*svp*dsvm*F15*F14i*c4*f1*u4*u5
     &  - 4.D0*i_*qq*svp*dsvm*F15*F13i*c3*f1*u3*u5
     &  - 4.D0*i_*qq*svp*dsvm*F15*F12i*c2*f1*u2*u5
     &  - 4.D0*i_*qq*svp*dsvm*F15*F11i*c1*f1*u1*u5
     &  - 4.D0*i_*qq*svp*dsvm*F15*F10i*c0*f1*u0*u5
     &  - 4.D0*i_*qq*svp*dsvm*F14*F15i*c5*f1*u4*u5
     &  - 4.D0*i_*qq*svp*dsvm*F14*F14i*c4*f1*u4**2
     &  - 4.D0*i_*qq*svp*dsvm*F14*F13i*c3*f1*u3*u4
     &
      traza1 = traza1 - 4.D0*i_*qq*svp*dsvm*F14*F12i*c2*f1*u2*u4
     &  - 4.D0*i_*qq*svp*dsvm*F14*F11i*c1*f1*u1*u4
     &  - 4.D0*i_*qq*svp*dsvm*F14*F10i*c0*f1*u0*u4
     &  - 4.D0*i_*qq*svp*dsvm*F13*F15i*c5*f1*u3*u5
     &  - 4.D0*i_*qq*svp*dsvm*F13*F14i*c4*f1*u3*u4
     &  - 4.D0*i_*qq*svp*dsvm*F13*F13i*c3*f1*u3**2
     &  - 4.D0*i_*qq*svp*dsvm*F13*F12i*c2*f1*u2*u3
     &  - 4.D0*i_*qq*svp*dsvm*F13*F11i*c1*f1*u1*u3
     &  - 4.D0*i_*qq*svp*dsvm*F13*F10i*c0*f1*u0*u3
     &  - 4.D0*i_*qq*svp*dsvm*F12*F15i*c5*f1*u2*u5
     &  - 4.D0*i_*qq*svp*dsvm*F12*F14i*c4*f1*u2*u4
     &  - 4.D0*i_*qq*svp*dsvm*F12*F13i*c3*f1*u2*u3
     &  - 4.D0*i_*qq*svp*dsvm*F12*F12i*c2*f1*u2**2
     &  - 4.D0*i_*qq*svp*dsvm*F12*F11i*c1*f1*u1*u2
     &  - 4.D0*i_*qq*svp*dsvm*F12*F10i*c0*f1*u0*u2
     &
      traza1 = traza1 - 4.D0*i_*qq*svp*dsvm*F11*F15i*c5*f1*u1*u5
     &  - 4.D0*i_*qq*svp*dsvm*F11*F14i*c4*f1*u1*u4
     &  - 4.D0*i_*qq*svp*dsvm*F11*F13i*c3*f1*u1*u3
     &  - 4.D0*i_*qq*svp*dsvm*F11*F12i*c2*f1*u1*u2
     &  - 4.D0*i_*qq*svp*dsvm*F11*F11i*c1*f1*u1**2
     &  - 4.D0*i_*qq*svp*dsvm*F11*F10i*c0*f1*u0*u1
     &  - 4.D0*i_*qq*svp*dsvm*F10*F15i*c5*f1*u0*u5
     &  - 4.D0*i_*qq*svp*dsvm*F10*F14i*c4*f1*u0*u4
     &  - 4.D0*i_*qq*svp*dsvm*F10*F13i*c3*f1*u0*u3
     &  - 4.D0*i_*qq*svp*dsvm*F10*F12i*c2*f1*u0*u2
     &  - 4.D0*i_*qq*svp*dsvm*F10*F11i*c1*f1*u0*u1
     &  - 4.D0*i_*qq*svp*dsvm*F10*F10i*c0*f1*u0**2
     &  + 4.D0*i_*PC2*svm*dsvp*F15*F15i*c5*f1*u5**2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F15*F14i*c4*f1*u4*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F15*F13i*c3*f1*u3*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svm*dsvp*F15*F12i*c2*f1*u2*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F15*F11i*c1*f1*u1*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F15*F10i*c0*f1*u0*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F14*F15i*c5*f1*u4*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F14*F14i*c4*f1*u4**2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F14*F13i*c3*f1*u3*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F14*F12i*c2*f1*u2*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F14*F11i*c1*f1*u1*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F14*F10i*c0*f1*u0*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F13*F15i*c5*f1*u3*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F13*F14i*c4*f1*u3*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F13*F13i*c3*f1*u3**2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F13*F12i*c2*f1*u2*u3*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F13*F11i*c1*f1*u1*u3*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F13*F10i*c0*f1*u0*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svm*dsvp*F12*F15i*c5*f1*u2*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F12*F14i*c4*f1*u2*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F12*F13i*c3*f1*u2*u3*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F12*F12i*c2*f1*u2**2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F12*F11i*c1*f1*u1*u2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F12*F10i*c0*f1*u0*u2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F11*F15i*c5*f1*u1*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F11*F14i*c4*f1*u1*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F11*F13i*c3*f1*u1*u3*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F11*F12i*c2*f1*u1*u2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F11*F11i*c1*f1*u1**2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F11*F10i*c0*f1*u0*u1*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F10*F15i*c5*f1*u0*u5*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F10*F14i*c4*f1*u0*u4*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F10*F13i*c3*f1*u0*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svm*dsvp*F10*F12i*c2*f1*u0*u2*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F10*F11i*c1*f1*u0*u1*w1*w2
     &  + 4.D0*i_*PC2*svm*dsvp*F10*F10i*c0*f1*u0**2*w1*w2
     &  - 4.D0*i_*PC2*svm*dssp*F25*F15i*c5*f1*u5**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F25*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F25*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F25*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F25*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F25*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F24*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F24*F14i*c4*f1*u4**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F24*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F24*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F24*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F24*F10i*c0*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*svm*dssp*F23*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F23*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F23*F13i*c3*f1*u3**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F23*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F23*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F23*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F22*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F22*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F22*F13i*c3*f1*u2*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F22*F12i*c2*f1*u2**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F22*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F22*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F21*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F21*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F21*F13i*c3*f1*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*svm*dssp*F21*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F21*F11i*c1*f1*u1**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F21*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*i_*PC2*svm*dssp*F20*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F20*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F20*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F20*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F20*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*i_*PC2*svm*dssp*F20*F10i*c0*f1*u0**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F15*F25i*c5*f2*u5**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F15*F24i*c4*f2*u4*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F15*F23i*c3*f2*u3*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F15*F22i*c2*f2*u2*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F15*F21i*c1*f2*u1*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F15*F20i*c0*f2*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*svm*dssp*F14*F25i*c5*f2*u4*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F14*F24i*c4*f2*u4**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F14*F23i*c3*f2*u3*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F14*F22i*c2*f2*u2*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F14*F21i*c1*f2*u1*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F14*F20i*c0*f2*u0*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F13*F25i*c5*f2*u3*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F13*F24i*c4*f2*u3*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F13*F23i*c3*f2*u3**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F13*F22i*c2*f2*u2*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F13*F21i*c1*f2*u1*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F13*F20i*c0*f2*u0*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F12*F25i*c5*f2*u2*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F12*F24i*c4*f2*u2*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F12*F23i*c3*f2*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*svm*dssp*F12*F22i*c2*f2*u2**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F12*F21i*c1*f2*u1*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F12*F20i*c0*f2*u0*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F11*F25i*c5*f2*u1*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F11*F24i*c4*f2*u1*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F11*F23i*c3*f2*u1*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F11*F22i*c2*f2*u1*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F11*F21i*c1*f2*u1**2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F11*F20i*c0*f2*u0*u1*w2
     &  - 4.D0*i_*PC2*svm*dssp*F10*F25i*c5*f2*u0*u5*w2
     &  - 4.D0*i_*PC2*svm*dssp*F10*F24i*c4*f2*u0*u4*w2
     &  - 4.D0*i_*PC2*svm*dssp*F10*F23i*c3*f2*u0*u3*w2
     &  - 4.D0*i_*PC2*svm*dssp*F10*F22i*c2*f2*u0*u2*w2
     &  - 4.D0*i_*PC2*svm*dssp*F10*F21i*c1*f2*u0*u1*w2
     &  - 4.D0*i_*PC2*svm*dssp*F10*F20i*c0*f2*u0**2*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svp*dsvm*F15*F15i*c5*f1*u5**2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F15*F14i*c4*f1*u4*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F15*F13i*c3*f1*u3*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F15*F12i*c2*f1*u2*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F15*F11i*c1*f1*u1*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F15*F10i*c0*f1*u0*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F14*F15i*c5*f1*u4*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F14*F14i*c4*f1*u4**2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F14*F13i*c3*f1*u3*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F14*F12i*c2*f1*u2*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F14*F11i*c1*f1*u1*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F14*F10i*c0*f1*u0*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F13*F15i*c5*f1*u3*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F13*F14i*c4*f1*u3*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F13*F13i*c3*f1*u3**2*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svp*dsvm*F13*F12i*c2*f1*u2*u3*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F13*F11i*c1*f1*u1*u3*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F13*F10i*c0*f1*u0*u3*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F12*F15i*c5*f1*u2*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F12*F14i*c4*f1*u2*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F12*F13i*c3*f1*u2*u3*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F12*F12i*c2*f1*u2**2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F12*F11i*c1*f1*u1*u2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F12*F10i*c0*f1*u0*u2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F11*F15i*c5*f1*u1*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F11*F14i*c4*f1*u1*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F11*F13i*c3*f1*u1*u3*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F11*F12i*c2*f1*u1*u2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F11*F11i*c1*f1*u1**2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F11*F10i*c0*f1*u0*u1*w1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*svp*dsvm*F10*F15i*c5*f1*u0*u5*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F10*F14i*c4*f1*u0*u4*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F10*F13i*c3*f1*u0*u3*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F10*F12i*c2*f1*u0*u2*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F10*F11i*c1*f1*u0*u1*w1*w2
     &  + 4.D0*i_*PC2*svp*dsvm*F10*F10i*c0*f1*u0**2*w1*w2
     &  - 4.D0*i_*PC2*svp*dssm*F25*F15i*c5*f1*u5**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F25*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F25*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F25*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F25*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F25*F10i*c0*f1*u0*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F24*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F24*F14i*c4*f1*u4**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F24*F13i*c3*f1*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*svp*dssm*F24*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F24*F11i*c1*f1*u1*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F24*F10i*c0*f1*u0*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F23*F15i*c5*f1*u3*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F23*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F23*F13i*c3*f1*u3**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F23*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F23*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F23*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F22*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F22*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F22*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F22*F12i*c2*f1*u2**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F22*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F22*F10i*c0*f1*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*svp*dssm*F21*F15i*c5*f1*u1*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F21*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F21*F13i*c3*f1*u1*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F21*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F21*F11i*c1*f1*u1**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F21*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*i_*PC2*svp*dssm*F20*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F20*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F20*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F20*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F20*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*i_*PC2*svp*dssm*F20*F10i*c0*f1*u0**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F15*F25i*c5*f2*u5**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F15*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F15*F23i*c3*f2*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*svp*dssm*F15*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F15*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F15*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F14*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F14*F24i*c4*f2*u4**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F14*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F14*F22i*c2*f2*u2*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F14*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F14*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F13*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F13*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F13*F23i*c3*f2*u3**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F13*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F13*F21i*c1*f2*u1*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F13*F20i*c0*f2*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*svp*dssm*F12*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F12*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F12*F23i*c3*f2*u2*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F12*F22i*c2*f2*u2**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F12*F21i*c1*f2*u1*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F12*F20i*c0*f2*u0*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F11*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F11*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F11*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*i_*PC2*svp*dssm*F11*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F11*F21i*c1*f2*u1**2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F11*F20i*c0*f2*u0*u1*w1
     &  - 4.D0*i_*PC2*svp*dssm*F10*F25i*c5*f2*u0*u5*w1
     &  - 4.D0*i_*PC2*svp*dssm*F10*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*i_*PC2*svp*dssm*F10*F23i*c3*f2*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*svp*dssm*F10*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*i_*PC2*svp*dssm*F10*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*i_*PC2*svp*dssm*F10*F20i*c0*f2*u0**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F25*F15i*c5*f1*u5**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F25*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F25*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F25*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F25*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F25*F10i*c0*f1*u0*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F24*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F24*F14i*c4*f1*u4**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F24*F13i*c3*f1*u3*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F24*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F24*F11i*c1*f1*u1*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F24*F10i*c0*f1*u0*u4*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*dsvp*F23*F15i*c5*f1*u3*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F23*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F23*F13i*c3*f1*u3**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F23*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F23*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F23*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F22*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F22*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F22*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F22*F12i*c2*f1*u2**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F22*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F22*F10i*c0*f1*u0*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F21*F15i*c5*f1*u1*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F21*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F21*F13i*c3*f1*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*dsvp*F21*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F21*F11i*c1*f1*u1**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F21*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F20*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F20*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F20*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F20*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F20*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F20*F10i*c0*f1*u0**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F15*F25i*c5*f2*u5**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F15*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F15*F23i*c3*f2*u3*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F15*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F15*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F15*F20i*c0*f2*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*dsvp*F14*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F14*F24i*c4*f2*u4**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F14*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F14*F22i*c2*f2*u2*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F14*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F14*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F13*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F13*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F13*F23i*c3*f2*u3**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F13*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F13*F21i*c1*f2*u1*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F13*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F12*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F12*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F12*F23i*c3*f2*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssm*dsvp*F12*F22i*c2*f2*u2**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F12*F21i*c1*f2*u1*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F12*F20i*c0*f2*u0*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F11*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F11*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F11*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F11*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F11*F21i*c1*f2*u1**2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F11*F20i*c0*f2*u0*u1*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F10*F25i*c5*f2*u0*u5*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F10*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F10*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F10*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F10*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*i_*PC2*ssm*dsvp*F10*F20i*c0*f2*u0**2*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssm*dssp*F25*F25i*c5*f2*u5**2
     &  + 4.D0*i_*PC2*ssm*dssp*F25*F24i*c4*f2*u4*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F25*F23i*c3*f2*u3*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F25*F22i*c2*f2*u2*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F25*F21i*c1*f2*u1*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F25*F20i*c0*f2*u0*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F24*F25i*c5*f2*u4*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F24*F24i*c4*f2*u4**2
     &  + 4.D0*i_*PC2*ssm*dssp*F24*F23i*c3*f2*u3*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F24*F22i*c2*f2*u2*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F24*F21i*c1*f2*u1*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F24*F20i*c0*f2*u0*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F23*F25i*c5*f2*u3*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F23*F24i*c4*f2*u3*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F23*F23i*c3*f2*u3**2
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssm*dssp*F23*F22i*c2*f2*u2*u3
     &  + 4.D0*i_*PC2*ssm*dssp*F23*F21i*c1*f2*u1*u3
     &  + 4.D0*i_*PC2*ssm*dssp*F23*F20i*c0*f2*u0*u3
     &  + 4.D0*i_*PC2*ssm*dssp*F22*F25i*c5*f2*u2*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F22*F24i*c4*f2*u2*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F22*F23i*c3*f2*u2*u3
     &  + 4.D0*i_*PC2*ssm*dssp*F22*F22i*c2*f2*u2**2
     &  + 4.D0*i_*PC2*ssm*dssp*F22*F21i*c1*f2*u1*u2
     &  + 4.D0*i_*PC2*ssm*dssp*F22*F20i*c0*f2*u0*u2
     &  + 4.D0*i_*PC2*ssm*dssp*F21*F25i*c5*f2*u1*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F21*F24i*c4*f2*u1*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F21*F23i*c3*f2*u1*u3
     &  + 4.D0*i_*PC2*ssm*dssp*F21*F22i*c2*f2*u1*u2
     &  + 4.D0*i_*PC2*ssm*dssp*F21*F21i*c1*f2*u1**2
     &  + 4.D0*i_*PC2*ssm*dssp*F21*F20i*c0*f2*u0*u1
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssm*dssp*F20*F25i*c5*f2*u0*u5
     &  + 4.D0*i_*PC2*ssm*dssp*F20*F24i*c4*f2*u0*u4
     &  + 4.D0*i_*PC2*ssm*dssp*F20*F23i*c3*f2*u0*u3
     &  + 4.D0*i_*PC2*ssm*dssp*F20*F22i*c2*f2*u0*u2
     &  + 4.D0*i_*PC2*ssm*dssp*F20*F21i*c1*f2*u0*u1
     &  + 4.D0*i_*PC2*ssm*dssp*F20*F20i*c0*f2*u0**2
     &  - 4.D0*i_*PC2*ssp*dsvm*F25*F15i*c5*f1*u5**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F25*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F25*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F25*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F25*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F25*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F24*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F24*F14i*c4*f1*u4**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F24*F13i*c3*f1*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*dsvm*F24*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F24*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F24*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F23*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F23*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F23*F13i*c3*f1*u3**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F23*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F23*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F23*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F22*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F22*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F22*F13i*c3*f1*u2*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F22*F12i*c2*f1*u2**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F22*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F22*F10i*c0*f1*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*dsvm*F21*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F21*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F21*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F21*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F21*F11i*c1*f1*u1**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F21*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F20*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F20*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F20*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F20*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F20*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F20*F10i*c0*f1*u0**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F15*F25i*c5*f2*u5**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F15*F24i*c4*f2*u4*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F15*F23i*c3*f2*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*dsvm*F15*F22i*c2*f2*u2*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F15*F21i*c1*f2*u1*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F15*F20i*c0*f2*u0*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F14*F25i*c5*f2*u4*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F14*F24i*c4*f2*u4**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F14*F23i*c3*f2*u3*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F14*F22i*c2*f2*u2*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F14*F21i*c1*f2*u1*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F14*F20i*c0*f2*u0*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F13*F25i*c5*f2*u3*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F13*F24i*c4*f2*u3*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F13*F23i*c3*f2*u3**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F13*F22i*c2*f2*u2*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F13*F21i*c1*f2*u1*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F13*F20i*c0*f2*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*dsvm*F12*F25i*c5*f2*u2*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F12*F24i*c4*f2*u2*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F12*F23i*c3*f2*u2*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F12*F22i*c2*f2*u2**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F12*F21i*c1*f2*u1*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F12*F20i*c0*f2*u0*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F11*F25i*c5*f2*u1*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F11*F24i*c4*f2*u1*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F11*F23i*c3*f2*u1*u3*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F11*F22i*c2*f2*u1*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F11*F21i*c1*f2*u1**2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F11*F20i*c0*f2*u0*u1*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F10*F25i*c5*f2*u0*u5*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F10*F24i*c4*f2*u0*u4*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F10*F23i*c3*f2*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2*ssp*dsvm*F10*F22i*c2*f2*u0*u2*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F10*F21i*c1*f2*u0*u1*w2
     &  - 4.D0*i_*PC2*ssp*dsvm*F10*F20i*c0*f2*u0**2*w2
     &  + 4.D0*i_*PC2*ssp*dssm*F25*F25i*c5*f2*u5**2
     &  + 4.D0*i_*PC2*ssp*dssm*F25*F24i*c4*f2*u4*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F25*F23i*c3*f2*u3*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F25*F22i*c2*f2*u2*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F25*F21i*c1*f2*u1*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F25*F20i*c0*f2*u0*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F24*F25i*c5*f2*u4*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F24*F24i*c4*f2*u4**2
     &  + 4.D0*i_*PC2*ssp*dssm*F24*F23i*c3*f2*u3*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F24*F22i*c2*f2*u2*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F24*F21i*c1*f2*u1*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F24*F20i*c0*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssp*dssm*F23*F25i*c5*f2*u3*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F23*F24i*c4*f2*u3*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F23*F23i*c3*f2*u3**2
     &  + 4.D0*i_*PC2*ssp*dssm*F23*F22i*c2*f2*u2*u3
     &  + 4.D0*i_*PC2*ssp*dssm*F23*F21i*c1*f2*u1*u3
     &  + 4.D0*i_*PC2*ssp*dssm*F23*F20i*c0*f2*u0*u3
     &  + 4.D0*i_*PC2*ssp*dssm*F22*F25i*c5*f2*u2*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F22*F24i*c4*f2*u2*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F22*F23i*c3*f2*u2*u3
     &  + 4.D0*i_*PC2*ssp*dssm*F22*F22i*c2*f2*u2**2
     &  + 4.D0*i_*PC2*ssp*dssm*F22*F21i*c1*f2*u1*u2
     &  + 4.D0*i_*PC2*ssp*dssm*F22*F20i*c0*f2*u0*u2
     &  + 4.D0*i_*PC2*ssp*dssm*F21*F25i*c5*f2*u1*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F21*F24i*c4*f2*u1*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F21*F23i*c3*f2*u1*u3
     &
      traza1 = traza1 + 4.D0*i_*PC2*ssp*dssm*F21*F22i*c2*f2*u1*u2
     &  + 4.D0*i_*PC2*ssp*dssm*F21*F21i*c1*f2*u1**2
     &  + 4.D0*i_*PC2*ssp*dssm*F21*F20i*c0*f2*u0*u1
     &  + 4.D0*i_*PC2*ssp*dssm*F20*F25i*c5*f2*u0*u5
     &  + 4.D0*i_*PC2*ssp*dssm*F20*F24i*c4*f2*u0*u4
     &  + 4.D0*i_*PC2*ssp*dssm*F20*F23i*c3*f2*u0*u3
     &  + 4.D0*i_*PC2*ssp*dssm*F20*F22i*c2*f2*u0*u2
     &  + 4.D0*i_*PC2*ssp*dssm*F20*F21i*c1*f2*u0*u1
     &  + 4.D0*i_*PC2*ssp*dssm*F20*F20i*c0*f2*u0**2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F15i*c5*f1*u5**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F15i*c5*f1*u5**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F14i*c4*f1*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F13i*c3*f1*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F13i*c3*f1*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F45*F12i*c2*f1*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F11i*c1*f1*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F10i*c0*f1*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F45*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F15i*c5*f1*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F14i*c4*f1*u4**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F14i*c4*f1*u4**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F13i*c3*f1*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F12i*c2*f1*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F12i*c2*f1*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F11i*c1*f1*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F44*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F10i*c0*f1*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F44*F10i*c0*f1*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F15i*c5*f1*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F15i*c5*f1*u3*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F14i*c4*f1*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F13i*c3*f1*u3**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F13i*c3*f1*u3**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F12i*c2*f1*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F11i*c1*f1*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F10i*c0*f1*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F43*F10i*c0*f1*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F42*F15i*c5*f1*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F14i*c4*f1*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F13i*c3*f1*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F12i*c2*f1*u2**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F12i*c2*f1*u2**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F11i*c1*f1*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F10i*c0*f1*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F42*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F15i*c5*f1*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F15i*c5*f1*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F14i*c4*f1*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F41*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F13i*c3*f1*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F13i*c3*f1*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F12i*c2*f1*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F12i*c2*f1*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F11i*c1*f1*u1**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F11i*c1*f1*u1**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F10i*c0*f1*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F41*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F15i*c5*f1*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F14i*c4*f1*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F13i*c3*f1*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F13i*c3*f1*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F40*F12i*c2*f1*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F11i*c1*f1*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F10i*c0*f1*u0**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F40*F10i*c0*f1*u0**2*w1
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F25*F25i*c5*f2*u5**2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F25*F24i*c4*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F25*F23i*c3*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F25*F22i*c2*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F25*F21i*c1*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F25*F20i*c0*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F24*F25i*c5*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F24*F24i*c4*f2*u4**2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F24*F23i*c3*f2*u3*u4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dsvp*F24*F22i*c2*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F24*F21i*c1*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F24*F20i*c0*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F23*F25i*c5*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F23*F24i*c4*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F23*F23i*c3*f2*u3**2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F23*F22i*c2*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F23*F21i*c1*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F23*F20i*c0*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F22*F25i*c5*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F22*F24i*c4*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F22*F23i*c3*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F22*F22i*c2*f2*u2**2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F22*F21i*c1*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F22*F20i*c0*f2*u0*u2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dsvp*F21*F25i*c5*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F21*F24i*c4*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F21*F23i*c3*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F21*F22i*c2*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F21*F21i*c1*f2*u1**2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F21*F20i*c0*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F20*F25i*c5*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F20*F24i*c4*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F20*F23i*c3*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F20*F22i*c2*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F20*F21i*c1*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svm*dsvp*F20*F20i*c0*f2*u0**2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F45i*c5*f4*u5**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F45i*c5*f4*u5**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F44i*c4*f4*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F15*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F15*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F44i*c4*f4*u4**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F44i*c4*f4*u4**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F43i*c3*f4*u3*u4*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F14*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F14*F40i*c0*f4*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F43i*c3*f4*u3**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F43i*c3*f4*u3**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F41i*c1*f4*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F13*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F13*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F42i*c2*f4*u2**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F42i*c2*f4*u2**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F12*F40i*c0*f4*u0*u2*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F11*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F41i*c1*f4*u1**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F41i*c1*f4*u1**2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F11*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F44i*c4*f4*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svm*dsvp*F10*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F43i*c3*f4*u0*u3*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F40i*c0*f4*u0**2*w2
     &  + 4.D0*i_*PC2*qq*svm*dsvp*F10*F40i*c0*f4*u0**2*w1
     &  - 4.D0*i_*PC2*qq*svm*dssp*F45*F25i*c5*f2*u5**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F45*F24i*c4*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F45*F23i*c3*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F45*F22i*c2*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F45*F21i*c1*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F45*F20i*c0*f2*u0*u5
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dssp*F44*F25i*c5*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F44*F24i*c4*f2*u4**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F44*F23i*c3*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F44*F22i*c2*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F44*F21i*c1*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F44*F20i*c0*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F43*F25i*c5*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F43*F24i*c4*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F43*F23i*c3*f2*u3**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F43*F22i*c2*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F43*F21i*c1*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F43*F20i*c0*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F42*F25i*c5*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F42*F24i*c4*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F42*F23i*c3*f2*u2*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dssp*F42*F22i*c2*f2*u2**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F42*F21i*c1*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F42*F20i*c0*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F41*F25i*c5*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F41*F24i*c4*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F41*F23i*c3*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F41*F22i*c2*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F41*F21i*c1*f2*u1**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F41*F20i*c0*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svm*dssp*F40*F25i*c5*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F40*F24i*c4*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F40*F23i*c3*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F40*F22i*c2*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F40*F21i*c1*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svm*dssp*F40*F20i*c0*f2*u0**2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dssp*F25*F45i*c5*f4*u5**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F25*F44i*c4*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F25*F43i*c3*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F25*F42i*c2*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F25*F41i*c1*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F25*F40i*c0*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F24*F45i*c5*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F24*F44i*c4*f4*u4**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F24*F43i*c3*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F24*F42i*c2*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F24*F41i*c1*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F24*F40i*c0*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F23*F45i*c5*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F23*F44i*c4*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F23*F43i*c3*f4*u3**2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dssp*F23*F42i*c2*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F23*F41i*c1*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F23*F40i*c0*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F22*F45i*c5*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F22*F44i*c4*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F22*F43i*c3*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F22*F42i*c2*f4*u2**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F22*F41i*c1*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F22*F40i*c0*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F21*F45i*c5*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F21*F44i*c4*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F21*F43i*c3*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F21*F42i*c2*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F21*F41i*c1*f4*u1**2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F21*F40i*c0*f4*u0*u1
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svm*dssp*F20*F45i*c5*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*svm*dssp*F20*F44i*c4*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*svm*dssp*F20*F43i*c3*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*svm*dssp*F20*F42i*c2*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*svm*dssp*F20*F41i*c1*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*svm*dssp*F20*F40i*c0*f4*u0**2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F15i*c5*f1*u5**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F15i*c5*f1*u5**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F14i*c4*f1*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F13i*c3*f1*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F13i*c3*f1*u3*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F12i*c2*f1*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F11i*c1*f1*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F45*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F10i*c0*f1*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F45*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F15i*c5*f1*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F14i*c4*f1*u4**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F14i*c4*f1*u4**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F13i*c3*f1*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F12i*c2*f1*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F12i*c2*f1*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F11i*c1*f1*u1*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F10i*c0*f1*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F44*F10i*c0*f1*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F43*F15i*c5*f1*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F15i*c5*f1*u3*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F14i*c4*f1*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F13i*c3*f1*u3**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F13i*c3*f1*u3**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F12i*c2*f1*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F11i*c1*f1*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F10i*c0*f1*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F43*F10i*c0*f1*u0*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F15i*c5*f1*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F14i*c4*f1*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F42*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F13i*c3*f1*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F12i*c2*f1*u2**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F12i*c2*f1*u2**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F11i*c1*f1*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F10i*c0*f1*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F42*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F15i*c5*f1*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F15i*c5*f1*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F14i*c4*f1*u1*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F13i*c3*f1*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F13i*c3*f1*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F41*F12i*c2*f1*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F12i*c2*f1*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F11i*c1*f1*u1**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F11i*c1*f1*u1**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F10i*c0*f1*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F41*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F15i*c5*f1*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F14i*c4*f1*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F13i*c3*f1*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F13i*c3*f1*u0*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F12i*c2*f1*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F11i*c1*f1*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F40*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F10i*c0*f1*u0**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F40*F10i*c0*f1*u0**2*w1
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F25*F25i*c5*f2*u5**2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F25*F24i*c4*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F25*F23i*c3*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F25*F22i*c2*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F25*F21i*c1*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F25*F20i*c0*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F24*F25i*c5*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F24*F24i*c4*f2*u4**2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F24*F23i*c3*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F24*F22i*c2*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F24*F21i*c1*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F24*F20i*c0*f2*u0*u4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dsvm*F23*F25i*c5*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F23*F24i*c4*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F23*F23i*c3*f2*u3**2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F23*F22i*c2*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F23*F21i*c1*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F23*F20i*c0*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F22*F25i*c5*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F22*F24i*c4*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F22*F23i*c3*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F22*F22i*c2*f2*u2**2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F22*F21i*c1*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F22*F20i*c0*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F21*F25i*c5*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F21*F24i*c4*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F21*F23i*c3*f2*u1*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dsvm*F21*F22i*c2*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F21*F21i*c1*f2*u1**2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F21*F20i*c0*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F20*F25i*c5*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F20*F24i*c4*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F20*F23i*c3*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F20*F22i*c2*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F20*F21i*c1*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svp*dsvm*F20*F20i*c0*f2*u0**2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F45i*c5*f4*u5**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F45i*c5*f4*u5**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F43i*c3*f4*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F15*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F15*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F44i*c4*f4*u4**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F44i*c4*f4*u4**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F41i*c1*f4*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F14*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F14*F40i*c0*f4*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F43i*c3*f4*u3**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F43i*c3*f4*u3**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F13*F40i*c0*f4*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F12*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F42i*c2*f4*u2**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F42i*c2*f4*u2**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F12*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F44i*c4*f4*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F11*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F41i*c1*f4*u1**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F41i*c1*f4*u1**2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F11*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F43i*c3*f4*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq*svp*dsvm*F10*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F40i*c0*f4*u0**2*w2
     &  + 4.D0*i_*PC2*qq*svp*dsvm*F10*F40i*c0*f4*u0**2*w1
     &  - 4.D0*i_*PC2*qq*svp*dssm*F45*F25i*c5*f2*u5**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F45*F24i*c4*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F45*F23i*c3*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F45*F22i*c2*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F45*F21i*c1*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F45*F20i*c0*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F44*F25i*c5*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F44*F24i*c4*f2*u4**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F44*F23i*c3*f2*u3*u4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dssm*F44*F22i*c2*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F44*F21i*c1*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F44*F20i*c0*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F43*F25i*c5*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F43*F24i*c4*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F43*F23i*c3*f2*u3**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F43*F22i*c2*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F43*F21i*c1*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F43*F20i*c0*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F42*F25i*c5*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F42*F24i*c4*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F42*F23i*c3*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F42*F22i*c2*f2*u2**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F42*F21i*c1*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F42*F20i*c0*f2*u0*u2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dssm*F41*F25i*c5*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F41*F24i*c4*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F41*F23i*c3*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F41*F22i*c2*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F41*F21i*c1*f2*u1**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F41*F20i*c0*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svp*dssm*F40*F25i*c5*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F40*F24i*c4*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F40*F23i*c3*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F40*F22i*c2*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F40*F21i*c1*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*svp*dssm*F40*F20i*c0*f2*u0**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F25*F45i*c5*f4*u5**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F25*F44i*c4*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F25*F43i*c3*f4*u3*u5
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dssm*F25*F42i*c2*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F25*F41i*c1*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F25*F40i*c0*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F24*F45i*c5*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F24*F44i*c4*f4*u4**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F24*F43i*c3*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F24*F42i*c2*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F24*F41i*c1*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F24*F40i*c0*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F23*F45i*c5*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F23*F44i*c4*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F23*F43i*c3*f4*u3**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F23*F42i*c2*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F23*F41i*c1*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F23*F40i*c0*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dssm*F22*F45i*c5*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F22*F44i*c4*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F22*F43i*c3*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F22*F42i*c2*f4*u2**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F22*F41i*c1*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F22*F40i*c0*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F21*F45i*c5*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F21*F44i*c4*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F21*F43i*c3*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*svp*dssm*F21*F42i*c2*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F21*F41i*c1*f4*u1**2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F21*F40i*c0*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*svp*dssm*F20*F45i*c5*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*svp*dssm*F20*F44i*c4*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*svp*dssm*F20*F43i*c3*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*svp*dssm*F20*F42i*c2*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*svp*dssm*F20*F41i*c1*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*svp*dssm*F20*F40i*c0*f4*u0**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F45*F25i*c5*f2*u5**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F45*F24i*c4*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F45*F23i*c3*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F45*F22i*c2*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F45*F21i*c1*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F45*F20i*c0*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F44*F25i*c5*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F44*F24i*c4*f2*u4**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F44*F23i*c3*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F44*F22i*c2*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F44*F21i*c1*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F44*F20i*c0*f2*u0*u4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dsvp*F43*F25i*c5*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F43*F24i*c4*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F43*F23i*c3*f2*u3**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F43*F22i*c2*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F43*F21i*c1*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F43*F20i*c0*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F42*F25i*c5*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F42*F24i*c4*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F42*F23i*c3*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F42*F22i*c2*f2*u2**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F42*F21i*c1*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F42*F20i*c0*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F41*F25i*c5*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F41*F24i*c4*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F41*F23i*c3*f2*u1*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dsvp*F41*F22i*c2*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F41*F21i*c1*f2*u1**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F41*F20i*c0*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F40*F25i*c5*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F40*F24i*c4*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F40*F23i*c3*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F40*F22i*c2*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F40*F21i*c1*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F40*F20i*c0*f2*u0**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F25*F45i*c5*f4*u5**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F25*F44i*c4*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F25*F43i*c3*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F25*F42i*c2*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F25*F41i*c1*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F25*F40i*c0*f4*u0*u5
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dsvp*F24*F45i*c5*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F24*F44i*c4*f4*u4**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F24*F43i*c3*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F24*F42i*c2*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F24*F41i*c1*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F24*F40i*c0*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F23*F45i*c5*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F23*F44i*c4*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F23*F43i*c3*f4*u3**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F23*F42i*c2*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F23*F41i*c1*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F23*F40i*c0*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F22*F45i*c5*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F22*F44i*c4*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F22*F43i*c3*f4*u2*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dsvp*F22*F42i*c2*f4*u2**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F22*F41i*c1*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F22*F40i*c0*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F21*F45i*c5*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F21*F44i*c4*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F21*F43i*c3*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F21*F42i*c2*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F21*F41i*c1*f4*u1**2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F21*F40i*c0*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F20*F45i*c5*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F20*F44i*c4*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F20*F43i*c3*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F20*F42i*c2*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F20*F41i*c1*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssm*dsvp*F20*F40i*c0*f4*u0**2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dssp*F45*F45i*c5*f4*u5**2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F45*F44i*c4*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F45*F43i*c3*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F45*F42i*c2*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F45*F41i*c1*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F45*F40i*c0*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F44*F45i*c5*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F44*F44i*c4*f4*u4**2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F44*F43i*c3*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F44*F42i*c2*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F44*F41i*c1*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F44*F40i*c0*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F43*F45i*c5*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F43*F44i*c4*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F43*F43i*c3*f4*u3**2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dssp*F43*F42i*c2*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F43*F41i*c1*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F43*F40i*c0*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F42*F45i*c5*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F42*F44i*c4*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F42*F43i*c3*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F42*F42i*c2*f4*u2**2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F42*F41i*c1*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F42*F40i*c0*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F41*F45i*c5*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F41*F44i*c4*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F41*F43i*c3*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F41*F42i*c2*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F41*F41i*c1*f4*u1**2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F41*F40i*c0*f4*u0*u1
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssm*dssp*F40*F45i*c5*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F40*F44i*c4*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F40*F43i*c3*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F40*F42i*c2*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F40*F41i*c1*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssm*dssp*F40*F40i*c0*f4*u0**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F45*F25i*c5*f2*u5**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F45*F24i*c4*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F45*F23i*c3*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F45*F22i*c2*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F45*F21i*c1*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F45*F20i*c0*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F44*F25i*c5*f2*u4*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F44*F24i*c4*f2*u4**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F44*F23i*c3*f2*u3*u4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dsvm*F44*F22i*c2*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F44*F21i*c1*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F44*F20i*c0*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F43*F25i*c5*f2*u3*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F43*F24i*c4*f2*u3*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F43*F23i*c3*f2*u3**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F43*F22i*c2*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F43*F21i*c1*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F43*F20i*c0*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F42*F25i*c5*f2*u2*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F42*F24i*c4*f2*u2*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F42*F23i*c3*f2*u2*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F42*F22i*c2*f2*u2**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F42*F21i*c1*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F42*F20i*c0*f2*u0*u2
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dsvm*F41*F25i*c5*f2*u1*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F41*F24i*c4*f2*u1*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F41*F23i*c3*f2*u1*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F41*F22i*c2*f2*u1*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F41*F21i*c1*f2*u1**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F41*F20i*c0*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F40*F25i*c5*f2*u0*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F40*F24i*c4*f2*u0*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F40*F23i*c3*f2*u0*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F40*F22i*c2*f2*u0*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F40*F21i*c1*f2*u0*u1
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F40*F20i*c0*f2*u0**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F25*F45i*c5*f4*u5**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F25*F44i*c4*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F25*F43i*c3*f4*u3*u5
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dsvm*F25*F42i*c2*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F25*F41i*c1*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F25*F40i*c0*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F24*F45i*c5*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F24*F44i*c4*f4*u4**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F24*F43i*c3*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F24*F42i*c2*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F24*F41i*c1*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F24*F40i*c0*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F23*F45i*c5*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F23*F44i*c4*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F23*F43i*c3*f4*u3**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F23*F42i*c2*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F23*F41i*c1*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F23*F40i*c0*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dsvm*F22*F45i*c5*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F22*F44i*c4*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F22*F43i*c3*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F22*F42i*c2*f4*u2**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F22*F41i*c1*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F22*F40i*c0*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F21*F45i*c5*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F21*F44i*c4*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F21*F43i*c3*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F21*F42i*c2*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F21*F41i*c1*f4*u1**2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F21*F40i*c0*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F20*F45i*c5*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F20*F44i*c4*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F20*F43i*c3*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dsvm*F20*F42i*c2*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F20*F41i*c1*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssp*dsvm*F20*F40i*c0*f4*u0**2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F45*F45i*c5*f4*u5**2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F45*F44i*c4*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F45*F43i*c3*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F45*F42i*c2*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F45*F41i*c1*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F45*F40i*c0*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F44*F45i*c5*f4*u4*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F44*F44i*c4*f4*u4**2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F44*F43i*c3*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F44*F42i*c2*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F44*F41i*c1*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F44*F40i*c0*f4*u0*u4
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dssm*F43*F45i*c5*f4*u3*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F43*F44i*c4*f4*u3*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F43*F43i*c3*f4*u3**2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F43*F42i*c2*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F43*F41i*c1*f4*u1*u3
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F43*F40i*c0*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F42*F45i*c5*f4*u2*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F42*F44i*c4*f4*u2*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F42*F43i*c3*f4*u2*u3
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F42*F42i*c2*f4*u2**2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F42*F41i*c1*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F42*F40i*c0*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F41*F45i*c5*f4*u1*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F41*F44i*c4*f4*u1*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F41*F43i*c3*f4*u1*u3
     &
      traza1 = traza1 - 4.D0*i_*PC2*qq*ssp*dssm*F41*F42i*c2*f4*u1*u2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F41*F41i*c1*f4*u1**2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F41*F40i*c0*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F40*F45i*c5*f4*u0*u5
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F40*F44i*c4*f4*u0*u4
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F40*F43i*c3*f4*u0*u3
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F40*F42i*c2*f4*u0*u2
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F40*F41i*c1*f4*u0*u1
     &  - 4.D0*i_*PC2*qq*ssp*dssm*F40*F40i*c0*f4*u0**2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F45*F45i*c5*f4*u5**2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F45*F44i*c4*f4*u4*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F45*F43i*c3*f4*u3*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F45*F42i*c2*f4*u2*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F45*F41i*c1*f4*u1*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F45*F40i*c0*f4*u0*u5
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svm*dsvp*F44*F45i*c5*f4*u4*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F44*F44i*c4*f4*u4**2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F44*F43i*c3*f4*u3*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F44*F42i*c2*f4*u2*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F44*F41i*c1*f4*u1*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F44*F40i*c0*f4*u0*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F43*F45i*c5*f4*u3*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F43*F44i*c4*f4*u3*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F43*F43i*c3*f4*u3**2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F43*F42i*c2*f4*u2*u3
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F43*F41i*c1*f4*u1*u3
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F43*F40i*c0*f4*u0*u3
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F42*F45i*c5*f4*u2*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F42*F44i*c4*f4*u2*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F42*F43i*c3*f4*u2*u3
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svm*dsvp*F42*F42i*c2*f4*u2**2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F42*F41i*c1*f4*u1*u2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F42*F40i*c0*f4*u0*u2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F41*F45i*c5*f4*u1*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F41*F44i*c4*f4*u1*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F41*F43i*c3*f4*u1*u3
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F41*F42i*c2*f4*u1*u2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F41*F41i*c1*f4*u1**2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F41*F40i*c0*f4*u0*u1
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F40*F45i*c5*f4*u0*u5
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F40*F44i*c4*f4*u0*u4
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F40*F43i*c3*f4*u0*u3
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F40*F42i*c2*f4*u0*u2
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F40*F41i*c1*f4*u0*u1
     &  + 4.D0*i_*PC2*qq**2*svm*dsvp*F40*F40i*c0*f4*u0**2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svp*dsvm*F45*F45i*c5*f4*u5**2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F45*F44i*c4*f4*u4*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F45*F43i*c3*f4*u3*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F45*F42i*c2*f4*u2*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F45*F41i*c1*f4*u1*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F45*F40i*c0*f4*u0*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F44*F45i*c5*f4*u4*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F44*F44i*c4*f4*u4**2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F44*F43i*c3*f4*u3*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F44*F42i*c2*f4*u2*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F44*F41i*c1*f4*u1*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F44*F40i*c0*f4*u0*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F43*F45i*c5*f4*u3*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F43*F44i*c4*f4*u3*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F43*F43i*c3*f4*u3**2
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svp*dsvm*F43*F42i*c2*f4*u2*u3
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F43*F41i*c1*f4*u1*u3
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F43*F40i*c0*f4*u0*u3
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F42*F45i*c5*f4*u2*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F42*F44i*c4*f4*u2*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F42*F43i*c3*f4*u2*u3
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F42*F42i*c2*f4*u2**2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F42*F41i*c1*f4*u1*u2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F42*F40i*c0*f4*u0*u2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F41*F45i*c5*f4*u1*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F41*F44i*c4*f4*u1*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F41*F43i*c3*f4*u1*u3
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F41*F42i*c2*f4*u1*u2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F41*F41i*c1*f4*u1**2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F41*F40i*c0*f4*u0*u1
     &
      traza1 = traza1 + 4.D0*i_*PC2*qq**2*svp*dsvm*F40*F45i*c5*f4*u0*u5
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F40*F44i*c4*f4*u0*u4
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F40*F43i*c3*f4*u0*u3
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F40*F42i*c2*f4*u0*u2
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F40*F41i*c1*f4*u0*u1
     &  + 4.D0*i_*PC2*qq**2*svp*dsvm*F40*F40i*c0*f4*u0**2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F25*F25i*c5*f2*u5**2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F25*F24i*c4*f2*u4*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F25*F23i*c3*f2*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F25*F22i*c2*f2*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F25*F21i*c1*f2*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F25*F20i*c0*f2*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F24*F25i*c5*f2*u4*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F24*F24i*c4*f2*u4**2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F24*F23i*c3*f2*u3*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svm*dsvp*F24*F22i*c2*f2*u2*u4*w1
     & *w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F24*F21i*c1*f2*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F24*F20i*c0*f2*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F23*F25i*c5*f2*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F23*F24i*c4*f2*u3*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F23*F23i*c3*f2*u3**2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F23*F22i*c2*f2*u2*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F23*F21i*c1*f2*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F23*F20i*c0*f2*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F22*F25i*c5*f2*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F22*F24i*c4*f2*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F22*F23i*c3*f2*u2*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F22*F22i*c2*f2*u2**2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F22*F21i*c1*f2*u1*u2*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svm*dsvp*F22*F20i*c0*f2*u0*u2*w1
     & *w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F21*F25i*c5*f2*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F21*F24i*c4*f2*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F21*F23i*c3*f2*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F21*F22i*c2*f2*u1*u2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F21*F21i*c1*f2*u1**2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F21*F20i*c0*f2*u0*u1*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F20*F25i*c5*f2*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F20*F24i*c4*f2*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F20*F23i*c3*f2*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F20*F22i*c2*f2*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F20*F21i*c1*f2*u0*u1*w1*w2
     &  - 4.D0*i_*PC2**2*svm*dsvp*F20*F20i*c0*f2*u0**2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F25*F25i*c5*f2*u5**2*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svp*dsvm*F25*F24i*c4*f2*u4*u5*w1
     & *w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F25*F23i*c3*f2*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F25*F22i*c2*f2*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F25*F21i*c1*f2*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F25*F20i*c0*f2*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F24*F25i*c5*f2*u4*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F24*F24i*c4*f2*u4**2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F24*F23i*c3*f2*u3*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F24*F22i*c2*f2*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F24*F21i*c1*f2*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F24*F20i*c0*f2*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F23*F25i*c5*f2*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F23*F24i*c4*f2*u3*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F23*F23i*c3*f2*u3**2*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svp*dsvm*F23*F22i*c2*f2*u2*u3*w1
     & *w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F23*F21i*c1*f2*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F23*F20i*c0*f2*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F22*F25i*c5*f2*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F22*F24i*c4*f2*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F22*F23i*c3*f2*u2*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F22*F22i*c2*f2*u2**2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F22*F21i*c1*f2*u1*u2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F22*F20i*c0*f2*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F21*F25i*c5*f2*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F21*F24i*c4*f2*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F21*F23i*c3*f2*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F21*F22i*c2*f2*u1*u2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F21*F21i*c1*f2*u1**2*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*svp*dsvm*F21*F20i*c0*f2*u0*u1*w1
     & *w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F20*F25i*c5*f2*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F20*F24i*c4*f2*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F20*F23i*c3*f2*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F20*F22i*c2*f2*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F20*F21i*c1*f2*u0*u1*w1*w2
     &  - 4.D0*i_*PC2**2*svp*dsvm*F20*F20i*c0*f2*u0**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F45*F45i*c5*f4*u5**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F45*F44i*c4*f4*u4*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F45*F43i*c3*f4*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F45*F42i*c2*f4*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F45*F41i*c1*f4*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F45*F40i*c0*f4*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F44*F45i*c5*f4*u4*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svm*dsvp*F44*F44i*c4*f4*u4**2
     & *w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F44*F43i*c3*f4*u3*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F44*F42i*c2*f4*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F44*F41i*c1*f4*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F44*F40i*c0*f4*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F43*F45i*c5*f4*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F43*F44i*c4*f4*u3*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F43*F43i*c3*f4*u3**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F43*F42i*c2*f4*u2*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F43*F41i*c1*f4*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F43*F40i*c0*f4*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F42*F45i*c5*f4*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F42*F44i*c4*f4*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F42*F43i*c3*f4*u2*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svm*dsvp*F42*F42i*c2*f4*u2**2
     & *w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F42*F41i*c1*f4*u1*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F42*F40i*c0*f4*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F41*F45i*c5*f4*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F41*F44i*c4*f4*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F41*F43i*c3*f4*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F41*F42i*c2*f4*u1*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F41*F41i*c1*f4*u1**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F41*F40i*c0*f4*u0*u1*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F40*F45i*c5*f4*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F40*F44i*c4*f4*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F40*F43i*c3*f4*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F40*F42i*c2*f4*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svm*dsvp*F40*F41i*c1*f4*u0*u1*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svm*dsvp*F40*F40i*c0*f4*u0**2
     & *w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F45*F45i*c5*f4*u5**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F45*F44i*c4*f4*u4*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F45*F43i*c3*f4*u3*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F45*F42i*c2*f4*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F45*F41i*c1*f4*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F45*F40i*c0*f4*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F44*F45i*c5*f4*u4*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F44*F44i*c4*f4*u4**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F44*F43i*c3*f4*u3*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F44*F42i*c2*f4*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F44*F41i*c1*f4*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F44*F40i*c0*f4*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F43*F45i*c5*f4*u3*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svp*dsvm*F43*F44i*c4*f4*u3*u4
     & *w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F43*F43i*c3*f4*u3**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F43*F42i*c2*f4*u2*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F43*F41i*c1*f4*u1*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F43*F40i*c0*f4*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F42*F45i*c5*f4*u2*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F42*F44i*c4*f4*u2*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F42*F43i*c3*f4*u2*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F42*F42i*c2*f4*u2**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F42*F41i*c1*f4*u1*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F42*F40i*c0*f4*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F41*F45i*c5*f4*u1*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F41*F44i*c4*f4*u1*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F41*F43i*c3*f4*u1*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*i_*PC2**2*qq*svp*dsvm*F41*F42i*c2*f4*u1*u2
     & *w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F41*F41i*c1*f4*u1**2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F41*F40i*c0*f4*u0*u1*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F40*F45i*c5*f4*u0*u5*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F40*F44i*c4*f4*u0*u4*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F40*F43i*c3*f4*u0*u3*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F40*F42i*c2*f4*u0*u2*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F40*F41i*c1*f4*u0*u1*w1*w2
     &  - 4.D0*i_*PC2**2*qq*svp*dsvm*F40*F40i*c0*f4*u0**2*w1*w2
     &  + 8.D0*PC_q*PC_uh*svp*svm*F25*F25i*c5*f2*u5**2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F25*F25i*c5*f2*u5**2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F25*F24i*c4*f2*u4*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F25*F24i*c4*f2*u4*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F25*F23i*c3*f2*u3*u5*w2
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*svp*svm*F25*F23i*c3*f2*u3*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F25*F22i*c2*f2*u2*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F25*F22i*c2*f2*u2*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F25*F21i*c1*f2*u1*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F25*F21i*c1*f2*u1*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F25*F20i*c0*f2*u0*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F25*F20i*c0*f2*u0*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F24*F25i*c5*f2*u4*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F24*F25i*c5*f2*u4*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F24*F24i*c4*f2*u4**2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F24*F24i*c4*f2*u4**2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F24*F23i*c3*f2*u3*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F24*F23i*c3*f2*u3*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F24*F22i*c2*f2*u2*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F24*F22i*c2*f2*u2*u4*w1
     &
      traza1 = traza1 + 8.D0*PC_q*PC_uh*svp*svm*F24*F21i*c1*f2*u1*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F24*F21i*c1*f2*u1*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F24*F20i*c0*f2*u0*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F24*F20i*c0*f2*u0*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F23*F25i*c5*f2*u3*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F23*F25i*c5*f2*u3*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F23*F24i*c4*f2*u3*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F23*F24i*c4*f2*u3*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F23*F23i*c3*f2*u3**2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F23*F23i*c3*f2*u3**2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F23*F22i*c2*f2*u2*u3*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F23*F22i*c2*f2*u2*u3*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F23*F21i*c1*f2*u1*u3*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F23*F21i*c1*f2*u1*u3*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F23*F20i*c0*f2*u0*u3*w2
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*svp*svm*F23*F20i*c0*f2*u0*u3*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F22*F25i*c5*f2*u2*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F22*F25i*c5*f2*u2*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F22*F24i*c4*f2*u2*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F22*F24i*c4*f2*u2*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F22*F23i*c3*f2*u2*u3*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F22*F23i*c3*f2*u2*u3*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F22*F22i*c2*f2*u2**2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F22*F22i*c2*f2*u2**2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F22*F21i*c1*f2*u1*u2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F22*F21i*c1*f2*u1*u2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F22*F20i*c0*f2*u0*u2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F22*F20i*c0*f2*u0*u2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F21*F25i*c5*f2*u1*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F21*F25i*c5*f2*u1*u5*w1
     &
      traza1 = traza1 + 8.D0*PC_q*PC_uh*svp*svm*F21*F24i*c4*f2*u1*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F21*F24i*c4*f2*u1*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F21*F23i*c3*f2*u1*u3*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F21*F23i*c3*f2*u1*u3*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F21*F22i*c2*f2*u1*u2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F21*F22i*c2*f2*u1*u2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F21*F21i*c1*f2*u1**2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F21*F21i*c1*f2*u1**2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F21*F20i*c0*f2*u0*u1*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F21*F20i*c0*f2*u0*u1*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F20*F25i*c5*f2*u0*u5*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F20*F25i*c5*f2*u0*u5*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F20*F24i*c4*f2*u0*u4*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F20*F24i*c4*f2*u0*u4*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F20*F23i*c3*f2*u0*u3*w2
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*svp*svm*F20*F23i*c3*f2*u0*u3*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F20*F22i*c2*f2*u0*u2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F20*F22i*c2*f2*u0*u2*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F20*F21i*c1*f2*u0*u1*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F20*F21i*c1*f2*u0*u1*w1
     &  + 8.D0*PC_q*PC_uh*svp*svm*F20*F20i*c0*f2*u0**2*w2
     &  - 8.D0*PC_q*PC_uh*svp*svm*F20*F20i*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F45*F25i*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F45*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F45*F23i*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F45*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F45*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F45*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F44*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F44*F24i*c4*f2*u4**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*ssm*svp*F44*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F44*F22i*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F44*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F44*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F43*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F43*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F43*F23i*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F43*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F43*F21i*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F43*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F42*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F42*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F42*F23i*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F42*F22i*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F42*F21i*c1*f2*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*ssm*svp*F42*F20i*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F41*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F41*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F41*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F41*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F41*F21i*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F41*F20i*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F40*F25i*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F40*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F40*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F40*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F40*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F40*F20i*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F25*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F25*F44i*c4*f4*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*ssm*svp*F25*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F25*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F25*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F25*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F24*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F24*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F24*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F24*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F24*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F24*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F23*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F23*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F23*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F23*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F23*F41i*c1*f4*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*ssm*svp*F23*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F22*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F22*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F22*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F22*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F22*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F22*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F21*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F21*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F21*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F21*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F21*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F21*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F20*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F20*F44i*c4*f4*u0*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*ssm*svp*F20*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F20*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F20*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*ssm*svp*F20*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F45*F25i*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F45*F24i*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F45*F23i*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F45*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F45*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F45*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F44*F25i*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F44*F24i*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F44*F23i*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F44*F22i*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F44*F21i*c1*f2*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*ssp*svm*F44*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F43*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F43*F24i*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F43*F23i*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F43*F22i*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F43*F21i*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F43*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F42*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F42*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F42*F23i*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F42*F22i*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F42*F21i*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F42*F20i*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F41*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F41*F24i*c4*f2*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*ssp*svm*F41*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F41*F22i*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F41*F21i*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F41*F20i*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F40*F25i*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F40*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F40*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F40*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F40*F21i*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F40*F20i*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F25*F45i*c5*f4*u5**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F25*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F25*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F25*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F25*F41i*c1*f4*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*ssp*svm*F25*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F24*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F24*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F24*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F24*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F24*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F24*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F23*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F23*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F23*F43i*c3*f4*u3**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F23*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F23*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F23*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F22*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F22*F44i*c4*f4*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*ssp*svm*F22*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F22*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F22*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F22*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F21*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F21*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F21*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F21*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F21*F41i*c1*f4*u1**2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F21*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F20*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F20*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F20*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F20*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*ssp*svm*F20*F41i*c1*f4*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*ssp*svm*F20*F40i*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F25i*c5*f2*u5**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F25i*c5*f2*u5**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F24i*c4*f2*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F24i*c4*f2*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F23i*c3*f2*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F23i*c3*f2*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F22i*c2*f2*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F22i*c2*f2*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F21i*c1*f2*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F21i*c1*f2*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F20i*c0*f2*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F35*F20i*c0*f2*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F25i*c5*f2*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F25i*c5*f2*u4*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F24i*c4*f2*u4**2
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F24i*c4*f2*u4**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F23i*c3*f2*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F23i*c3*f2*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F22i*c2*f2*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F22i*c2*f2*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F21i*c1*f2*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F21i*c1*f2*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F20i*c0*f2*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F34*F20i*c0*f2*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F25i*c5*f2*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F25i*c5*f2*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F24i*c4*f2*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F24i*c4*f2*u3*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F23i*c3*f2*u3**2
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F23i*c3*f2*u3**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F22i*c2*f2*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F22i*c2*f2*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F21i*c1*f2*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F21i*c1*f2*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F20i*c0*f2*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F33*F20i*c0*f2*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F25i*c5*f2*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F25i*c5*f2*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F24i*c4*f2*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F24i*c4*f2*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F23i*c3*f2*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F23i*c3*f2*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F22i*c2*f2*u2**2
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F22i*c2*f2*u2**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F21i*c1*f2*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F21i*c1*f2*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F20i*c0*f2*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F32*F20i*c0*f2*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F25i*c5*f2*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F25i*c5*f2*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F24i*c4*f2*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F24i*c4*f2*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F23i*c3*f2*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F23i*c3*f2*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F22i*c2*f2*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F22i*c2*f2*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F21i*c1*f2*u1**2
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F21i*c1*f2*u1**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F20i*c0*f2*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F31*F20i*c0*f2*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F25i*c5*f2*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F25i*c5*f2*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F24i*c4*f2*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F24i*c4*f2*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F23i*c3*f2*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F23i*c3*f2*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F22i*c2*f2*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F22i*c2*f2*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F21i*c1*f2*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F21i*c1*f2*u0*u1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F20i*c0*f2*u0**2
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F30*F20i*c0*f2*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F35i*c5*f3*u5**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F35i*c5*f3*u5**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F34i*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F34i*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F33i*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F32i*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F31i*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F30i*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F25*F30i*c0*f3*u0*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F35i*c5*f3*u4*u5
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F35i*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F34i*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F33i*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F33i*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F32i*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F32i*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F31i*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F31i*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F30i*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F24*F30i*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F35i*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F35i*c5*f3*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F34i*c4*f3*u3*u4
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F34i*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F33i*c3*f3*u3**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F33i*c3*f3*u3**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F32i*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F32i*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F31i*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F23*F30i*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F35i*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F34i*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F34i*c4*f3*u2*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F33i*c3*f3*u2*u3
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F33i*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F32i*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F31i*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F31i*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F30i*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F22*F30i*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F35i*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F35i*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F34i*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F34i*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F33i*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F33i*c3*f3*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F32i*c2*f3*u1*u2
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F32i*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F31i*c1*f3*u1**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F31i*c1*f3*u1**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F30i*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F21*F30i*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F35i*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F34i*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F33i*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F32i*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F32i*c2*f3*u0*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F31i*c1*f3*u0*u1
     & *w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F31i*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F30i*c0*f3*u0**2*w2
     &  - 4.D0*PC_q*PC_uh*qq*svp*svm*F20*F30i*c0*f3*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F45*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F45*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F45*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F45*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F45*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F45*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F44*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F44*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F44*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F44*F32i*c2*f3*u2*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*qq*ssm*svp*F44*F31i*c1*f3*u1*u4
     & *w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F44*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F43*F35i*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F43*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F43*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F43*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F43*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F43*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F42*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F42*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F42*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F42*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F42*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F42*F30i*c0*f3*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*qq*ssm*svp*F41*F35i*c5*f3*u1*u5
     & *w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F41*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F41*F33i*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F41*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F41*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F41*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F40*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F40*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F40*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F40*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F40*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F40*F30i*c0*f3*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F35*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F35*F44i*c4*f4*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*qq*ssm*svp*F35*F43i*c3*f4*u3*u5
     & *w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F35*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F35*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F35*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F34*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F34*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F34*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F34*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F34*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F34*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F33*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F33*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F33*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F33*F42i*c2*f4*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*qq*ssm*svp*F33*F41i*c1*f4*u1*u3
     & *w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F33*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F32*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F32*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F32*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F32*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F32*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F32*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F31*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F31*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F31*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F31*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F31*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F31*F40i*c0*f4*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*qq*ssm*svp*F30*F45i*c5*f4*u0*u5
     & *w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F30*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F30*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F30*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F30*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*qq*ssm*svp*F30*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F45*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F45*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F45*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F45*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F45*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F45*F30i*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F44*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F44*F34i*c4*f3*u4**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*ssp*svm*F44*F33i*c3*f3*u3*u4
     & *w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F44*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F44*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F44*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F43*F35i*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F43*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F43*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F43*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F43*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F43*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F42*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F42*F34i*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F42*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F42*F32i*c2*f3*u2**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*ssp*svm*F42*F31i*c1*f3*u1*u2
     & *w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F42*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F41*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F41*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F41*F33i*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F41*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F41*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F41*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F40*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F40*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F40*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F40*F32i*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F40*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F40*F30i*c0*f3*u0**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*ssp*svm*F35*F45i*c5*f4*u5**2
     & *w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F35*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F35*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F35*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F35*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F35*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F34*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F34*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F34*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F34*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F34*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F34*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F33*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F33*F44i*c4*f4*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*ssp*svm*F33*F43i*c3*f4*u3**2
     & *w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F33*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F33*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F33*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F32*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F32*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F32*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F32*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F32*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F32*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F31*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F31*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F31*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F31*F42i*c2*f4*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*qq*ssp*svm*F31*F41i*c1*f4*u1**2
     & *w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F31*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F30*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F30*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F30*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F30*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F30*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*qq*ssp*svm*F30*F40i*c0*f4*u0**2*w2
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F25r*c5*f2*u5**2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F25r*c5*f2*u5**2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F24r*c4*f2*u4*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F24r*c4*f2*u4*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F23r*c3*f2*u3*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F23r*c3*f2*u3*u5*w1
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F22r*c2*f2*u2*u5
     & *w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F22r*c2*f2*u2*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F21r*c1*f2*u1*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F21r*c1*f2*u1*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F20r*c0*f2*u0*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F25*F20r*c0*f2*u0*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F25r*c5*f2*u4*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F25r*c5*f2*u4*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F24r*c4*f2*u4**2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F24r*c4*f2*u4**2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F23r*c3*f2*u3*u4*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F23r*c3*f2*u3*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F22r*c2*f2*u2*u4*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F22r*c2*f2*u2*u4*w1
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F21r*c1*f2*u1*u4
     & *w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F21r*c1*f2*u1*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F20r*c0*f2*u0*u4*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F24*F20r*c0*f2*u0*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F25r*c5*f2*u3*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F25r*c5*f2*u3*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F24r*c4*f2*u3*u4*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F24r*c4*f2*u3*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F23r*c3*f2*u3**2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F23r*c3*f2*u3**2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F22r*c2*f2*u2*u3*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F22r*c2*f2*u2*u3*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F21r*c1*f2*u1*u3*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F21r*c1*f2*u1*u3*w1
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F20r*c0*f2*u0*u3
     & *w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F23*F20r*c0*f2*u0*u3*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F25r*c5*f2*u2*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F25r*c5*f2*u2*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F24r*c4*f2*u2*u4*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F24r*c4*f2*u2*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F23r*c3*f2*u2*u3*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F23r*c3*f2*u2*u3*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F22r*c2*f2*u2**2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F22r*c2*f2*u2**2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F21r*c1*f2*u1*u2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F21r*c1*f2*u1*u2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F20r*c0*f2*u0*u2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F22*F20r*c0*f2*u0*u2*w1
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F25r*c5*f2*u1*u5
     & *w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F25r*c5*f2*u1*u5*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F24r*c4*f2*u1*u4*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F24r*c4*f2*u1*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F23r*c3*f2*u1*u3*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F23r*c3*f2*u1*u3*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F22r*c2*f2*u1*u2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F22r*c2*f2*u1*u2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F21r*c1*f2*u1**2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F21r*c1*f2*u1**2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F20r*c0*f2*u0*u1*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F21*F20r*c0*f2*u0*u1*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F25r*c5*f2*u0*u5*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F25r*c5*f2*u0*u5*w1
     &
      traza1 = traza1 - 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F24r*c4*f2*u0*u4
     & *w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F24r*c4*f2*u0*u4*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F23r*c3*f2*u0*u3*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F23r*c3*f2*u0*u3*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F22r*c2*f2*u0*u2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F22r*c2*f2*u0*u2*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F21r*c1*f2*u0*u1*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F21r*c1*f2*u0*u1*w1
     &  - 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F20r*c0*f2*u0**2*w2
     &  + 8.D0*PC_q*PC_uh*i_*svp*svm*F20*F20r*c0*f2*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F45*F25r*c5*f2*u5**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F45*F24r*c4*f2*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F45*F23r*c3*f2*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F45*F22r*c2*f2*u2*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*ssm*svp*F45*F21r*c1*f2*u1*u5
     & *w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F45*F20r*c0*f2*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F44*F25r*c5*f2*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F44*F24r*c4*f2*u4**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F44*F23r*c3*f2*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F44*F22r*c2*f2*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F44*F21r*c1*f2*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F44*F20r*c0*f2*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F43*F25r*c5*f2*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F43*F24r*c4*f2*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F43*F23r*c3*f2*u3**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F43*F22r*c2*f2*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F43*F21r*c1*f2*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F43*F20r*c0*f2*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*ssm*svp*F42*F25r*c5*f2*u2*u5
     & *w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F42*F24r*c4*f2*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F42*F23r*c3*f2*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F42*F22r*c2*f2*u2**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F42*F21r*c1*f2*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F42*F20r*c0*f2*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F41*F25r*c5*f2*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F41*F24r*c4*f2*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F41*F23r*c3*f2*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F41*F22r*c2*f2*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F41*F21r*c1*f2*u1**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F41*F20r*c0*f2*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F40*F25r*c5*f2*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F40*F24r*c4*f2*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*ssm*svp*F40*F23r*c3*f2*u0*u3
     & *w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F40*F22r*c2*f2*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F40*F21r*c1*f2*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F40*F20r*c0*f2*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F25*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F25*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F25*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F25*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F25*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F25*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F24*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F24*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F24*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F24*F42r*c2*f4*u2*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*ssm*svp*F24*F41r*c1*f4*u1*u4
     & *w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F24*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F23*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F23*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F23*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F23*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F23*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F23*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F22*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F22*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F22*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F22*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F22*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F22*F40r*c0*f4*u0*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*ssm*svp*F21*F45r*c5*f4*u1*u5
     & *w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F21*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F21*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F21*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F21*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F21*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F20*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F20*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F20*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F20*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F20*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*ssm*svp*F20*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F45*F25r*c5*f2*u5**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F45*F24r*c4*f2*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*ssp*svm*F45*F23r*c3*f2*u3*u5
     & *w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F45*F22r*c2*f2*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F45*F21r*c1*f2*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F45*F20r*c0*f2*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F44*F25r*c5*f2*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F44*F24r*c4*f2*u4**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F44*F23r*c3*f2*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F44*F22r*c2*f2*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F44*F21r*c1*f2*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F44*F20r*c0*f2*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F43*F25r*c5*f2*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F43*F24r*c4*f2*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F43*F23r*c3*f2*u3**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F43*F22r*c2*f2*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*ssp*svm*F43*F21r*c1*f2*u1*u3
     & *w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F43*F20r*c0*f2*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F42*F25r*c5*f2*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F42*F24r*c4*f2*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F42*F23r*c3*f2*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F42*F22r*c2*f2*u2**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F42*F21r*c1*f2*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F42*F20r*c0*f2*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F41*F25r*c5*f2*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F41*F24r*c4*f2*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F41*F23r*c3*f2*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F41*F22r*c2*f2*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F41*F21r*c1*f2*u1**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F41*F20r*c0*f2*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*ssp*svm*F40*F25r*c5*f2*u0*u5
     & *w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F40*F24r*c4*f2*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F40*F23r*c3*f2*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F40*F22r*c2*f2*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F40*F21r*c1*f2*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F40*F20r*c0*f2*u0**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F25*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F25*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F25*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F25*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F25*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F25*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F24*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F24*F44r*c4*f4*u4**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*ssp*svm*F24*F43r*c3*f4*u3*u4
     & *w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F24*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F24*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F24*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F23*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F23*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F23*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F23*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F23*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F23*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F22*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F22*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F22*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F22*F42r*c2*f4*u2**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*ssp*svm*F22*F41r*c1*f4*u1*u2
     & *w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F22*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F21*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F21*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F21*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F21*F42r*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F21*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F21*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F20*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F20*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F20*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F20*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F20*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*ssp*svm*F20*F40r*c0*f4*u0**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F25r*c5*f2*
     & u5**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F24r*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F23r*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F22r*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F21r*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F20r*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F35*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F25r*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F25r*c5*f2*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F24r*c4*f2*
     & u4**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F23r*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F22r*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F21r*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F20r*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F34*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F25r*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F24r*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F24r*c4*f2*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F23r*c3*f2*
     & u3**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F22r*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F21r*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F20r*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F33*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F25r*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F24r*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F23r*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F23r*c3*f2*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F22r*c2*f2*
     & u2**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F22r*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F21r*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F20r*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F32*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F25r*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F24r*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F23r*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F22r*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F22r*c2*f2*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F21r*c1*f2*
     & u1**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F20r*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F31*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F25r*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F24r*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F23r*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F22r*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F21r*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F21r*c1*f2*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F20r*c0*f2*
     & u0**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F30*F20r*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F25*F30r*c0*f3*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F35r*c5*f3*u4
     & *u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F24*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F35r*c5*f3*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F34r*c4*f3*u3
     & *u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F23*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F34r*c4*f3*u2*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F33r*c3*f3*u2
     & *u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F22*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F33r*c3*f3*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F32r*c2*f3*u1
     & *u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F21*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F32r*c2*f3*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F31r*c1*f3*u0
     & *u1*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q*PC_uh*i_*qq*svp*svm*F20*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F45*F35r*c5*f3*u5**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F45*F34r*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F45*F33r*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F45*F32r*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F45*F31r*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F45*F30r*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F44*F35r*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F44*F34r*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F44*F33r*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F44*F32r*c2*f3*u2*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F44*F31r*c1*f3*u1
     & *u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F44*F30r*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F43*F35r*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F43*F34r*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F43*F33r*c3*f3*u3**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F43*F32r*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F43*F31r*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F43*F30r*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F42*F35r*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F42*F34r*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F42*F33r*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F42*F32r*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F42*F31r*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F42*F30r*c0*f3*u0*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F41*F35r*c5*f3*u1
     & *u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F41*F34r*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F41*F33r*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F41*F32r*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F41*F31r*c1*f3*u1**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F41*F30r*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F40*F35r*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F40*F34r*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F40*F33r*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F40*F32r*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F40*F31r*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F40*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F35*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F35*F44r*c4*f4*u4*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F35*F43r*c3*f4*u3
     & *u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F35*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F35*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F35*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F34*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F34*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F34*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F34*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F34*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F34*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F33*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F33*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F33*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F33*F42r*c2*f4*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F33*F41r*c1*f4*u1
     & *u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F33*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F32*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F32*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F32*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F32*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F32*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F32*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F31*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F31*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F31*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F31*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F31*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F31*F40r*c0*f4*u0*u1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F30*F45r*c5*f4*u0
     & *u5*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F30*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F30*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F30*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F30*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q*PC_uh*i_*qq*ssm*svp*F30*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F45*F35r*c5*f3*u5**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F45*F34r*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F45*F33r*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F45*F32r*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F45*F31r*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F45*F30r*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F44*F35r*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F44*F34r*c4*f3*u4**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F44*F33r*c3*f3*u3
     & *u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F44*F32r*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F44*F31r*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F44*F30r*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F43*F35r*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F43*F34r*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F43*F33r*c3*f3*u3**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F43*F32r*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F43*F31r*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F43*F30r*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F42*F35r*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F42*F34r*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F42*F33r*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F42*F32r*c2*f3*u2**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F42*F31r*c1*f3*u1
     & *u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F42*F30r*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F41*F35r*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F41*F34r*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F41*F33r*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F41*F32r*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F41*F31r*c1*f3*u1**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F41*F30r*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F40*F35r*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F40*F34r*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F40*F33r*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F40*F32r*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F40*F31r*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F40*F30r*c0*f3*u0**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F35*F45r*c5*f4*
     & u5**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F35*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F35*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F35*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F35*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F35*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F34*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F34*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F34*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F34*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F34*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F34*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F33*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F33*F44r*c4*f4*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F33*F43r*c3*f4*
     & u3**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F33*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F33*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F33*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F32*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F32*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F32*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F32*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F32*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F32*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F31*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F31*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F31*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F31*F42r*c2*f4*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F31*F41r*c1*f4*
     & u1**2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F31*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F30*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F30*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F30*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F30*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F30*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC_uh*i_*qq*ssp*svm*F30*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F15i*c5*f1*u5**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F15i*c5*f1*u5**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F14i*c4*f1*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F13i*c3*f1*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F13i*c3*f1*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F45*F12i*c2*f1*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F11i*c1*f1*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F10i*c0*f1*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F45*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F15i*c5*f1*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F14i*c4*f1*u4**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F14i*c4*f1*u4**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F13i*c3*f1*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F12i*c2*f1*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F12i*c2*f1*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F11i*c1*f1*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F44*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F10i*c0*f1*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F44*F10i*c0*f1*u0*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F15i*c5*f1*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F15i*c5*f1*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F14i*c4*f1*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F13i*c3*f1*u3**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F13i*c3*f1*u3**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F12i*c2*f1*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F11i*c1*f1*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F10i*c0*f1*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F43*F10i*c0*f1*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F42*F15i*c5*f1*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F14i*c4*f1*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F13i*c3*f1*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F12i*c2*f1*u2**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F12i*c2*f1*u2**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F11i*c1*f1*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F10i*c0*f1*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F42*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F15i*c5*f1*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F15i*c5*f1*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F14i*c4*f1*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F41*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F13i*c3*f1*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F13i*c3*f1*u1*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F12i*c2*f1*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F12i*c2*f1*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F11i*c1*f1*u1**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F11i*c1*f1*u1**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F10i*c0*f1*u0*u1*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F41*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F15i*c5*f1*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F14i*c4*f1*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F13i*c3*f1*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F13i*c3*f1*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F40*F12i*c2*f1*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F11i*c1*f1*u0*u1*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F10i*c0*f1*u0**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F40*F10i*c0*f1*u0**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F45i*c5*f4*u5**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F45i*c5*f4*u5**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F41i*c1*f4*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F15*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F15*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F44i*c4*f4*u4**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F14*F40i*c0*f4*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F13*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F43i*c3*f4*u3**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F43i*c3*f4*u3**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F13*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F44i*c4*f4*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F12*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F42i*c2*f4*u2**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F12*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F43i*c3*f4*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F11*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F41i*c1*f4*u1**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F41i*c1*f4*u1**2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F11*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F43i*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F41i*c1*f4*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*svp*svm*F10*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F40i*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*q_uh*svp*svm*F10*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F35*F15i*c5*f1*u5**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F35*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F35*F13i*c3*f1*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F35*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F35*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F35*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F34*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F34*F14i*c4*f1*u4**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F34*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F34*F12i*c2*f1*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F34*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F34*F10i*c0*f1*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssm*svp*F33*F15i*c5*f1*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F33*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F33*F13i*c3*f1*u3**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F33*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F33*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F33*F10i*c0*f1*u0*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F32*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F32*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F32*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F32*F12i*c2*f1*u2**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F32*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F32*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F31*F15i*c5*f1*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F31*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F31*F13i*c3*f1*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssm*svp*F31*F12i*c2*f1*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F31*F11i*c1*f1*u1**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F31*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F30*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F30*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F30*F13i*c3*f1*u0*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F30*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F30*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F30*F10i*c0*f1*u0**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F15*F35i*c5*f3*u5**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F15*F34i*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F15*F33i*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F15*F32i*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F15*F31i*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F15*F30i*c0*f3*u0*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssm*svp*F14*F35i*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F14*F34i*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F14*F33i*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F14*F32i*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F14*F31i*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F14*F30i*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F13*F35i*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F13*F34i*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F13*F33i*c3*f3*u3**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F13*F32i*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F13*F31i*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F13*F30i*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F12*F35i*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F12*F34i*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F12*F33i*c3*f3*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssm*svp*F12*F32i*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F12*F31i*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F12*F30i*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F11*F35i*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F11*F34i*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F11*F33i*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F11*F32i*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F11*F31i*c1*f3*u1**2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F11*F30i*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F10*F35i*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F10*F34i*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F10*F33i*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F10*F32i*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F10*F31i*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*q_uh*ssm*svp*F10*F30i*c0*f3*u0**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssp*svm*F35*F15i*c5*f1*u5**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F35*F14i*c4*f1*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F35*F13i*c3*f1*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F35*F12i*c2*f1*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F35*F11i*c1*f1*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F35*F10i*c0*f1*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F34*F15i*c5*f1*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F34*F14i*c4*f1*u4**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F34*F13i*c3*f1*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F34*F12i*c2*f1*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F34*F11i*c1*f1*u1*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F34*F10i*c0*f1*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F33*F15i*c5*f1*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F33*F14i*c4*f1*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F33*F13i*c3*f1*u3**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssp*svm*F33*F12i*c2*f1*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F33*F11i*c1*f1*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F33*F10i*c0*f1*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F32*F15i*c5*f1*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F32*F14i*c4*f1*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F32*F13i*c3*f1*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F32*F12i*c2*f1*u2**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F32*F11i*c1*f1*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F32*F10i*c0*f1*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F31*F15i*c5*f1*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F31*F14i*c4*f1*u1*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F31*F13i*c3*f1*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F31*F12i*c2*f1*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F31*F11i*c1*f1*u1**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F31*F10i*c0*f1*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssp*svm*F30*F15i*c5*f1*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F30*F14i*c4*f1*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F30*F13i*c3*f1*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F30*F12i*c2*f1*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F30*F11i*c1*f1*u0*u1*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F30*F10i*c0*f1*u0**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F15*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F15*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F15*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F15*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F15*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F15*F30i*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F14*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F14*F34i*c4*f3*u4**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F14*F33i*c3*f3*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssp*svm*F14*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F14*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F14*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F13*F35i*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F13*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F13*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F13*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F13*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F13*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F12*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F12*F34i*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F12*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F12*F32i*c2*f3*u2**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F12*F31i*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F12*F30i*c0*f3*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*q_uh*ssp*svm*F11*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F11*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F11*F33i*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F11*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F11*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F11*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F10*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F10*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F10*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F10*F32i*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F10*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q*q_uh*ssp*svm*F10*F30i*c0*f3*u0**2*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F35*F25i*c5*f2*u5**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F35*F24i*c4*f2*u4*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F35*F23i*c3*f2*u3*u5*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q*q_uh*PC2*svp*svm*F35*F22i*c2*f2*u2*u5
     & *w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F35*F21i*c1*f2*u1*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F35*F20i*c0*f2*u0*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F34*F25i*c5*f2*u4*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F34*F24i*c4*f2*u4**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F34*F23i*c3*f2*u3*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F34*F22i*c2*f2*u2*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F34*F21i*c1*f2*u1*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F34*F20i*c0*f2*u0*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F33*F25i*c5*f2*u3*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F33*F24i*c4*f2*u3*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F33*F23i*c3*f2*u3**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F33*F22i*c2*f2*u2*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F33*F21i*c1*f2*u1*u3*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q*q_uh*PC2*svp*svm*F33*F20i*c0*f2*u0*u3
     & *w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F32*F25i*c5*f2*u2*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F32*F24i*c4*f2*u2*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F32*F23i*c3*f2*u2*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F32*F22i*c2*f2*u2**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F32*F21i*c1*f2*u1*u2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F32*F20i*c0*f2*u0*u2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F31*F25i*c5*f2*u1*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F31*F24i*c4*f2*u1*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F31*F23i*c3*f2*u1*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F31*F22i*c2*f2*u1*u2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F31*F21i*c1*f2*u1**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F31*F20i*c0*f2*u0*u1*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F30*F25i*c5*f2*u0*u5*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q*q_uh*PC2*svp*svm*F30*F24i*c4*f2*u0*u4
     & *w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F30*F23i*c3*f2*u0*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F30*F22i*c2*f2*u0*u2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F30*F21i*c1*f2*u0*u1*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F30*F20i*c0*f2*u0**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F25*F35i*c5*f3*u5**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F25*F34i*c4*f3*u4*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F25*F33i*c3*f3*u3*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F25*F32i*c2*f3*u2*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F25*F31i*c1*f3*u1*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F25*F30i*c0*f3*u0*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F24*F35i*c5*f3*u4*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F24*F34i*c4*f3*u4**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F24*F33i*c3*f3*u3*u4*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q*q_uh*PC2*svp*svm*F24*F32i*c2*f3*u2*u4
     & *w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F24*F31i*c1*f3*u1*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F24*F30i*c0*f3*u0*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F23*F35i*c5*f3*u3*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F23*F34i*c4*f3*u3*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F23*F33i*c3*f3*u3**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F23*F32i*c2*f3*u2*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F23*F31i*c1*f3*u1*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F23*F30i*c0*f3*u0*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F22*F35i*c5*f3*u2*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F22*F34i*c4*f3*u2*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F22*F33i*c3*f3*u2*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F22*F32i*c2*f3*u2**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F22*F31i*c1*f3*u1*u2*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q*q_uh*PC2*svp*svm*F22*F30i*c0*f3*u0*u2
     & *w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F21*F35i*c5*f3*u1*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F21*F34i*c4*f3*u1*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F21*F33i*c3*f3*u1*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F21*F32i*c2*f3*u1*u2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F21*F31i*c1*f3*u1**2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F21*F30i*c0*f3*u0*u1*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F20*F35i*c5*f3*u0*u5*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F20*F34i*c4*f3*u0*u4*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F20*F33i*c3*f3*u0*u3*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F20*F32i*c2*f3*u0*u2*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F20*F31i*c1*f3*u0*u1*w1*w2
     &  + 8.D0*PC_q*q_uh*PC2*svp*svm*F20*F30i*c0*f3*u0**2*w1*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F15r*c5*f1*u5**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F15r*c5*f1*u5**2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F45*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F14r*c4*f1*u4**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F14r*c4*f1*u4**2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F44*F10r*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F13r*c3*f1*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F13r*c3*f1*u3**2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F43*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F12r*c2*f1*u2**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F12r*c2*f1*u2**2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F42*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F11r*c1*f1*u1**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F11r*c1*f1*u1**2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F41*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F10r*c0*f1*u0**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F40*F10r*c0*f1*u0**2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F45r*c5*f4*u5**2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F15*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F45r*c5*f4*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F45r*c5*f4*u4*u5*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F14*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F44r*c4*f4*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F44r*c4*f4*u3*u4*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F43r*c3*f4*u3**2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F42r*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F13*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F43r*c3*f4*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F43r*c3*f4*u2*u3*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F12*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F42r*c2*f4*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F42r*c2*f4*u1*u2*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F11*F40r*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F41r*c1*f4*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F41r*c1*f4*u0*u1*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q*q_uh*i_*svp*svm*F10*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F35*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F35*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F35*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F35*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F35*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F35*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F34*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F34*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F34*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F34*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F34*F11r*c1*f1*u1*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssm*svp*F34*F10r*c0*f1*u0*u4*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F33*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F33*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F33*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F33*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F33*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F33*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F32*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F32*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F32*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F32*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F32*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F32*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F31*F15r*c5*f1*u1*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssm*svp*F31*F14r*c4*f1*u1*u4*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F31*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F31*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F31*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F31*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F30*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F30*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F30*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F30*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F30*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F30*F10r*c0*f1*u0**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F15*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F15*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F15*F33r*c3*f3*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssm*svp*F15*F32r*c2*f3*u2*u5*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F15*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F15*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F14*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F14*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F14*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F14*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F14*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F14*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F13*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F13*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F13*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F13*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F13*F31r*c1*f3*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssm*svp*F13*F30r*c0*f3*u0*u3*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F12*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F12*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F12*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F12*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F12*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F12*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F11*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F11*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F11*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F11*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F11*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F11*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F10*F35r*c5*f3*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssm*svp*F10*F34r*c4*f3*u0*u4*
     & w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F10*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F10*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F10*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q*q_uh*i_*ssm*svp*F10*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F35*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F35*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F35*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F35*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F35*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F35*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F34*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F34*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F34*F13r*c3*f1*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssp*svm*F34*F12r*c2*f1*u2*u4*
     & w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F34*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F34*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F33*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F33*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F33*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F33*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F33*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F33*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F32*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F32*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F32*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F32*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F32*F11r*c1*f1*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssp*svm*F32*F10r*c0*f1*u0*u2*
     & w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F31*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F31*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F31*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F31*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F31*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F31*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F30*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F30*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F30*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F30*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F30*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F30*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F15*F35r*c5*f3*u5**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssp*svm*F15*F34r*c4*f3*u4*u5*
     & w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F15*F33r*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F15*F32r*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F15*F31r*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F15*F30r*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F14*F35r*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F14*F34r*c4*f3*u4**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F14*F33r*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F14*F32r*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F14*F31r*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F14*F30r*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F13*F35r*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F13*F34r*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F13*F33r*c3*f3*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssp*svm*F13*F32r*c2*f3*u2*u3*
     & w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F13*F31r*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F13*F30r*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F12*F35r*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F12*F34r*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F12*F33r*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F12*F32r*c2*f3*u2**2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F12*F31r*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F12*F30r*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F11*F35r*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F11*F34r*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F11*F33r*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F11*F32r*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F11*F31r*c1*f3*u1**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*q_uh*i_*ssp*svm*F11*F30r*c0*f3*u0*u1*
     & w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F10*F35r*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F10*F34r*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F10*F33r*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F10*F32r*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F10*F31r*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q*q_uh*i_*ssp*svm*F10*F30r*c0*f3*u0**2*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F35*F25r*c5*f2*u5**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F35*F24r*c4*f2*u4*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F35*F23r*c3*f2*u3*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F35*F22r*c2*f2*u2*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F35*F21r*c1*f2*u1*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F35*F20r*c0*f2*u0*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F34*F25r*c5*f2*u4*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F34*F24r*c4*f2*
     & u4**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F34*F23r*c3*f2*u3*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F34*F22r*c2*f2*u2*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F34*F21r*c1*f2*u1*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F34*F20r*c0*f2*u0*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F33*F25r*c5*f2*u3*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F33*F24r*c4*f2*u3*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F33*F23r*c3*f2*u3**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F33*F22r*c2*f2*u2*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F33*F21r*c1*f2*u1*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F33*F20r*c0*f2*u0*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F32*F25r*c5*f2*u2*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F32*F24r*c4*f2*u2*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F32*F23r*c3*f2*u2*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F32*F22r*c2*f2*
     & u2**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F32*F21r*c1*f2*u1*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F32*F20r*c0*f2*u0*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F31*F25r*c5*f2*u1*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F31*F24r*c4*f2*u1*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F31*F23r*c3*f2*u1*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F31*F22r*c2*f2*u1*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F31*F21r*c1*f2*u1**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F31*F20r*c0*f2*u0*u1*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F30*F25r*c5*f2*u0*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F30*F24r*c4*f2*u0*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F30*F23r*c3*f2*u0*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F30*F22r*c2*f2*u0*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F30*F21r*c1*f2*u0*u1*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F30*F20r*c0*f2*
     & u0**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F25*F35r*c5*f3*u5**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F25*F34r*c4*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F25*F33r*c3*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F25*F32r*c2*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F25*F31r*c1*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F25*F30r*c0*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F24*F35r*c5*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F24*F34r*c4*f3*u4**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F24*F33r*c3*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F24*F32r*c2*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F24*F31r*c1*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F24*F30r*c0*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F23*F35r*c5*f3*u3*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F23*F34r*c4*f3*u3
     & *u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F23*F33r*c3*f3*u3**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F23*F32r*c2*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F23*F31r*c1*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F23*F30r*c0*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F22*F35r*c5*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F22*F34r*c4*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F22*F33r*c3*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F22*F32r*c2*f3*u2**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F22*F31r*c1*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F22*F30r*c0*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F21*F35r*c5*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F21*F34r*c4*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F21*F33r*c3*f3*u1*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F21*F32r*c2*f3*u1
     & *u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F21*F31r*c1*f3*u1**2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F21*F30r*c0*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F20*F35r*c5*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F20*F34r*c4*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F20*F33r*c3*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F20*F32r*c2*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F20*F31r*c1*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q*q_uh*i_*PC2*svp*svm*F20*F30r*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q*svm*dsvp*F15*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q*svm*dsvp*F15*F15r*c5*f1*u5**2*w1
     &  + 4.D0*PC_q*svm*dsvp*F15*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F15*F14r*c4*f1*u4*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F15*F13r*c3*f1*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svm*dsvp*F15*F13r*c3*f1*u3*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F15*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F15*F12r*c2*f1*u2*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F15*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F15*F11r*c1*f1*u1*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F15*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F15*F10r*c0*f1*u0*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F14*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F14*F15r*c5*f1*u4*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F14*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q*svm*dsvp*F14*F14r*c4*f1*u4**2*w1
     &  + 4.D0*PC_q*svm*dsvp*F14*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F14*F13r*c3*f1*u3*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F14*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F14*F12r*c2*f1*u2*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dsvp*F14*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F14*F11r*c1*f1*u1*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F14*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F14*F10r*c0*f1*u0*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F13*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F13*F15r*c5*f1*u3*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F13*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F13*F14r*c4*f1*u3*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F13*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q*svm*dsvp*F13*F13r*c3*f1*u3**2*w1
     &  + 4.D0*PC_q*svm*dsvp*F13*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q*svm*dsvp*F13*F12r*c2*f1*u2*u3*w1
     &  + 4.D0*PC_q*svm*dsvp*F13*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q*svm*dsvp*F13*F11r*c1*f1*u1*u3*w1
     &  + 4.D0*PC_q*svm*dsvp*F13*F10r*c0*f1*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svm*dsvp*F13*F10r*c0*f1*u0*u3*w1
     &  + 4.D0*PC_q*svm*dsvp*F12*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F12*F15r*c5*f1*u2*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F12*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F12*F14r*c4*f1*u2*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F12*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q*svm*dsvp*F12*F13r*c3*f1*u2*u3*w1
     &  + 4.D0*PC_q*svm*dsvp*F12*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q*svm*dsvp*F12*F12r*c2*f1*u2**2*w1
     &  + 4.D0*PC_q*svm*dsvp*F12*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q*svm*dsvp*F12*F11r*c1*f1*u1*u2*w1
     &  + 4.D0*PC_q*svm*dsvp*F12*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q*svm*dsvp*F12*F10r*c0*f1*u0*u2*w1
     &  + 4.D0*PC_q*svm*dsvp*F11*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F11*F15r*c5*f1*u1*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dsvp*F11*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F11*F14r*c4*f1*u1*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F11*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q*svm*dsvp*F11*F13r*c3*f1*u1*u3*w1
     &  + 4.D0*PC_q*svm*dsvp*F11*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q*svm*dsvp*F11*F12r*c2*f1*u1*u2*w1
     &  + 4.D0*PC_q*svm*dsvp*F11*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q*svm*dsvp*F11*F11r*c1*f1*u1**2*w1
     &  + 4.D0*PC_q*svm*dsvp*F11*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q*svm*dsvp*F11*F10r*c0*f1*u0*u1*w1
     &  + 4.D0*PC_q*svm*dsvp*F10*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q*svm*dsvp*F10*F15r*c5*f1*u0*u5*w1
     &  + 4.D0*PC_q*svm*dsvp*F10*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q*svm*dsvp*F10*F14r*c4*f1*u0*u4*w1
     &  + 4.D0*PC_q*svm*dsvp*F10*F13r*c3*f1*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svm*dsvp*F10*F13r*c3*f1*u0*u3*w1
     &  + 4.D0*PC_q*svm*dsvp*F10*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q*svm*dsvp*F10*F12r*c2*f1*u0*u2*w1
     &  + 4.D0*PC_q*svm*dsvp*F10*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q*svm*dsvp*F10*F11r*c1*f1*u0*u1*w1
     &  + 4.D0*PC_q*svm*dsvp*F10*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC_q*svm*dsvp*F10*F10r*c0*f1*u0**2*w1
     &  + 4.D0*PC_q*svm*dssp*F25*F15r*c5*f1*u5**2
     &  + 4.D0*PC_q*svm*dssp*F25*F14r*c4*f1*u4*u5
     &  + 4.D0*PC_q*svm*dssp*F25*F13r*c3*f1*u3*u5
     &  + 4.D0*PC_q*svm*dssp*F25*F12r*c2*f1*u2*u5
     &  + 4.D0*PC_q*svm*dssp*F25*F11r*c1*f1*u1*u5
     &  + 4.D0*PC_q*svm*dssp*F25*F10r*c0*f1*u0*u5
     &  + 4.D0*PC_q*svm*dssp*F24*F15r*c5*f1*u4*u5
     &  + 4.D0*PC_q*svm*dssp*F24*F14r*c4*f1*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dssp*F24*F13r*c3*f1*u3*u4
     &  + 4.D0*PC_q*svm*dssp*F24*F12r*c2*f1*u2*u4
     &  + 4.D0*PC_q*svm*dssp*F24*F11r*c1*f1*u1*u4
     &  + 4.D0*PC_q*svm*dssp*F24*F10r*c0*f1*u0*u4
     &  + 4.D0*PC_q*svm*dssp*F23*F15r*c5*f1*u3*u5
     &  + 4.D0*PC_q*svm*dssp*F23*F14r*c4*f1*u3*u4
     &  + 4.D0*PC_q*svm*dssp*F23*F13r*c3*f1*u3**2
     &  + 4.D0*PC_q*svm*dssp*F23*F12r*c2*f1*u2*u3
     &  + 4.D0*PC_q*svm*dssp*F23*F11r*c1*f1*u1*u3
     &  + 4.D0*PC_q*svm*dssp*F23*F10r*c0*f1*u0*u3
     &  + 4.D0*PC_q*svm*dssp*F22*F15r*c5*f1*u2*u5
     &  + 4.D0*PC_q*svm*dssp*F22*F14r*c4*f1*u2*u4
     &  + 4.D0*PC_q*svm*dssp*F22*F13r*c3*f1*u2*u3
     &  + 4.D0*PC_q*svm*dssp*F22*F12r*c2*f1*u2**2
     &  + 4.D0*PC_q*svm*dssp*F22*F11r*c1*f1*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dssp*F22*F10r*c0*f1*u0*u2
     &  + 4.D0*PC_q*svm*dssp*F21*F15r*c5*f1*u1*u5
     &  + 4.D0*PC_q*svm*dssp*F21*F14r*c4*f1*u1*u4
     &  + 4.D0*PC_q*svm*dssp*F21*F13r*c3*f1*u1*u3
     &  + 4.D0*PC_q*svm*dssp*F21*F12r*c2*f1*u1*u2
     &  + 4.D0*PC_q*svm*dssp*F21*F11r*c1*f1*u1**2
     &  + 4.D0*PC_q*svm*dssp*F21*F10r*c0*f1*u0*u1
     &  + 4.D0*PC_q*svm*dssp*F20*F15r*c5*f1*u0*u5
     &  + 4.D0*PC_q*svm*dssp*F20*F14r*c4*f1*u0*u4
     &  + 4.D0*PC_q*svm*dssp*F20*F13r*c3*f1*u0*u3
     &  + 4.D0*PC_q*svm*dssp*F20*F12r*c2*f1*u0*u2
     &  + 4.D0*PC_q*svm*dssp*F20*F11r*c1*f1*u0*u1
     &  + 4.D0*PC_q*svm*dssp*F20*F10r*c0*f1*u0**2
     &  + 4.D0*PC_q*svm*dssp*F15*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q*svm*dssp*F15*F24r*c4*f2*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dssp*F15*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q*svm*dssp*F15*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q*svm*dssp*F15*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q*svm*dssp*F15*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q*svm*dssp*F14*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q*svm*dssp*F14*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q*svm*dssp*F14*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q*svm*dssp*F14*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q*svm*dssp*F14*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q*svm*dssp*F14*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q*svm*dssp*F13*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q*svm*dssp*F13*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q*svm*dssp*F13*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q*svm*dssp*F13*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q*svm*dssp*F13*F21r*c1*f2*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dssp*F13*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q*svm*dssp*F12*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q*svm*dssp*F12*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q*svm*dssp*F12*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q*svm*dssp*F12*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q*svm*dssp*F12*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q*svm*dssp*F12*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q*svm*dssp*F11*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q*svm*dssp*F11*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q*svm*dssp*F11*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q*svm*dssp*F11*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q*svm*dssp*F11*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q*svm*dssp*F11*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q*svm*dssp*F10*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q*svm*dssp*F10*F24r*c4*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q*svm*dssp*F10*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q*svm*dssp*F10*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q*svm*dssp*F10*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q*svm*dssp*F10*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q*svp*dsvm*F15*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q*svp*dsvm*F15*F15r*c5*f1*u5**2*w1
     &  + 4.D0*PC_q*svp*dsvm*F15*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F15*F14r*c4*f1*u4*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F15*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F15*F13r*c3*f1*u3*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F15*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F15*F12r*c2*f1*u2*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F15*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F15*F11r*c1*f1*u1*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F15*F10r*c0*f1*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dsvm*F15*F10r*c0*f1*u0*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F14*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F14*F15r*c5*f1*u4*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F14*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q*svp*dsvm*F14*F14r*c4*f1*u4**2*w1
     &  + 4.D0*PC_q*svp*dsvm*F14*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F14*F13r*c3*f1*u3*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F14*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F14*F12r*c2*f1*u2*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F14*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F14*F11r*c1*f1*u1*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F14*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F14*F10r*c0*f1*u0*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F13*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F13*F15r*c5*f1*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svp*dsvm*F13*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F13*F14r*c4*f1*u3*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F13*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q*svp*dsvm*F13*F13r*c3*f1*u3**2*w1
     &  + 4.D0*PC_q*svp*dsvm*F13*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q*svp*dsvm*F13*F12r*c2*f1*u2*u3*w1
     &  + 4.D0*PC_q*svp*dsvm*F13*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q*svp*dsvm*F13*F11r*c1*f1*u1*u3*w1
     &  + 4.D0*PC_q*svp*dsvm*F13*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q*svp*dsvm*F13*F10r*c0*f1*u0*u3*w1
     &  + 4.D0*PC_q*svp*dsvm*F12*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F12*F15r*c5*f1*u2*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F12*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F12*F14r*c4*f1*u2*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F12*F13r*c3*f1*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dsvm*F12*F13r*c3*f1*u2*u3*w1
     &  + 4.D0*PC_q*svp*dsvm*F12*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q*svp*dsvm*F12*F12r*c2*f1*u2**2*w1
     &  + 4.D0*PC_q*svp*dsvm*F12*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q*svp*dsvm*F12*F11r*c1*f1*u1*u2*w1
     &  + 4.D0*PC_q*svp*dsvm*F12*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q*svp*dsvm*F12*F10r*c0*f1*u0*u2*w1
     &  + 4.D0*PC_q*svp*dsvm*F11*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F11*F15r*c5*f1*u1*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F11*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F11*F14r*c4*f1*u1*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F11*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q*svp*dsvm*F11*F13r*c3*f1*u1*u3*w1
     &  + 4.D0*PC_q*svp*dsvm*F11*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q*svp*dsvm*F11*F12r*c2*f1*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*svp*dsvm*F11*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q*svp*dsvm*F11*F11r*c1*f1*u1**2*w1
     &  + 4.D0*PC_q*svp*dsvm*F11*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q*svp*dsvm*F11*F10r*c0*f1*u0*u1*w1
     &  + 4.D0*PC_q*svp*dsvm*F10*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q*svp*dsvm*F10*F15r*c5*f1*u0*u5*w1
     &  + 4.D0*PC_q*svp*dsvm*F10*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q*svp*dsvm*F10*F14r*c4*f1*u0*u4*w1
     &  + 4.D0*PC_q*svp*dsvm*F10*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q*svp*dsvm*F10*F13r*c3*f1*u0*u3*w1
     &  + 4.D0*PC_q*svp*dsvm*F10*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q*svp*dsvm*F10*F12r*c2*f1*u0*u2*w1
     &  + 4.D0*PC_q*svp*dsvm*F10*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q*svp*dsvm*F10*F11r*c1*f1*u0*u1*w1
     &  + 4.D0*PC_q*svp*dsvm*F10*F10r*c0*f1*u0**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dsvm*F10*F10r*c0*f1*u0**2*w1
     &  - 4.D0*PC_q*svp*dssm*F25*F15r*c5*f1*u5**2
     &  - 4.D0*PC_q*svp*dssm*F25*F14r*c4*f1*u4*u5
     &  - 4.D0*PC_q*svp*dssm*F25*F13r*c3*f1*u3*u5
     &  - 4.D0*PC_q*svp*dssm*F25*F12r*c2*f1*u2*u5
     &  - 4.D0*PC_q*svp*dssm*F25*F11r*c1*f1*u1*u5
     &  - 4.D0*PC_q*svp*dssm*F25*F10r*c0*f1*u0*u5
     &  - 4.D0*PC_q*svp*dssm*F24*F15r*c5*f1*u4*u5
     &  - 4.D0*PC_q*svp*dssm*F24*F14r*c4*f1*u4**2
     &  - 4.D0*PC_q*svp*dssm*F24*F13r*c3*f1*u3*u4
     &  - 4.D0*PC_q*svp*dssm*F24*F12r*c2*f1*u2*u4
     &  - 4.D0*PC_q*svp*dssm*F24*F11r*c1*f1*u1*u4
     &  - 4.D0*PC_q*svp*dssm*F24*F10r*c0*f1*u0*u4
     &  - 4.D0*PC_q*svp*dssm*F23*F15r*c5*f1*u3*u5
     &  - 4.D0*PC_q*svp*dssm*F23*F14r*c4*f1*u3*u4
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dssm*F23*F13r*c3*f1*u3**2
     &  - 4.D0*PC_q*svp*dssm*F23*F12r*c2*f1*u2*u3
     &  - 4.D0*PC_q*svp*dssm*F23*F11r*c1*f1*u1*u3
     &  - 4.D0*PC_q*svp*dssm*F23*F10r*c0*f1*u0*u3
     &  - 4.D0*PC_q*svp*dssm*F22*F15r*c5*f1*u2*u5
     &  - 4.D0*PC_q*svp*dssm*F22*F14r*c4*f1*u2*u4
     &  - 4.D0*PC_q*svp*dssm*F22*F13r*c3*f1*u2*u3
     &  - 4.D0*PC_q*svp*dssm*F22*F12r*c2*f1*u2**2
     &  - 4.D0*PC_q*svp*dssm*F22*F11r*c1*f1*u1*u2
     &  - 4.D0*PC_q*svp*dssm*F22*F10r*c0*f1*u0*u2
     &  - 4.D0*PC_q*svp*dssm*F21*F15r*c5*f1*u1*u5
     &  - 4.D0*PC_q*svp*dssm*F21*F14r*c4*f1*u1*u4
     &  - 4.D0*PC_q*svp*dssm*F21*F13r*c3*f1*u1*u3
     &  - 4.D0*PC_q*svp*dssm*F21*F12r*c2*f1*u1*u2
     &  - 4.D0*PC_q*svp*dssm*F21*F11r*c1*f1*u1**2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dssm*F21*F10r*c0*f1*u0*u1
     &  - 4.D0*PC_q*svp*dssm*F20*F15r*c5*f1*u0*u5
     &  - 4.D0*PC_q*svp*dssm*F20*F14r*c4*f1*u0*u4
     &  - 4.D0*PC_q*svp*dssm*F20*F13r*c3*f1*u0*u3
     &  - 4.D0*PC_q*svp*dssm*F20*F12r*c2*f1*u0*u2
     &  - 4.D0*PC_q*svp*dssm*F20*F11r*c1*f1*u0*u1
     &  - 4.D0*PC_q*svp*dssm*F20*F10r*c0*f1*u0**2
     &  - 4.D0*PC_q*svp*dssm*F15*F25r*c5*f2*u5**2
     &  - 4.D0*PC_q*svp*dssm*F15*F24r*c4*f2*u4*u5
     &  - 4.D0*PC_q*svp*dssm*F15*F23r*c3*f2*u3*u5
     &  - 4.D0*PC_q*svp*dssm*F15*F22r*c2*f2*u2*u5
     &  - 4.D0*PC_q*svp*dssm*F15*F21r*c1*f2*u1*u5
     &  - 4.D0*PC_q*svp*dssm*F15*F20r*c0*f2*u0*u5
     &  - 4.D0*PC_q*svp*dssm*F14*F25r*c5*f2*u4*u5
     &  - 4.D0*PC_q*svp*dssm*F14*F24r*c4*f2*u4**2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dssm*F14*F23r*c3*f2*u3*u4
     &  - 4.D0*PC_q*svp*dssm*F14*F22r*c2*f2*u2*u4
     &  - 4.D0*PC_q*svp*dssm*F14*F21r*c1*f2*u1*u4
     &  - 4.D0*PC_q*svp*dssm*F14*F20r*c0*f2*u0*u4
     &  - 4.D0*PC_q*svp*dssm*F13*F25r*c5*f2*u3*u5
     &  - 4.D0*PC_q*svp*dssm*F13*F24r*c4*f2*u3*u4
     &  - 4.D0*PC_q*svp*dssm*F13*F23r*c3*f2*u3**2
     &  - 4.D0*PC_q*svp*dssm*F13*F22r*c2*f2*u2*u3
     &  - 4.D0*PC_q*svp*dssm*F13*F21r*c1*f2*u1*u3
     &  - 4.D0*PC_q*svp*dssm*F13*F20r*c0*f2*u0*u3
     &  - 4.D0*PC_q*svp*dssm*F12*F25r*c5*f2*u2*u5
     &  - 4.D0*PC_q*svp*dssm*F12*F24r*c4*f2*u2*u4
     &  - 4.D0*PC_q*svp*dssm*F12*F23r*c3*f2*u2*u3
     &  - 4.D0*PC_q*svp*dssm*F12*F22r*c2*f2*u2**2
     &  - 4.D0*PC_q*svp*dssm*F12*F21r*c1*f2*u1*u2
     &
      traza1 = traza1 - 4.D0*PC_q*svp*dssm*F12*F20r*c0*f2*u0*u2
     &  - 4.D0*PC_q*svp*dssm*F11*F25r*c5*f2*u1*u5
     &  - 4.D0*PC_q*svp*dssm*F11*F24r*c4*f2*u1*u4
     &  - 4.D0*PC_q*svp*dssm*F11*F23r*c3*f2*u1*u3
     &  - 4.D0*PC_q*svp*dssm*F11*F22r*c2*f2*u1*u2
     &  - 4.D0*PC_q*svp*dssm*F11*F21r*c1*f2*u1**2
     &  - 4.D0*PC_q*svp*dssm*F11*F20r*c0*f2*u0*u1
     &  - 4.D0*PC_q*svp*dssm*F10*F25r*c5*f2*u0*u5
     &  - 4.D0*PC_q*svp*dssm*F10*F24r*c4*f2*u0*u4
     &  - 4.D0*PC_q*svp*dssm*F10*F23r*c3*f2*u0*u3
     &  - 4.D0*PC_q*svp*dssm*F10*F22r*c2*f2*u0*u2
     &  - 4.D0*PC_q*svp*dssm*F10*F21r*c1*f2*u0*u1
     &  - 4.D0*PC_q*svp*dssm*F10*F20r*c0*f2*u0**2
     &  - 4.D0*PC_q*ssm*dsvp*F25*F15r*c5*f1*u5**2
     &  - 4.D0*PC_q*ssm*dsvp*F25*F14r*c4*f1*u4*u5
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*dsvp*F25*F13r*c3*f1*u3*u5
     &  - 4.D0*PC_q*ssm*dsvp*F25*F12r*c2*f1*u2*u5
     &  - 4.D0*PC_q*ssm*dsvp*F25*F11r*c1*f1*u1*u5
     &  - 4.D0*PC_q*ssm*dsvp*F25*F10r*c0*f1*u0*u5
     &  - 4.D0*PC_q*ssm*dsvp*F24*F15r*c5*f1*u4*u5
     &  - 4.D0*PC_q*ssm*dsvp*F24*F14r*c4*f1*u4**2
     &  - 4.D0*PC_q*ssm*dsvp*F24*F13r*c3*f1*u3*u4
     &  - 4.D0*PC_q*ssm*dsvp*F24*F12r*c2*f1*u2*u4
     &  - 4.D0*PC_q*ssm*dsvp*F24*F11r*c1*f1*u1*u4
     &  - 4.D0*PC_q*ssm*dsvp*F24*F10r*c0*f1*u0*u4
     &  - 4.D0*PC_q*ssm*dsvp*F23*F15r*c5*f1*u3*u5
     &  - 4.D0*PC_q*ssm*dsvp*F23*F14r*c4*f1*u3*u4
     &  - 4.D0*PC_q*ssm*dsvp*F23*F13r*c3*f1*u3**2
     &  - 4.D0*PC_q*ssm*dsvp*F23*F12r*c2*f1*u2*u3
     &  - 4.D0*PC_q*ssm*dsvp*F23*F11r*c1*f1*u1*u3
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*dsvp*F23*F10r*c0*f1*u0*u3
     &  - 4.D0*PC_q*ssm*dsvp*F22*F15r*c5*f1*u2*u5
     &  - 4.D0*PC_q*ssm*dsvp*F22*F14r*c4*f1*u2*u4
     &  - 4.D0*PC_q*ssm*dsvp*F22*F13r*c3*f1*u2*u3
     &  - 4.D0*PC_q*ssm*dsvp*F22*F12r*c2*f1*u2**2
     &  - 4.D0*PC_q*ssm*dsvp*F22*F11r*c1*f1*u1*u2
     &  - 4.D0*PC_q*ssm*dsvp*F22*F10r*c0*f1*u0*u2
     &  - 4.D0*PC_q*ssm*dsvp*F21*F15r*c5*f1*u1*u5
     &  - 4.D0*PC_q*ssm*dsvp*F21*F14r*c4*f1*u1*u4
     &  - 4.D0*PC_q*ssm*dsvp*F21*F13r*c3*f1*u1*u3
     &  - 4.D0*PC_q*ssm*dsvp*F21*F12r*c2*f1*u1*u2
     &  - 4.D0*PC_q*ssm*dsvp*F21*F11r*c1*f1*u1**2
     &  - 4.D0*PC_q*ssm*dsvp*F21*F10r*c0*f1*u0*u1
     &  - 4.D0*PC_q*ssm*dsvp*F20*F15r*c5*f1*u0*u5
     &  - 4.D0*PC_q*ssm*dsvp*F20*F14r*c4*f1*u0*u4
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*dsvp*F20*F13r*c3*f1*u0*u3
     &  - 4.D0*PC_q*ssm*dsvp*F20*F12r*c2*f1*u0*u2
     &  - 4.D0*PC_q*ssm*dsvp*F20*F11r*c1*f1*u0*u1
     &  - 4.D0*PC_q*ssm*dsvp*F20*F10r*c0*f1*u0**2
     &  - 4.D0*PC_q*ssm*dsvp*F15*F25r*c5*f2*u5**2
     &  - 4.D0*PC_q*ssm*dsvp*F15*F24r*c4*f2*u4*u5
     &  - 4.D0*PC_q*ssm*dsvp*F15*F23r*c3*f2*u3*u5
     &  - 4.D0*PC_q*ssm*dsvp*F15*F22r*c2*f2*u2*u5
     &  - 4.D0*PC_q*ssm*dsvp*F15*F21r*c1*f2*u1*u5
     &  - 4.D0*PC_q*ssm*dsvp*F15*F20r*c0*f2*u0*u5
     &  - 4.D0*PC_q*ssm*dsvp*F14*F25r*c5*f2*u4*u5
     &  - 4.D0*PC_q*ssm*dsvp*F14*F24r*c4*f2*u4**2
     &  - 4.D0*PC_q*ssm*dsvp*F14*F23r*c3*f2*u3*u4
     &  - 4.D0*PC_q*ssm*dsvp*F14*F22r*c2*f2*u2*u4
     &  - 4.D0*PC_q*ssm*dsvp*F14*F21r*c1*f2*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*dsvp*F14*F20r*c0*f2*u0*u4
     &  - 4.D0*PC_q*ssm*dsvp*F13*F25r*c5*f2*u3*u5
     &  - 4.D0*PC_q*ssm*dsvp*F13*F24r*c4*f2*u3*u4
     &  - 4.D0*PC_q*ssm*dsvp*F13*F23r*c3*f2*u3**2
     &  - 4.D0*PC_q*ssm*dsvp*F13*F22r*c2*f2*u2*u3
     &  - 4.D0*PC_q*ssm*dsvp*F13*F21r*c1*f2*u1*u3
     &  - 4.D0*PC_q*ssm*dsvp*F13*F20r*c0*f2*u0*u3
     &  - 4.D0*PC_q*ssm*dsvp*F12*F25r*c5*f2*u2*u5
     &  - 4.D0*PC_q*ssm*dsvp*F12*F24r*c4*f2*u2*u4
     &  - 4.D0*PC_q*ssm*dsvp*F12*F23r*c3*f2*u2*u3
     &  - 4.D0*PC_q*ssm*dsvp*F12*F22r*c2*f2*u2**2
     &  - 4.D0*PC_q*ssm*dsvp*F12*F21r*c1*f2*u1*u2
     &  - 4.D0*PC_q*ssm*dsvp*F12*F20r*c0*f2*u0*u2
     &  - 4.D0*PC_q*ssm*dsvp*F11*F25r*c5*f2*u1*u5
     &  - 4.D0*PC_q*ssm*dsvp*F11*F24r*c4*f2*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*ssm*dsvp*F11*F23r*c3*f2*u1*u3
     &  - 4.D0*PC_q*ssm*dsvp*F11*F22r*c2*f2*u1*u2
     &  - 4.D0*PC_q*ssm*dsvp*F11*F21r*c1*f2*u1**2
     &  - 4.D0*PC_q*ssm*dsvp*F11*F20r*c0*f2*u0*u1
     &  - 4.D0*PC_q*ssm*dsvp*F10*F25r*c5*f2*u0*u5
     &  - 4.D0*PC_q*ssm*dsvp*F10*F24r*c4*f2*u0*u4
     &  - 4.D0*PC_q*ssm*dsvp*F10*F23r*c3*f2*u0*u3
     &  - 4.D0*PC_q*ssm*dsvp*F10*F22r*c2*f2*u0*u2
     &  - 4.D0*PC_q*ssm*dsvp*F10*F21r*c1*f2*u0*u1
     &  - 4.D0*PC_q*ssm*dsvp*F10*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q*ssp*dsvm*F25*F15r*c5*f1*u5**2
     &  + 4.D0*PC_q*ssp*dsvm*F25*F14r*c4*f1*u4*u5
     &  + 4.D0*PC_q*ssp*dsvm*F25*F13r*c3*f1*u3*u5
     &  + 4.D0*PC_q*ssp*dsvm*F25*F12r*c2*f1*u2*u5
     &  + 4.D0*PC_q*ssp*dsvm*F25*F11r*c1*f1*u1*u5
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*dsvm*F25*F10r*c0*f1*u0*u5
     &  + 4.D0*PC_q*ssp*dsvm*F24*F15r*c5*f1*u4*u5
     &  + 4.D0*PC_q*ssp*dsvm*F24*F14r*c4*f1*u4**2
     &  + 4.D0*PC_q*ssp*dsvm*F24*F13r*c3*f1*u3*u4
     &  + 4.D0*PC_q*ssp*dsvm*F24*F12r*c2*f1*u2*u4
     &  + 4.D0*PC_q*ssp*dsvm*F24*F11r*c1*f1*u1*u4
     &  + 4.D0*PC_q*ssp*dsvm*F24*F10r*c0*f1*u0*u4
     &  + 4.D0*PC_q*ssp*dsvm*F23*F15r*c5*f1*u3*u5
     &  + 4.D0*PC_q*ssp*dsvm*F23*F14r*c4*f1*u3*u4
     &  + 4.D0*PC_q*ssp*dsvm*F23*F13r*c3*f1*u3**2
     &  + 4.D0*PC_q*ssp*dsvm*F23*F12r*c2*f1*u2*u3
     &  + 4.D0*PC_q*ssp*dsvm*F23*F11r*c1*f1*u1*u3
     &  + 4.D0*PC_q*ssp*dsvm*F23*F10r*c0*f1*u0*u3
     &  + 4.D0*PC_q*ssp*dsvm*F22*F15r*c5*f1*u2*u5
     &  + 4.D0*PC_q*ssp*dsvm*F22*F14r*c4*f1*u2*u4
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*dsvm*F22*F13r*c3*f1*u2*u3
     &  + 4.D0*PC_q*ssp*dsvm*F22*F12r*c2*f1*u2**2
     &  + 4.D0*PC_q*ssp*dsvm*F22*F11r*c1*f1*u1*u2
     &  + 4.D0*PC_q*ssp*dsvm*F22*F10r*c0*f1*u0*u2
     &  + 4.D0*PC_q*ssp*dsvm*F21*F15r*c5*f1*u1*u5
     &  + 4.D0*PC_q*ssp*dsvm*F21*F14r*c4*f1*u1*u4
     &  + 4.D0*PC_q*ssp*dsvm*F21*F13r*c3*f1*u1*u3
     &  + 4.D0*PC_q*ssp*dsvm*F21*F12r*c2*f1*u1*u2
     &  + 4.D0*PC_q*ssp*dsvm*F21*F11r*c1*f1*u1**2
     &  + 4.D0*PC_q*ssp*dsvm*F21*F10r*c0*f1*u0*u1
     &  + 4.D0*PC_q*ssp*dsvm*F20*F15r*c5*f1*u0*u5
     &  + 4.D0*PC_q*ssp*dsvm*F20*F14r*c4*f1*u0*u4
     &  + 4.D0*PC_q*ssp*dsvm*F20*F13r*c3*f1*u0*u3
     &  + 4.D0*PC_q*ssp*dsvm*F20*F12r*c2*f1*u0*u2
     &  + 4.D0*PC_q*ssp*dsvm*F20*F11r*c1*f1*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*dsvm*F20*F10r*c0*f1*u0**2
     &  + 4.D0*PC_q*ssp*dsvm*F15*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q*ssp*dsvm*F15*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q*ssp*dsvm*F15*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q*ssp*dsvm*F15*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q*ssp*dsvm*F15*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q*ssp*dsvm*F15*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q*ssp*dsvm*F14*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q*ssp*dsvm*F14*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q*ssp*dsvm*F14*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q*ssp*dsvm*F14*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q*ssp*dsvm*F14*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q*ssp*dsvm*F14*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q*ssp*dsvm*F13*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q*ssp*dsvm*F13*F24r*c4*f2*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*dsvm*F13*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q*ssp*dsvm*F13*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q*ssp*dsvm*F13*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q*ssp*dsvm*F13*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q*ssp*dsvm*F12*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q*ssp*dsvm*F12*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q*ssp*dsvm*F12*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q*ssp*dsvm*F12*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q*ssp*dsvm*F12*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q*ssp*dsvm*F12*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q*ssp*dsvm*F11*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q*ssp*dsvm*F11*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q*ssp*dsvm*F11*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q*ssp*dsvm*F11*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q*ssp*dsvm*F11*F21r*c1*f2*u1**2
     &
      traza1 = traza1 + 4.D0*PC_q*ssp*dsvm*F11*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q*ssp*dsvm*F10*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q*ssp*dsvm*F10*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q*ssp*dsvm*F10*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q*ssp*dsvm*F10*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q*ssp*dsvm*F10*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q*ssp*dsvm*F10*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q*qq*svm*dssp*F35*F15r*c5*f1*u5**2
     &  + 4.D0*PC_q*qq*svm*dssp*F35*F14r*c4*f1*u4*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F35*F13r*c3*f1*u3*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F35*F12r*c2*f1*u2*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F35*F11r*c1*f1*u1*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F35*F10r*c0*f1*u0*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F34*F15r*c5*f1*u4*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F34*F14r*c4*f1*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q*qq*svm*dssp*F34*F13r*c3*f1*u3*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F34*F12r*c2*f1*u2*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F34*F11r*c1*f1*u1*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F34*F10r*c0*f1*u0*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F33*F15r*c5*f1*u3*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F33*F14r*c4*f1*u3*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F33*F13r*c3*f1*u3**2
     &  + 4.D0*PC_q*qq*svm*dssp*F33*F12r*c2*f1*u2*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F33*F11r*c1*f1*u1*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F33*F10r*c0*f1*u0*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F32*F15r*c5*f1*u2*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F32*F14r*c4*f1*u2*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F32*F13r*c3*f1*u2*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F32*F12r*c2*f1*u2**2
     &  + 4.D0*PC_q*qq*svm*dssp*F32*F11r*c1*f1*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q*qq*svm*dssp*F32*F10r*c0*f1*u0*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F31*F15r*c5*f1*u1*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F31*F14r*c4*f1*u1*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F31*F13r*c3*f1*u1*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F31*F12r*c2*f1*u1*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F31*F11r*c1*f1*u1**2
     &  + 4.D0*PC_q*qq*svm*dssp*F31*F10r*c0*f1*u0*u1
     &  + 4.D0*PC_q*qq*svm*dssp*F30*F15r*c5*f1*u0*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F30*F14r*c4*f1*u0*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F30*F13r*c3*f1*u0*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F30*F12r*c2*f1*u0*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F30*F11r*c1*f1*u0*u1
     &  + 4.D0*PC_q*qq*svm*dssp*F30*F10r*c0*f1*u0**2
     &  + 4.D0*PC_q*qq*svm*dssp*F15*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q*qq*svm*dssp*F15*F34r*c4*f3*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q*qq*svm*dssp*F15*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F15*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F15*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F15*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F14*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F14*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q*qq*svm*dssp*F14*F33r*c3*f3*u3*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F14*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F14*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F14*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F13*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F13*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F13*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q*qq*svm*dssp*F13*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F13*F31r*c1*f3*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q*qq*svm*dssp*F13*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F12*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F12*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F12*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F12*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q*qq*svm*dssp*F12*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F12*F30r*c0*f3*u0*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F11*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F11*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q*qq*svm*dssp*F11*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F11*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F11*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q*qq*svm*dssp*F11*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q*qq*svm*dssp*F10*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q*qq*svm*dssp*F10*F34r*c4*f3*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q*qq*svm*dssp*F10*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q*qq*svm*dssp*F10*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q*qq*svm*dssp*F10*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q*qq*svm*dssp*F10*F30r*c0*f3*u0**2
     &  - 4.D0*PC_q*qq*svp*dssm*F35*F15r*c5*f1*u5**2
     &  - 4.D0*PC_q*qq*svp*dssm*F35*F14r*c4*f1*u4*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F35*F13r*c3*f1*u3*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F35*F12r*c2*f1*u2*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F35*F11r*c1*f1*u1*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F35*F10r*c0*f1*u0*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F34*F15r*c5*f1*u4*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F34*F14r*c4*f1*u4**2
     &  - 4.D0*PC_q*qq*svp*dssm*F34*F13r*c3*f1*u3*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F34*F12r*c2*f1*u2*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F34*F11r*c1*f1*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*qq*svp*dssm*F34*F10r*c0*f1*u0*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F33*F15r*c5*f1*u3*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F33*F14r*c4*f1*u3*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F33*F13r*c3*f1*u3**2
     &  - 4.D0*PC_q*qq*svp*dssm*F33*F12r*c2*f1*u2*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F33*F11r*c1*f1*u1*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F33*F10r*c0*f1*u0*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F32*F15r*c5*f1*u2*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F32*F14r*c4*f1*u2*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F32*F13r*c3*f1*u2*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F32*F12r*c2*f1*u2**2
     &  - 4.D0*PC_q*qq*svp*dssm*F32*F11r*c1*f1*u1*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F32*F10r*c0*f1*u0*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F31*F15r*c5*f1*u1*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F31*F14r*c4*f1*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*qq*svp*dssm*F31*F13r*c3*f1*u1*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F31*F12r*c2*f1*u1*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F31*F11r*c1*f1*u1**2
     &  - 4.D0*PC_q*qq*svp*dssm*F31*F10r*c0*f1*u0*u1
     &  - 4.D0*PC_q*qq*svp*dssm*F30*F15r*c5*f1*u0*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F30*F14r*c4*f1*u0*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F30*F13r*c3*f1*u0*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F30*F12r*c2*f1*u0*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F30*F11r*c1*f1*u0*u1
     &  - 4.D0*PC_q*qq*svp*dssm*F30*F10r*c0*f1*u0**2
     &  - 4.D0*PC_q*qq*svp*dssm*F15*F35r*c5*f3*u5**2
     &  - 4.D0*PC_q*qq*svp*dssm*F15*F34r*c4*f3*u4*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F15*F33r*c3*f3*u3*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F15*F32r*c2*f3*u2*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F15*F31r*c1*f3*u1*u5
     &
      traza1 = traza1 - 4.D0*PC_q*qq*svp*dssm*F15*F30r*c0*f3*u0*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F14*F35r*c5*f3*u4*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F14*F34r*c4*f3*u4**2
     &  - 4.D0*PC_q*qq*svp*dssm*F14*F33r*c3*f3*u3*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F14*F32r*c2*f3*u2*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F14*F31r*c1*f3*u1*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F14*F30r*c0*f3*u0*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F13*F35r*c5*f3*u3*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F13*F34r*c4*f3*u3*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F13*F33r*c3*f3*u3**2
     &  - 4.D0*PC_q*qq*svp*dssm*F13*F32r*c2*f3*u2*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F13*F31r*c1*f3*u1*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F13*F30r*c0*f3*u0*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F12*F35r*c5*f3*u2*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F12*F34r*c4*f3*u2*u4
     &
      traza1 = traza1 - 4.D0*PC_q*qq*svp*dssm*F12*F33r*c3*f3*u2*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F12*F32r*c2*f3*u2**2
     &  - 4.D0*PC_q*qq*svp*dssm*F12*F31r*c1*f3*u1*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F12*F30r*c0*f3*u0*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F11*F35r*c5*f3*u1*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F11*F34r*c4*f3*u1*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F11*F33r*c3*f3*u1*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F11*F32r*c2*f3*u1*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F11*F31r*c1*f3*u1**2
     &  - 4.D0*PC_q*qq*svp*dssm*F11*F30r*c0*f3*u0*u1
     &  - 4.D0*PC_q*qq*svp*dssm*F10*F35r*c5*f3*u0*u5
     &  - 4.D0*PC_q*qq*svp*dssm*F10*F34r*c4*f3*u0*u4
     &  - 4.D0*PC_q*qq*svp*dssm*F10*F33r*c3*f3*u0*u3
     &  - 4.D0*PC_q*qq*svp*dssm*F10*F32r*c2*f3*u0*u2
     &  - 4.D0*PC_q*qq*svp*dssm*F10*F31r*c1*f3*u0*u1
     &
      traza1 = traza1 - 4.D0*PC_q*qq*svp*dssm*F10*F30r*c0*f3*u0**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F35*F15r*c5*f1*u5**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F35*F14r*c4*f1*u4*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F35*F13r*c3*f1*u3*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F35*F12r*c2*f1*u2*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F35*F11r*c1*f1*u1*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F35*F10r*c0*f1*u0*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F34*F15r*c5*f1*u4*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F34*F14r*c4*f1*u4**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F34*F13r*c3*f1*u3*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F34*F12r*c2*f1*u2*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F34*F11r*c1*f1*u1*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F34*F10r*c0*f1*u0*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F33*F15r*c5*f1*u3*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F33*F14r*c4*f1*u3*u4
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*dsvp*F33*F13r*c3*f1*u3**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F33*F12r*c2*f1*u2*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F33*F11r*c1*f1*u1*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F33*F10r*c0*f1*u0*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F32*F15r*c5*f1*u2*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F32*F14r*c4*f1*u2*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F32*F13r*c3*f1*u2*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F32*F12r*c2*f1*u2**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F32*F11r*c1*f1*u1*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F32*F10r*c0*f1*u0*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F31*F15r*c5*f1*u1*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F31*F14r*c4*f1*u1*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F31*F13r*c3*f1*u1*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F31*F12r*c2*f1*u1*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F31*F11r*c1*f1*u1**2
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*dsvp*F31*F10r*c0*f1*u0*u1
     &  - 4.D0*PC_q*qq*ssm*dsvp*F30*F15r*c5*f1*u0*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F30*F14r*c4*f1*u0*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F30*F13r*c3*f1*u0*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F30*F12r*c2*f1*u0*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F30*F11r*c1*f1*u0*u1
     &  - 4.D0*PC_q*qq*ssm*dsvp*F30*F10r*c0*f1*u0**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F15*F35r*c5*f3*u5**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F15*F34r*c4*f3*u4*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F15*F33r*c3*f3*u3*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F15*F32r*c2*f3*u2*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F15*F31r*c1*f3*u1*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F15*F30r*c0*f3*u0*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F14*F35r*c5*f3*u4*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F14*F34r*c4*f3*u4**2
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*dsvp*F14*F33r*c3*f3*u3*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F14*F32r*c2*f3*u2*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F14*F31r*c1*f3*u1*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F14*F30r*c0*f3*u0*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F13*F35r*c5*f3*u3*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F13*F34r*c4*f3*u3*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F13*F33r*c3*f3*u3**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F13*F32r*c2*f3*u2*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F13*F31r*c1*f3*u1*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F13*F30r*c0*f3*u0*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F12*F35r*c5*f3*u2*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F12*F34r*c4*f3*u2*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F12*F33r*c3*f3*u2*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F12*F32r*c2*f3*u2**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F12*F31r*c1*f3*u1*u2
     &
      traza1 = traza1 - 4.D0*PC_q*qq*ssm*dsvp*F12*F30r*c0*f3*u0*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F11*F35r*c5*f3*u1*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F11*F34r*c4*f3*u1*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F11*F33r*c3*f3*u1*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F11*F32r*c2*f3*u1*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F11*F31r*c1*f3*u1**2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F11*F30r*c0*f3*u0*u1
     &  - 4.D0*PC_q*qq*ssm*dsvp*F10*F35r*c5*f3*u0*u5
     &  - 4.D0*PC_q*qq*ssm*dsvp*F10*F34r*c4*f3*u0*u4
     &  - 4.D0*PC_q*qq*ssm*dsvp*F10*F33r*c3*f3*u0*u3
     &  - 4.D0*PC_q*qq*ssm*dsvp*F10*F32r*c2*f3*u0*u2
     &  - 4.D0*PC_q*qq*ssm*dsvp*F10*F31r*c1*f3*u0*u1
     &  - 4.D0*PC_q*qq*ssm*dsvp*F10*F30r*c0*f3*u0**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F35*F15r*c5*f1*u5**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F35*F14r*c4*f1*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*dsvm*F35*F13r*c3*f1*u3*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F35*F12r*c2*f1*u2*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F35*F11r*c1*f1*u1*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F35*F10r*c0*f1*u0*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F34*F15r*c5*f1*u4*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F34*F14r*c4*f1*u4**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F34*F13r*c3*f1*u3*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F34*F12r*c2*f1*u2*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F34*F11r*c1*f1*u1*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F34*F10r*c0*f1*u0*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F33*F15r*c5*f1*u3*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F33*F14r*c4*f1*u3*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F33*F13r*c3*f1*u3**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F33*F12r*c2*f1*u2*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F33*F11r*c1*f1*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*dsvm*F33*F10r*c0*f1*u0*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F32*F15r*c5*f1*u2*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F32*F14r*c4*f1*u2*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F32*F13r*c3*f1*u2*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F32*F12r*c2*f1*u2**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F32*F11r*c1*f1*u1*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F32*F10r*c0*f1*u0*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F31*F15r*c5*f1*u1*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F31*F14r*c4*f1*u1*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F31*F13r*c3*f1*u1*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F31*F12r*c2*f1*u1*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F31*F11r*c1*f1*u1**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F31*F10r*c0*f1*u0*u1
     &  + 4.D0*PC_q*qq*ssp*dsvm*F30*F15r*c5*f1*u0*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F30*F14r*c4*f1*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*dsvm*F30*F13r*c3*f1*u0*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F30*F12r*c2*f1*u0*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F30*F11r*c1*f1*u0*u1
     &  + 4.D0*PC_q*qq*ssp*dsvm*F30*F10r*c0*f1*u0**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F15*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F15*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F15*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F15*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F15*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F15*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F14*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F14*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F14*F33r*c3*f3*u3*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F14*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F14*F31r*c1*f3*u1*u4
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*dsvm*F14*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F13*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F13*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F13*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F13*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F13*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F13*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F12*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F12*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F12*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F12*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F12*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F12*F30r*c0*f3*u0*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F11*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F11*F34r*c4*f3*u1*u4
     &
      traza1 = traza1 + 4.D0*PC_q*qq*ssp*dsvm*F11*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F11*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F11*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F11*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q*qq*ssp*dsvm*F10*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q*qq*ssp*dsvm*F10*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q*qq*ssp*dsvm*F10*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q*qq*ssp*dsvm*F10*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q*qq*ssp*dsvm*F10*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q*qq*ssp*dsvm*F10*F30r*c0*f3*u0**2
     &  - 4.D0*PC_q*PC2*svm*dsvp*F25*F25r*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F25*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F25*F24r*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F25*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F25*F23r*c3*f2*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svm*dsvp*F25*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F25*F22r*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F25*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F25*F21r*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F25*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F25*F20r*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F25*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F24*F25r*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F24*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F24*F24r*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F24*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F24*F23r*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F24*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F24*F22r*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F24*F22r*c2*f2*u2*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svm*dsvp*F24*F21r*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F24*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F24*F20r*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F24*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F23*F25r*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F23*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F23*F24r*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F23*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F23*F23r*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F23*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F23*F22r*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F23*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F23*F21r*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F23*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F23*F20r*c0*f2*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svm*dsvp*F23*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F22*F25r*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F22*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F22*F24r*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F22*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F22*F23r*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F22*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F22*F22r*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F22*F22r*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F22*F21r*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F22*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F22*F20r*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F22*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F21*F25r*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F21*F25r*c5*f2*u1*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svm*dsvp*F21*F24r*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F21*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F21*F23r*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F21*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F21*F22r*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F21*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F21*F21r*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F21*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F21*F20r*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F21*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F20*F25r*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F20*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F20*F24r*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F20*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F20*F23r*c3*f2*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svm*dsvp*F20*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F20*F22r*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F20*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F20*F21r*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F20*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*svm*dsvp*F20*F20r*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*PC2*svm*dsvp*F20*F20r*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F25*F25r*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F25*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F25*F24r*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F25*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F25*F23r*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F25*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F25*F22r*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F25*F22r*c2*f2*u2*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svp*dsvm*F25*F21r*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F25*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F25*F20r*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F25*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F24*F25r*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F24*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F24*F24r*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F24*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F24*F23r*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F24*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F24*F22r*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F24*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F24*F21r*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F24*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F24*F20r*c0*f2*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svp*dsvm*F24*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F23*F25r*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F23*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F23*F24r*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F23*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F23*F23r*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F23*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F23*F22r*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F23*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F23*F21r*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F23*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F23*F20r*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F23*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F22*F25r*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F22*F25r*c5*f2*u2*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svp*dsvm*F22*F24r*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F22*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F22*F23r*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F22*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F22*F22r*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F22*F22r*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F22*F21r*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F22*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F22*F20r*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F22*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F21*F25r*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F21*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F21*F24r*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F21*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F21*F23r*c3*f2*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*svp*dsvm*F21*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F21*F22r*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F21*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F21*F21r*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F21*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F21*F20r*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F21*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F20*F25r*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F20*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F20*F24r*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F20*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F20*F23r*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F20*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F20*F22r*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F20*F22r*c2*f2*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*svp*dsvm*F20*F21r*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F20*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*svp*dsvm*F20*F20r*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*PC2*svp*dsvm*F20*F20r*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F45r*c5*f4*u5**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F40r*c0*f4*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F45*F40r*c0*f4*u0*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F44*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F45r*c5*f4*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F45r*c5*f4*u3*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F43r*c3*f4*u3**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F42r*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F43*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F44r*c4*f4*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F44r*c4*f4*u2*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F42*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F43r*c3*f4*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F43r*c3*f4*u1*u3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F41*F40r*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F42r*c2*f4*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F42r*c2*f4*u0*u2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F40*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F25r*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F24r*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F23r*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F22r*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F21r*c1*f2*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F21r*c1*f2*u1*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F20r*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F35*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F25r*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F24r*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F23r*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F22r*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F21r*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F20r*c0*f2*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F34*F20r*c0*f2*u0*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F25r*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F24r*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F23r*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F22r*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F21r*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F20r*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F33*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F25r*c5*f2*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F25r*c5*f2*u2*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F24r*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F23r*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F22r*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F22r*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F21r*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F20r*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F32*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F25r*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F24r*c4*f2*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F24r*c4*f2*u1*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F23r*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F22r*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F21r*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F20r*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F31*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F25r*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F24r*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F23r*c3*f2*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F23r*c3*f2*u0*u3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F22r*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F21r*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F20r*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F30*F20r*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F32r*c2*f3*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F32r*c2*f3*u2*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F25*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F35r*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F31r*c1*f3*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F31r*c1*f3*u1*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F24*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F34r*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F30r*c0*f3*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F23*F30r*c0*f3*u0*u3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F33r*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F22*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F35r*c5*f3*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F35r*c5*f3*u1*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F32r*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F21*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F34r*c4*f3*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F34r*c4*f3*u0*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F31r*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q*PC2*qq*svm*dsvp*F20*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F45*F35r*c5*f3*u5**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F45*F34r*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F45*F33r*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F45*F32r*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F45*F31r*c1*f3*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*svm*dssp*F45*F30r*c0*f3*u0*u5*
     & w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F44*F35r*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F44*F34r*c4*f3*u4**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F44*F33r*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F44*F32r*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F44*F31r*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F44*F30r*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F43*F35r*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F43*F34r*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F43*F33r*c3*f3*u3**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F43*F32r*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F43*F31r*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F43*F30r*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F42*F35r*c5*f3*u2*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*svm*dssp*F42*F34r*c4*f3*u2*u4*
     & w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F42*F33r*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F42*F32r*c2*f3*u2**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F42*F31r*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F42*F30r*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F41*F35r*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F41*F34r*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F41*F33r*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F41*F32r*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F41*F31r*c1*f3*u1**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F41*F30r*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F40*F35r*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F40*F34r*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F40*F33r*c3*f3*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*svm*dssp*F40*F32r*c2*f3*u0*u2*
     & w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F40*F31r*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F40*F30r*c0*f3*u0**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F35*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F35*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F35*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F35*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F35*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F35*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F34*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F34*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F34*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F34*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F34*F41r*c1*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*svm*dssp*F34*F40r*c0*f4*u0*u4*
     & w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F33*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F33*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F33*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F33*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F33*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F33*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F32*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F32*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F32*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F32*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F32*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F32*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F31*F45r*c5*f4*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*svm*dssp*F31*F44r*c4*f4*u1*u4*
     & w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F31*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F31*F42r*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F31*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F31*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F30*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F30*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F30*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F30*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F30*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*svm*dssp*F30*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F45r*c5*f4*u5**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F44r*c4*f4*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F44r*c4*f4*u4*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F45*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F43r*c3*f4*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F43r*c3*f4*u3*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F44*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F43r*c3*f4*u3**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F42r*c2*f4*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F42r*c2*f4*u2*u3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F43*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F41r*c1*f4*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F41r*c1*f4*u1*u2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F42*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F40r*c0*f4*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F41*F40r*c0*f4*u0*u1*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F40*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F25r*c5*f2*u5**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F25r*c5*f2*u5**2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F24r*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F23r*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F22r*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F21r*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F20r*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F35*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F25r*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F24r*c4*f2*u4**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F24r*c4*f2*u4**2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F23r*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F22r*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F21r*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F20r*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F34*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F25r*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F24r*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F23r*c3*f2*u3**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F23r*c3*f2*u3**2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F22r*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F21r*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F20r*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F33*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F25r*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F24r*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F23r*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F22r*c2*f2*u2**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F22r*c2*f2*u2**2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F21r*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F20r*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F32*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F25r*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F24r*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F23r*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F22r*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F21r*c1*f2*u1**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F21r*c1*f2*u1**2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F20r*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F31*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F25r*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F24r*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F23r*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F22r*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F21r*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F20r*c0*f2*u0**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F30*F20r*c0*f2*u0**2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F25*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F35r*c5*f3*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F35r*c5*f3*u4*u5*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F24*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F34r*c4*f3*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F34r*c4*f3*u3*u4*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F23*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F33r*c3*f3*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F33r*c3*f3*u2*u3*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F22*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F32r*c2*f3*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F32r*c2*f3*u1*u2*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F21*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F31r*c1*f3*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F31r*c1*f3*u0*u1*
     & w1
     &  - 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q*PC2*qq*svp*dsvm*F20*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F45*F35r*c5*f3*u5**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F45*F34r*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F45*F33r*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F45*F32r*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F45*F31r*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F45*F30r*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F44*F35r*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F44*F34r*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F44*F33r*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F44*F32r*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F44*F31r*c1*f3*u1*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dssm*F44*F30r*c0*f3*u0*u4*
     & w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F43*F35r*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F43*F34r*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F43*F33r*c3*f3*u3**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F43*F32r*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F43*F31r*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F43*F30r*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F42*F35r*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F42*F34r*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F42*F33r*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F42*F32r*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F42*F31r*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F42*F30r*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F41*F35r*c5*f3*u1*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dssm*F41*F34r*c4*f3*u1*u4*
     & w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F41*F33r*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F41*F32r*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F41*F31r*c1*f3*u1**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F41*F30r*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F40*F35r*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F40*F34r*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F40*F33r*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F40*F32r*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F40*F31r*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F40*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F35*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F35*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F35*F43r*c3*f4*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dssm*F35*F42r*c2*f4*u2*u5*
     & w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F35*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F35*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F34*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F34*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F34*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F34*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F34*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F34*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F33*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F33*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F33*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F33*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F33*F41r*c1*f4*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dssm*F33*F40r*c0*f4*u0*u3*
     & w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F32*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F32*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F32*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F32*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F32*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F32*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F31*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F31*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F31*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F31*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F31*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F31*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F30*F45r*c5*f4*u0*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*svp*dssm*F30*F44r*c4*f4*u0*u4*
     & w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F30*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F30*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F30*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*svp*dssm*F30*F40r*c0*f4*u0**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F45*F35r*c5*f3*u5**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F45*F34r*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F45*F33r*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F45*F32r*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F45*F31r*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F45*F30r*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F44*F35r*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F44*F34r*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F44*F33r*c3*f3*u3*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*dsvp*F44*F32r*c2*f3*u2*u4*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F44*F31r*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F44*F30r*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F43*F35r*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F43*F34r*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F43*F33r*c3*f3*u3**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F43*F32r*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F43*F31r*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F43*F30r*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F42*F35r*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F42*F34r*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F42*F33r*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F42*F32r*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F42*F31r*c1*f3*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*dsvp*F42*F30r*c0*f3*u0*u2*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F41*F35r*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F41*F34r*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F41*F33r*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F41*F32r*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F41*F31r*c1*f3*u1**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F41*F30r*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F40*F35r*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F40*F34r*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F40*F33r*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F40*F32r*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F40*F31r*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F40*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F35*F45r*c5*f4*u5**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*dsvp*F35*F44r*c4*f4*u4*u5*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F35*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F35*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F35*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F35*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F34*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F34*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F34*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F34*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F34*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F34*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F33*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F33*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F33*F43r*c3*f4*u3**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*dsvp*F33*F42r*c2*f4*u2*u3*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F33*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F33*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F32*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F32*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F32*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F32*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F32*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F32*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F31*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F31*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F31*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F31*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F31*F41r*c1*f4*u1**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*PC2*qq*ssm*dsvp*F31*F40r*c0*f4*u0*u1*
     & w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F30*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F30*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F30*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F30*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F30*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q*PC2*qq*ssm*dsvp*F30*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F45*F35r*c5*f3*u5**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F45*F34r*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F45*F33r*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F45*F32r*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F45*F31r*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F45*F30r*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F44*F35r*c5*f3*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*dsvm*F44*F34r*c4*f3*u4**2*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F44*F33r*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F44*F32r*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F44*F31r*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F44*F30r*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F43*F35r*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F43*F34r*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F43*F33r*c3*f3*u3**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F43*F32r*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F43*F31r*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F43*F30r*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F42*F35r*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F42*F34r*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F42*F33r*c3*f3*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*dsvm*F42*F32r*c2*f3*u2**2*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F42*F31r*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F42*F30r*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F41*F35r*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F41*F34r*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F41*F33r*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F41*F32r*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F41*F31r*c1*f3*u1**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F41*F30r*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F40*F35r*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F40*F34r*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F40*F33r*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F40*F32r*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F40*F31r*c1*f3*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*dsvm*F40*F30r*c0*f3*u0**2*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F35*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F35*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F35*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F35*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F35*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F35*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F34*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F34*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F34*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F34*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F34*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F34*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F33*F45r*c5*f4*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*dsvm*F33*F44r*c4*f4*u3*u4*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F33*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F33*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F33*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F33*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F32*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F32*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F32*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F32*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F32*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F32*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F31*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F31*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F31*F43r*c3*f4*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*PC2*qq*ssp*dsvm*F31*F42r*c2*f4*u1*u2*
     & w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F31*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F31*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F30*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F30*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F30*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F30*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F30*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q*PC2*qq*ssp*dsvm*F30*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*i_*svm*dsvp*F15*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F15*F15i*c5*f1*u5**2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F15*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F15*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F15*F13i*c3*f1*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svm*dsvp*F15*F13i*c3*f1*u3*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F15*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F15*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F15*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F15*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F15*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F15*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F14*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F14*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F14*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F14*F14i*c4*f1*u4**2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F14*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F14*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F14*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F14*F12i*c2*f1*u2*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dsvp*F14*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F14*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F14*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F14*F10i*c0*f1*u0*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F13*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F13*F15i*c5*f1*u3*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F13*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F13*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F13*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F13*F13i*c3*f1*u3**2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F13*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F13*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F13*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F13*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F13*F10i*c0*f1*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svm*dsvp*F13*F10i*c0*f1*u0*u3*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F12*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F12*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F12*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F12*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F12*F13i*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F12*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F12*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F12*F12i*c2*f1*u2**2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F12*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F12*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F12*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F12*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F11*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F11*F15i*c5*f1*u1*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dsvp*F11*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F11*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F11*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F11*F13i*c3*f1*u1*u3*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F11*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F11*F12i*c2*f1*u1*u2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F11*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F11*F11i*c1*f1*u1**2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F11*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F11*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F10*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F10*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F10*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F10*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F10*F13i*c3*f1*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svm*dsvp*F10*F13i*c3*f1*u0*u3*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F10*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F10*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F10*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F10*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*PC_q*i_*svm*dsvp*F10*F10i*c0*f1*u0**2*w2
     &  - 4.D0*PC_q*i_*svm*dsvp*F10*F10i*c0*f1*u0**2*w1
     &  + 4.D0*PC_q*i_*svm*dssp*F25*F15i*c5*f1*u5**2
     &  + 4.D0*PC_q*i_*svm*dssp*F25*F14i*c4*f1*u4*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F25*F13i*c3*f1*u3*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F25*F12i*c2*f1*u2*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F25*F11i*c1*f1*u1*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F25*F10i*c0*f1*u0*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F24*F15i*c5*f1*u4*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F24*F14i*c4*f1*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dssp*F24*F13i*c3*f1*u3*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F24*F12i*c2*f1*u2*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F24*F11i*c1*f1*u1*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F24*F10i*c0*f1*u0*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F23*F15i*c5*f1*u3*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F23*F14i*c4*f1*u3*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F23*F13i*c3*f1*u3**2
     &  + 4.D0*PC_q*i_*svm*dssp*F23*F12i*c2*f1*u2*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F23*F11i*c1*f1*u1*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F23*F10i*c0*f1*u0*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F22*F15i*c5*f1*u2*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F22*F14i*c4*f1*u2*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F22*F13i*c3*f1*u2*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F22*F12i*c2*f1*u2**2
     &  + 4.D0*PC_q*i_*svm*dssp*F22*F11i*c1*f1*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dssp*F22*F10i*c0*f1*u0*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F21*F15i*c5*f1*u1*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F21*F14i*c4*f1*u1*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F21*F13i*c3*f1*u1*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F21*F12i*c2*f1*u1*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F21*F11i*c1*f1*u1**2
     &  + 4.D0*PC_q*i_*svm*dssp*F21*F10i*c0*f1*u0*u1
     &  + 4.D0*PC_q*i_*svm*dssp*F20*F15i*c5*f1*u0*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F20*F14i*c4*f1*u0*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F20*F13i*c3*f1*u0*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F20*F12i*c2*f1*u0*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F20*F11i*c1*f1*u0*u1
     &  + 4.D0*PC_q*i_*svm*dssp*F20*F10i*c0*f1*u0**2
     &  + 4.D0*PC_q*i_*svm*dssp*F15*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q*i_*svm*dssp*F15*F24i*c4*f2*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dssp*F15*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F15*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F15*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F15*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F14*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F14*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q*i_*svm*dssp*F14*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F14*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F14*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F14*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F13*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F13*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F13*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q*i_*svm*dssp*F13*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F13*F21i*c1*f2*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dssp*F13*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F12*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F12*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F12*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F12*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q*i_*svm*dssp*F12*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F12*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F11*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F11*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q*i_*svm*dssp*F11*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F11*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F11*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q*i_*svm*dssp*F11*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q*i_*svm*dssp*F10*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q*i_*svm*dssp*F10*F24i*c4*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svm*dssp*F10*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q*i_*svm*dssp*F10*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q*i_*svm*dssp*F10*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q*i_*svm*dssp*F10*F20i*c0*f2*u0**2
     &  + 4.D0*PC_q*i_*svp*dsvm*F15*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F15*F15i*c5*f1*u5**2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F15*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F15*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F15*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F15*F13i*c3*f1*u3*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F15*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F15*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F15*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F15*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F15*F10i*c0*f1*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dsvm*F15*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F14*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F14*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F14*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F14*F14i*c4*f1*u4**2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F14*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F14*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F14*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F14*F12i*c2*f1*u2*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F14*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F14*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F14*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F14*F10i*c0*f1*u0*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F13*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F13*F15i*c5*f1*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svp*dsvm*F13*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F13*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F13*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F13*F13i*c3*f1*u3**2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F13*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F13*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F13*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F13*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F13*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F13*F10i*c0*f1*u0*u3*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F12*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F12*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F12*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F12*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F12*F13i*c3*f1*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dsvm*F12*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F12*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F12*F12i*c2*f1*u2**2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F12*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F12*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F12*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F12*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F11*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F11*F15i*c5*f1*u1*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F11*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F11*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F11*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F11*F13i*c3*f1*u1*u3*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F11*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F11*F12i*c2*f1*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*svp*dsvm*F11*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F11*F11i*c1*f1*u1**2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F11*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F11*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F10*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F10*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F10*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F10*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F10*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F10*F13i*c3*f1*u0*u3*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F10*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F10*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F10*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q*i_*svp*dsvm*F10*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*PC_q*i_*svp*dsvm*F10*F10i*c0*f1*u0**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dsvm*F10*F10i*c0*f1*u0**2*w1
     &  - 4.D0*PC_q*i_*svp*dssm*F25*F15i*c5*f1*u5**2
     &  - 4.D0*PC_q*i_*svp*dssm*F25*F14i*c4*f1*u4*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F25*F13i*c3*f1*u3*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F25*F12i*c2*f1*u2*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F25*F11i*c1*f1*u1*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F25*F10i*c0*f1*u0*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F24*F15i*c5*f1*u4*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F24*F14i*c4*f1*u4**2
     &  - 4.D0*PC_q*i_*svp*dssm*F24*F13i*c3*f1*u3*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F24*F12i*c2*f1*u2*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F24*F11i*c1*f1*u1*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F24*F10i*c0*f1*u0*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F23*F15i*c5*f1*u3*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F23*F14i*c4*f1*u3*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dssm*F23*F13i*c3*f1*u3**2
     &  - 4.D0*PC_q*i_*svp*dssm*F23*F12i*c2*f1*u2*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F23*F11i*c1*f1*u1*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F23*F10i*c0*f1*u0*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F22*F15i*c5*f1*u2*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F22*F14i*c4*f1*u2*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F22*F13i*c3*f1*u2*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F22*F12i*c2*f1*u2**2
     &  - 4.D0*PC_q*i_*svp*dssm*F22*F11i*c1*f1*u1*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F22*F10i*c0*f1*u0*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F21*F15i*c5*f1*u1*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F21*F14i*c4*f1*u1*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F21*F13i*c3*f1*u1*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F21*F12i*c2*f1*u1*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F21*F11i*c1*f1*u1**2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dssm*F21*F10i*c0*f1*u0*u1
     &  - 4.D0*PC_q*i_*svp*dssm*F20*F15i*c5*f1*u0*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F20*F14i*c4*f1*u0*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F20*F13i*c3*f1*u0*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F20*F12i*c2*f1*u0*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F20*F11i*c1*f1*u0*u1
     &  - 4.D0*PC_q*i_*svp*dssm*F20*F10i*c0*f1*u0**2
     &  - 4.D0*PC_q*i_*svp*dssm*F15*F25i*c5*f2*u5**2
     &  - 4.D0*PC_q*i_*svp*dssm*F15*F24i*c4*f2*u4*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F15*F23i*c3*f2*u3*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F15*F22i*c2*f2*u2*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F15*F21i*c1*f2*u1*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F15*F20i*c0*f2*u0*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F14*F25i*c5*f2*u4*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F14*F24i*c4*f2*u4**2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dssm*F14*F23i*c3*f2*u3*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F14*F22i*c2*f2*u2*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F14*F21i*c1*f2*u1*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F14*F20i*c0*f2*u0*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F13*F25i*c5*f2*u3*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F13*F24i*c4*f2*u3*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F13*F23i*c3*f2*u3**2
     &  - 4.D0*PC_q*i_*svp*dssm*F13*F22i*c2*f2*u2*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F13*F21i*c1*f2*u1*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F13*F20i*c0*f2*u0*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F12*F25i*c5*f2*u2*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F12*F24i*c4*f2*u2*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F12*F23i*c3*f2*u2*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F12*F22i*c2*f2*u2**2
     &  - 4.D0*PC_q*i_*svp*dssm*F12*F21i*c1*f2*u1*u2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*svp*dssm*F12*F20i*c0*f2*u0*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F11*F25i*c5*f2*u1*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F11*F24i*c4*f2*u1*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F11*F23i*c3*f2*u1*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F11*F22i*c2*f2*u1*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F11*F21i*c1*f2*u1**2
     &  - 4.D0*PC_q*i_*svp*dssm*F11*F20i*c0*f2*u0*u1
     &  - 4.D0*PC_q*i_*svp*dssm*F10*F25i*c5*f2*u0*u5
     &  - 4.D0*PC_q*i_*svp*dssm*F10*F24i*c4*f2*u0*u4
     &  - 4.D0*PC_q*i_*svp*dssm*F10*F23i*c3*f2*u0*u3
     &  - 4.D0*PC_q*i_*svp*dssm*F10*F22i*c2*f2*u0*u2
     &  - 4.D0*PC_q*i_*svp*dssm*F10*F21i*c1*f2*u0*u1
     &  - 4.D0*PC_q*i_*svp*dssm*F10*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F25*F15i*c5*f1*u5**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F25*F14i*c4*f1*u4*u5
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*dsvp*F25*F13i*c3*f1*u3*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F25*F12i*c2*f1*u2*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F25*F11i*c1*f1*u1*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F25*F10i*c0*f1*u0*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F24*F15i*c5*f1*u4*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F24*F14i*c4*f1*u4**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F24*F13i*c3*f1*u3*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F24*F12i*c2*f1*u2*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F24*F11i*c1*f1*u1*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F24*F10i*c0*f1*u0*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F23*F15i*c5*f1*u3*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F23*F14i*c4*f1*u3*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F23*F13i*c3*f1*u3**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F23*F12i*c2*f1*u2*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F23*F11i*c1*f1*u1*u3
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*dsvp*F23*F10i*c0*f1*u0*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F22*F15i*c5*f1*u2*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F22*F14i*c4*f1*u2*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F22*F13i*c3*f1*u2*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F22*F12i*c2*f1*u2**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F22*F11i*c1*f1*u1*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F22*F10i*c0*f1*u0*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F21*F15i*c5*f1*u1*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F21*F14i*c4*f1*u1*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F21*F13i*c3*f1*u1*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F21*F12i*c2*f1*u1*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F21*F11i*c1*f1*u1**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F21*F10i*c0*f1*u0*u1
     &  - 4.D0*PC_q*i_*ssm*dsvp*F20*F15i*c5*f1*u0*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F20*F14i*c4*f1*u0*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*dsvp*F20*F13i*c3*f1*u0*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F20*F12i*c2*f1*u0*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F20*F11i*c1*f1*u0*u1
     &  - 4.D0*PC_q*i_*ssm*dsvp*F20*F10i*c0*f1*u0**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F15*F25i*c5*f2*u5**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F15*F24i*c4*f2*u4*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F15*F23i*c3*f2*u3*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F15*F22i*c2*f2*u2*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F15*F21i*c1*f2*u1*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F15*F20i*c0*f2*u0*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F14*F25i*c5*f2*u4*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F14*F24i*c4*f2*u4**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F14*F23i*c3*f2*u3*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F14*F22i*c2*f2*u2*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F14*F21i*c1*f2*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*dsvp*F14*F20i*c0*f2*u0*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F13*F25i*c5*f2*u3*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F13*F24i*c4*f2*u3*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F13*F23i*c3*f2*u3**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F13*F22i*c2*f2*u2*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F13*F21i*c1*f2*u1*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F13*F20i*c0*f2*u0*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F12*F25i*c5*f2*u2*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F12*F24i*c4*f2*u2*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F12*F23i*c3*f2*u2*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F12*F22i*c2*f2*u2**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F12*F21i*c1*f2*u1*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F12*F20i*c0*f2*u0*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F11*F25i*c5*f2*u1*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F11*F24i*c4*f2*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*ssm*dsvp*F11*F23i*c3*f2*u1*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F11*F22i*c2*f2*u1*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F11*F21i*c1*f2*u1**2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F11*F20i*c0*f2*u0*u1
     &  - 4.D0*PC_q*i_*ssm*dsvp*F10*F25i*c5*f2*u0*u5
     &  - 4.D0*PC_q*i_*ssm*dsvp*F10*F24i*c4*f2*u0*u4
     &  - 4.D0*PC_q*i_*ssm*dsvp*F10*F23i*c3*f2*u0*u3
     &  - 4.D0*PC_q*i_*ssm*dsvp*F10*F22i*c2*f2*u0*u2
     &  - 4.D0*PC_q*i_*ssm*dsvp*F10*F21i*c1*f2*u0*u1
     &  - 4.D0*PC_q*i_*ssm*dsvp*F10*F20i*c0*f2*u0**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F25*F15i*c5*f1*u5**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F25*F14i*c4*f1*u4*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F25*F13i*c3*f1*u3*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F25*F12i*c2*f1*u2*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F25*F11i*c1*f1*u1*u5
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*dsvm*F25*F10i*c0*f1*u0*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F24*F15i*c5*f1*u4*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F24*F14i*c4*f1*u4**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F24*F13i*c3*f1*u3*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F24*F12i*c2*f1*u2*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F24*F11i*c1*f1*u1*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F24*F10i*c0*f1*u0*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F23*F15i*c5*f1*u3*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F23*F14i*c4*f1*u3*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F23*F13i*c3*f1*u3**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F23*F12i*c2*f1*u2*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F23*F11i*c1*f1*u1*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F23*F10i*c0*f1*u0*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F22*F15i*c5*f1*u2*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F22*F14i*c4*f1*u2*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*dsvm*F22*F13i*c3*f1*u2*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F22*F12i*c2*f1*u2**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F22*F11i*c1*f1*u1*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F22*F10i*c0*f1*u0*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F21*F15i*c5*f1*u1*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F21*F14i*c4*f1*u1*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F21*F13i*c3*f1*u1*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F21*F12i*c2*f1*u1*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F21*F11i*c1*f1*u1**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F21*F10i*c0*f1*u0*u1
     &  + 4.D0*PC_q*i_*ssp*dsvm*F20*F15i*c5*f1*u0*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F20*F14i*c4*f1*u0*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F20*F13i*c3*f1*u0*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F20*F12i*c2*f1*u0*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F20*F11i*c1*f1*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*dsvm*F20*F10i*c0*f1*u0**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F15*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F15*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F15*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F15*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F15*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F15*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F14*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F14*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F14*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F14*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F14*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F14*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F13*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F13*F24i*c4*f2*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*dsvm*F13*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F13*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F13*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F13*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F12*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F12*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F12*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F12*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F12*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F12*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F11*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F11*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F11*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F11*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F11*F21i*c1*f2*u1**2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*ssp*dsvm*F11*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q*i_*ssp*dsvm*F10*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q*i_*ssp*dsvm*F10*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q*i_*ssp*dsvm*F10*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q*i_*ssp*dsvm*F10*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q*i_*ssp*dsvm*F10*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q*i_*ssp*dsvm*F10*F20i*c0*f2*u0**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F35*F15i*c5*f1*u5**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F35*F14i*c4*f1*u4*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F35*F13i*c3*f1*u3*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F35*F12i*c2*f1*u2*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F35*F11i*c1*f1*u1*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F35*F10i*c0*f1*u0*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F34*F15i*c5*f1*u4*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F34*F14i*c4*f1*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*svm*dssp*F34*F13i*c3*f1*u3*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F34*F12i*c2*f1*u2*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F34*F11i*c1*f1*u1*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F34*F10i*c0*f1*u0*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F33*F15i*c5*f1*u3*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F33*F14i*c4*f1*u3*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F33*F13i*c3*f1*u3**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F33*F12i*c2*f1*u2*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F33*F11i*c1*f1*u1*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F33*F10i*c0*f1*u0*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F32*F15i*c5*f1*u2*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F32*F14i*c4*f1*u2*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F32*F13i*c3*f1*u2*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F32*F12i*c2*f1*u2**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F32*F11i*c1*f1*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*svm*dssp*F32*F10i*c0*f1*u0*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F31*F15i*c5*f1*u1*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F31*F14i*c4*f1*u1*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F31*F13i*c3*f1*u1*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F31*F12i*c2*f1*u1*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F31*F11i*c1*f1*u1**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F31*F10i*c0*f1*u0*u1
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F30*F15i*c5*f1*u0*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F30*F14i*c4*f1*u0*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F30*F13i*c3*f1*u0*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F30*F12i*c2*f1*u0*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F30*F11i*c1*f1*u0*u1
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F30*F10i*c0*f1*u0**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F15*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F15*F34i*c4*f3*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*svm*dssp*F15*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F15*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F15*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F15*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F14*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F14*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F14*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F14*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F14*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F14*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F13*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F13*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F13*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F13*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F13*F31i*c1*f3*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*svm*dssp*F13*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F12*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F12*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F12*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F12*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F12*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F12*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F11*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F11*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F11*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F11*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F11*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F11*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F10*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F10*F34i*c4*f3*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*svm*dssp*F10*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F10*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F10*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q*i_*qq*svm*dssp*F10*F30i*c0*f3*u0**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F35*F15i*c5*f1*u5**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F35*F14i*c4*f1*u4*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F35*F13i*c3*f1*u3*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F35*F12i*c2*f1*u2*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F35*F11i*c1*f1*u1*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F35*F10i*c0*f1*u0*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F34*F15i*c5*f1*u4*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F34*F14i*c4*f1*u4**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F34*F13i*c3*f1*u3*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F34*F12i*c2*f1*u2*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F34*F11i*c1*f1*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*svp*dssm*F34*F10i*c0*f1*u0*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F33*F15i*c5*f1*u3*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F33*F14i*c4*f1*u3*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F33*F13i*c3*f1*u3**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F33*F12i*c2*f1*u2*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F33*F11i*c1*f1*u1*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F33*F10i*c0*f1*u0*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F32*F15i*c5*f1*u2*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F32*F14i*c4*f1*u2*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F32*F13i*c3*f1*u2*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F32*F12i*c2*f1*u2**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F32*F11i*c1*f1*u1*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F32*F10i*c0*f1*u0*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F31*F15i*c5*f1*u1*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F31*F14i*c4*f1*u1*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*svp*dssm*F31*F13i*c3*f1*u1*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F31*F12i*c2*f1*u1*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F31*F11i*c1*f1*u1**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F31*F10i*c0*f1*u0*u1
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F30*F15i*c5*f1*u0*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F30*F14i*c4*f1*u0*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F30*F13i*c3*f1*u0*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F30*F12i*c2*f1*u0*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F30*F11i*c1*f1*u0*u1
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F30*F10i*c0*f1*u0**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F15*F35i*c5*f3*u5**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F15*F34i*c4*f3*u4*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F15*F33i*c3*f3*u3*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F15*F32i*c2*f3*u2*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F15*F31i*c1*f3*u1*u5
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*svp*dssm*F15*F30i*c0*f3*u0*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F14*F35i*c5*f3*u4*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F14*F34i*c4*f3*u4**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F14*F33i*c3*f3*u3*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F14*F32i*c2*f3*u2*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F14*F31i*c1*f3*u1*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F14*F30i*c0*f3*u0*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F13*F35i*c5*f3*u3*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F13*F34i*c4*f3*u3*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F13*F33i*c3*f3*u3**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F13*F32i*c2*f3*u2*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F13*F31i*c1*f3*u1*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F13*F30i*c0*f3*u0*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F12*F35i*c5*f3*u2*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F12*F34i*c4*f3*u2*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*svp*dssm*F12*F33i*c3*f3*u2*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F12*F32i*c2*f3*u2**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F12*F31i*c1*f3*u1*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F12*F30i*c0*f3*u0*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F11*F35i*c5*f3*u1*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F11*F34i*c4*f3*u1*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F11*F33i*c3*f3*u1*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F11*F32i*c2*f3*u1*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F11*F31i*c1*f3*u1**2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F11*F30i*c0*f3*u0*u1
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F10*F35i*c5*f3*u0*u5
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F10*F34i*c4*f3*u0*u4
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F10*F33i*c3*f3*u0*u3
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F10*F32i*c2*f3*u0*u2
     &  - 4.D0*PC_q*i_*qq*svp*dssm*F10*F31i*c1*f3*u0*u1
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*svp*dssm*F10*F30i*c0*f3*u0**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F35*F15i*c5*f1*u5**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F35*F14i*c4*f1*u4*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F35*F13i*c3*f1*u3*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F35*F12i*c2*f1*u2*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F35*F11i*c1*f1*u1*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F35*F10i*c0*f1*u0*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F34*F15i*c5*f1*u4*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F34*F14i*c4*f1*u4**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F34*F13i*c3*f1*u3*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F34*F12i*c2*f1*u2*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F34*F11i*c1*f1*u1*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F34*F10i*c0*f1*u0*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F33*F15i*c5*f1*u3*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F33*F14i*c4*f1*u3*u4
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*dsvp*F33*F13i*c3*f1*u3**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F33*F12i*c2*f1*u2*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F33*F11i*c1*f1*u1*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F33*F10i*c0*f1*u0*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F32*F15i*c5*f1*u2*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F32*F14i*c4*f1*u2*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F32*F13i*c3*f1*u2*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F32*F12i*c2*f1*u2**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F32*F11i*c1*f1*u1*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F32*F10i*c0*f1*u0*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F31*F15i*c5*f1*u1*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F31*F14i*c4*f1*u1*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F31*F13i*c3*f1*u1*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F31*F12i*c2*f1*u1*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F31*F11i*c1*f1*u1**2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*dsvp*F31*F10i*c0*f1*u0*u1
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F30*F15i*c5*f1*u0*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F30*F14i*c4*f1*u0*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F30*F13i*c3*f1*u0*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F30*F12i*c2*f1*u0*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F30*F11i*c1*f1*u0*u1
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F30*F10i*c0*f1*u0**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F15*F35i*c5*f3*u5**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F15*F34i*c4*f3*u4*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F15*F33i*c3*f3*u3*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F15*F32i*c2*f3*u2*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F15*F31i*c1*f3*u1*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F15*F30i*c0*f3*u0*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F14*F35i*c5*f3*u4*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F14*F34i*c4*f3*u4**2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*dsvp*F14*F33i*c3*f3*u3*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F14*F32i*c2*f3*u2*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F14*F31i*c1*f3*u1*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F14*F30i*c0*f3*u0*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F13*F35i*c5*f3*u3*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F13*F34i*c4*f3*u3*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F13*F33i*c3*f3*u3**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F13*F32i*c2*f3*u2*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F13*F31i*c1*f3*u1*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F13*F30i*c0*f3*u0*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F12*F35i*c5*f3*u2*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F12*F34i*c4*f3*u2*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F12*F33i*c3*f3*u2*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F12*F32i*c2*f3*u2**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F12*F31i*c1*f3*u1*u2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*qq*ssm*dsvp*F12*F30i*c0*f3*u0*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F11*F35i*c5*f3*u1*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F11*F34i*c4*f3*u1*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F11*F33i*c3*f3*u1*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F11*F32i*c2*f3*u1*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F11*F31i*c1*f3*u1**2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F11*F30i*c0*f3*u0*u1
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F10*F35i*c5*f3*u0*u5
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F10*F34i*c4*f3*u0*u4
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F10*F33i*c3*f3*u0*u3
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F10*F32i*c2*f3*u0*u2
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F10*F31i*c1*f3*u0*u1
     &  - 4.D0*PC_q*i_*qq*ssm*dsvp*F10*F30i*c0*f3*u0**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F35*F15i*c5*f1*u5**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F35*F14i*c4*f1*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*dsvm*F35*F13i*c3*f1*u3*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F35*F12i*c2*f1*u2*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F35*F11i*c1*f1*u1*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F35*F10i*c0*f1*u0*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F34*F15i*c5*f1*u4*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F34*F14i*c4*f1*u4**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F34*F13i*c3*f1*u3*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F34*F12i*c2*f1*u2*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F34*F11i*c1*f1*u1*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F34*F10i*c0*f1*u0*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F33*F15i*c5*f1*u3*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F33*F14i*c4*f1*u3*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F33*F13i*c3*f1*u3**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F33*F12i*c2*f1*u2*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F33*F11i*c1*f1*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*dsvm*F33*F10i*c0*f1*u0*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F32*F15i*c5*f1*u2*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F32*F14i*c4*f1*u2*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F32*F13i*c3*f1*u2*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F32*F12i*c2*f1*u2**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F32*F11i*c1*f1*u1*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F32*F10i*c0*f1*u0*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F31*F15i*c5*f1*u1*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F31*F14i*c4*f1*u1*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F31*F13i*c3*f1*u1*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F31*F12i*c2*f1*u1*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F31*F11i*c1*f1*u1**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F31*F10i*c0*f1*u0*u1
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F30*F15i*c5*f1*u0*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F30*F14i*c4*f1*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*dsvm*F30*F13i*c3*f1*u0*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F30*F12i*c2*f1*u0*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F30*F11i*c1*f1*u0*u1
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F30*F10i*c0*f1*u0**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F15*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F15*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F15*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F15*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F15*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F15*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F14*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F14*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F14*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F14*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F14*F31i*c1*f3*u1*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*dsvm*F14*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F13*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F13*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F13*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F13*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F13*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F13*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F12*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F12*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F12*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F12*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F12*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F12*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F11*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F11*F34i*c4*f3*u1*u4
     &
      traza1 = traza1 + 4.D0*PC_q*i_*qq*ssp*dsvm*F11*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F11*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F11*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F11*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F10*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F10*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F10*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F10*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F10*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q*i_*qq*ssp*dsvm*F10*F30i*c0*f3*u0**2
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F25i*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F25i*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F24i*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F23i*c3*f2*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F23i*c3*f2*u3*u5*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F25*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F25i*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F24i*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F24i*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F23i*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F22i*c2*f2*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F22i*c2*f2*u2*u4*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F21i*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F24*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F24i*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F23i*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F23i*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F22i*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F21i*c1*f2*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F21i*c1*f2*u1*u3*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F23*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F23i*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F23i*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F22i*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F22i*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F21i*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F21i*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F20i*c0*f2*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svm*dsvp*F22*F20i*c0*f2*u0*u2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F24i*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F22i*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F21i*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F21i*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F20i*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F21*F20i*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F25i*c5*f2*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F25i*c5*f2*u0*u5*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F21i*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F20i*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*svm*dsvp*F20*F20i*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F25i*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F25i*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F24i*c4*f2*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F24i*c4*f2*u4*u5*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F23i*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F23i*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F25*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F25i*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F24i*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F24i*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F23i*c3*f2*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F23i*c3*f2*u3*u4*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F22i*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F22i*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F21i*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F24*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F24i*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F23i*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F23i*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F22i*c2*f2*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F22i*c2*f2*u2*u3*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F21i*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F21i*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F23*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F23i*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F23i*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F22i*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F22i*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F21i*c1*f2*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F21i*c1*f2*u1*u2*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F20i*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F22*F20i*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F24i*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F22i*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F21i*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F21i*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F20i*c0*f2*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*svp*dsvm*F21*F20i*c0*f2*u0*u1*
     & w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F25i*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F25i*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F21i*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F20i*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*svp*dsvm*F20*F20i*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F45i*c5*f4*u5**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F45i*c5*f4*
     & u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F45*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F44i*c4*f4*u4**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F44i*c4*f4*
     & u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F44*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F43i*c3*f4*u3**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F43i*c3*f4*
     & u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F43*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F42i*c2*f4*u2**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F42i*c2*f4*
     & u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F42*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F41i*c1*f4*u1**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F41i*c1*f4*
     & u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F41*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F40i*c0*f4*u0**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F40*F40i*c0*f4*
     & u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F25i*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F25i*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F24i*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F23i*c3*f2*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F23i*c3*f2*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F35*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F25i*c5*f2*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F25i*c5*f2*u4*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F24i*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F24i*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F23i*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F22i*c2*f2*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F22i*c2*f2*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F21i*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F34*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F24i*c4*f2*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F24i*c4*f2*u3*
     & u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F23i*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F23i*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F22i*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F21i*c1*f2*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F21i*c1*f2*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F33*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F23i*c3*f2*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F23i*c3*f2*u2*
     & u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F22i*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F22i*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F21i*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F21i*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F20i*c0*f2*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F32*F20i*c0*f2*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F24i*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F22i*c2*f2*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F22i*c2*f2*u1*
     & u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F21i*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F21i*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F20i*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F31*F20i*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F25i*c5*f2*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F25i*c5*f2*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F21i*c1*f2*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F21i*c1*f2*u0*
     & u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F20i*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F30*F20i*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F30i*c0*f3*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F25*F30i*c0*f3*u0*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F34i*c4*f3*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F33i*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F24*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F35i*c5*f3*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F35i*c5*f3*u3*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F23*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F34i*c4*f3*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F34i*c4*f3*u2*
     & u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F32i*c2*f3*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F31i*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F22*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F33i*c3*f3*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F33i*c3*f3*u1*
     & u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F21*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F32i*c2*f3*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F32i*c2*f3*u0*
     & u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F30i*c0*f3*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svm*dsvp*F20*F30i*c0*f3*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F45*F35i*c5*f3*u5**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F45*F34i*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F45*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F45*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F45*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F45*F30i*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F44*F35i*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F44*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F44*F33i*c3*f3*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F44*F32i*c2*f3*u2*
     & u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F44*F31i*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F44*F30i*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F43*F35i*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F43*F34i*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F43*F33i*c3*f3*u3**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F43*F32i*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F43*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F43*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F42*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F42*F34i*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F42*F33i*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F42*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F42*F31i*c1*f3*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F42*F30i*c0*f3*u0*
     & u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F41*F35i*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F41*F34i*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F41*F33i*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F41*F32i*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F41*F31i*c1*f3*u1**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F41*F30i*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F40*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F40*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F40*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F40*F32i*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F40*F31i*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F40*F30i*c0*f3*u0**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F35*F45i*c5*f4*u5**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F35*F44i*c4*f4*u4*
     & u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F35*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F35*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F35*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F35*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F34*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F34*F44i*c4*f4*u4**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F34*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F34*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F34*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F34*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F33*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F33*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F33*F43i*c3*f4*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F33*F42i*c2*f4*u2*
     & u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F33*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F33*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F32*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F32*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F32*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F32*F42i*c2*f4*u2**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F32*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F32*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F31*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F31*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F31*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F31*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F31*F41i*c1*f4*u1**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F31*F40i*c0*f4*u0*
     & u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F30*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F30*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F30*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F30*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F30*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svm*dssp*F30*F40i*c0*f4*u0**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F45i*c5*f4*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F42i*c2*f4*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F42i*c2*f4*u2*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F45*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F41i*c1*f4*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F41i*c1*f4*u1*
     & u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F44*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F43i*c3*f4*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F40i*c0*f4*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F43*F40i*c0*f4*u0*
     & u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F42*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F45i*c5*f4*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F45i*c5*f4*u1*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F41i*c1*f4*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F41*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F44i*c4*f4*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F44i*c4*f4*u0*
     & u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F40i*c0*f4*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F40*F40i*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F25i*c5*f2*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F25i*c5*f2*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F24i*c4*f2*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F23i*c3*f2*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F23i*c3*f2*u3*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F35*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F25i*c5*f2*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F24i*c4*f2*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F24i*c4*f2*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F23i*c3*f2*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F22i*c2*f2*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F22i*c2*f2*u2*
     & u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F21i*c1*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F34*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F24i*c4*f2*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F23i*c3*f2*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F23i*c3*f2*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F22i*c2*f2*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F21i*c1*f2*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F21i*c1*f2*u1*
     & u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F33*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F23i*c3*f2*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F23i*c3*f2*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F22i*c2*f2*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F22i*c2*f2*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F21i*c1*f2*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F21i*c1*f2*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F20i*c0*f2*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F32*F20i*c0*f2*u0*
     & u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F24i*c4*f2*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F22i*c2*f2*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F21i*c1*f2*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F21i*c1*f2*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F20i*c0*f2*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F31*F20i*c0*f2*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F25i*c5*f2*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F25i*c5*f2*u0*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F21i*c1*f2*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F20i*c0*f2*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F30*F20i*c0*f2*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F34i*c4*f3*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F34i*c4*f3*u4*
     & u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F30i*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F25*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F34i*c4*f3*u4**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F33i*c3*f3*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F33i*c3*f3*u3*
     & u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F24*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F35i*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F35i*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F32i*c2*f3*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F32i*c2*f3*u2*
     & u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F23*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F34i*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F32i*c2*f3*u2**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F31i*c1*f3*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F31i*c1*f3*u1*
     & u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F22*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F33i*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F33i*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F30i*c0*f3*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F21*F30i*c0*f3*u0*
     & u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F32i*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F30i*c0*f3*u0**2*w2
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dsvm*F20*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F45*F35i*c5*f3*u5**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F45*F34i*c4*f3*u4*
     & u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F45*F33i*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F45*F32i*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F45*F31i*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F45*F30i*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F44*F35i*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F44*F34i*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F44*F33i*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F44*F32i*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F44*F31i*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F44*F30i*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F43*F35i*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F43*F34i*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F43*F33i*c3*f3*u3**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F43*F32i*c2*f3*u2*
     & u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F43*F31i*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F43*F30i*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F42*F35i*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F42*F34i*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F42*F33i*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F42*F32i*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F42*F31i*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F42*F30i*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F41*F35i*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F41*F34i*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F41*F33i*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F41*F32i*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F41*F31i*c1*f3*u1**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F41*F30i*c0*f3*u0*
     & u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F40*F35i*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F40*F34i*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F40*F33i*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F40*F32i*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F40*F31i*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F40*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F35*F45i*c5*f4*u5**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F35*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F35*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F35*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F35*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F35*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F34*F45i*c5*f4*u4*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F34*F44i*c4*f4*
     & u4**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F34*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F34*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F34*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F34*F40i*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F33*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F33*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F33*F43i*c3*f4*u3**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F33*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F33*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F33*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F32*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F32*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F32*F43i*c3*f4*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F32*F42i*c2*f4*
     & u2**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F32*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F32*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F31*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F31*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F31*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F31*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F31*F41i*c1*f4*u1**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F31*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F30*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F30*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F30*F43i*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F30*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F30*F41i*c1*f4*u0*u1*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*svp*dssm*F30*F40i*c0*f4*
     & u0**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F45*F35i*c5*f3*u5**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F45*F34i*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F45*F33i*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F45*F32i*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F45*F31i*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F45*F30i*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F44*F35i*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F44*F34i*c4*f3*u4**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F44*F33i*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F44*F32i*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F44*F31i*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F44*F30i*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F43*F35i*c5*f3*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F43*F34i*c4*f3*u3*
     & u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F43*F33i*c3*f3*u3**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F43*F32i*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F43*F31i*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F43*F30i*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F42*F35i*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F42*F34i*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F42*F33i*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F42*F32i*c2*f3*u2**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F42*F31i*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F42*F30i*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F41*F35i*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F41*F34i*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F41*F33i*c3*f3*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F41*F32i*c2*f3*u1*
     & u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F41*F31i*c1*f3*u1**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F41*F30i*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F40*F35i*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F40*F34i*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F40*F33i*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F40*F32i*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F40*F31i*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F40*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F35*F45i*c5*f4*u5**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F35*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F35*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F35*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F35*F41i*c1*f4*u1*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F35*F40i*c0*f4*u0*
     & u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F34*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F34*F44i*c4*f4*u4**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F34*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F34*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F34*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F34*F40i*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F33*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F33*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F33*F43i*c3*f4*u3**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F33*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F33*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F33*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F32*F45i*c5*f4*u2*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F32*F44i*c4*f4*u2*
     & u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F32*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F32*F42i*c2*f4*u2**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F32*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F32*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F31*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F31*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F31*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F31*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F31*F41i*c1*f4*u1**2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F31*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F30*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F30*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F30*F43i*c3*f4*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F30*F42i*c2*f4*u0*
     & u2*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F30*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q*i_*PC2*qq*ssm*dsvp*F30*F40i*c0*f4*u0**2*w1
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F45*F35i*c5*f3*u5**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F45*F34i*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F45*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F45*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F45*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F45*F30i*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F44*F35i*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F44*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F44*F33i*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F44*F32i*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F44*F31i*c1*f3*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F44*F30i*c0*f3*u0*
     & u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F43*F35i*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F43*F34i*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F43*F33i*c3*f3*u3**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F43*F32i*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F43*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F43*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F42*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F42*F34i*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F42*F33i*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F42*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F42*F31i*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F42*F30i*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F41*F35i*c5*f3*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F41*F34i*c4*f3*u1*
     & u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F41*F33i*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F41*F32i*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F41*F31i*c1*f3*u1**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F41*F30i*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F40*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F40*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F40*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F40*F32i*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F40*F31i*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F40*F30i*c0*f3*u0**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F35*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F35*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F35*F43i*c3*f4*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F35*F42i*c2*f4*u2*
     & u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F35*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F35*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F34*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F34*F44i*c4*f4*u4**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F34*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F34*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F34*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F34*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F33*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F33*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F33*F43i*c3*f4*u3**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F33*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F33*F41i*c1*f4*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F33*F40i*c0*f4*u0*
     & u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F32*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F32*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F32*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F32*F42i*c2*f4*u2**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F32*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F32*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F31*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F31*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F31*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F31*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F31*F41i*c1*f4*u1**2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F31*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F30*F45i*c5*f4*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F30*F44i*c4*f4*u0*
     & u4*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F30*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F30*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F30*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q*i_*PC2*qq*ssp*dsvm*F30*F40i*c0*f4*u0**2*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F45*F45i*c5*f4*u5**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F45*F44i*c4*f4*u4*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F45*F43i*c3*f4*u3*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F45*F42i*c2*f4*u2*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F45*F41i*c1*f4*u1*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F45*F40i*c0*f4*u0*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F44*F45i*c5*f4*u4*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F44*F44i*c4*f4*u4**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F44*F43i*c3*f4*u3*u4*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**2*PC_uh*svp*svm*F44*F42i*c2*f4*u2*u4
     & *w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F44*F41i*c1*f4*u1*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F44*F40i*c0*f4*u0*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F43*F45i*c5*f4*u3*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F43*F44i*c4*f4*u3*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F43*F43i*c3*f4*u3**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F43*F42i*c2*f4*u2*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F43*F41i*c1*f4*u1*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F43*F40i*c0*f4*u0*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F42*F45i*c5*f4*u2*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F42*F44i*c4*f4*u2*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F42*F43i*c3*f4*u2*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F42*F42i*c2*f4*u2**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F42*F41i*c1*f4*u1*u2*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**2*PC_uh*svp*svm*F42*F40i*c0*f4*u0*u2
     & *w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F41*F45i*c5*f4*u1*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F41*F44i*c4*f4*u1*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F41*F43i*c3*f4*u1*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F41*F42i*c2*f4*u1*u2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F41*F41i*c1*f4*u1**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F41*F40i*c0*f4*u0*u1*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F40*F45i*c5*f4*u0*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F40*F44i*c4*f4*u0*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F40*F43i*c3*f4*u0*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F40*F42i*c2*f4*u0*u2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F40*F41i*c1*f4*u0*u1*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*svp*svm*F40*F40i*c0*f4*u0**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F35*F35i*c5*f3*u5**2*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F35*F34i*c4*f3*u4
     & *u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F35*F33i*c3*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F35*F32i*c2*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F35*F31i*c1*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F35*F30i*c0*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F34*F35i*c5*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F34*F34i*c4*f3*u4**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F34*F33i*c3*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F34*F32i*c2*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F34*F31i*c1*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F34*F30i*c0*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F33*F35i*c5*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F33*F34i*c4*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F33*F33i*c3*f3*u3**2*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F33*F32i*c2*f3*u2
     & *u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F33*F31i*c1*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F33*F30i*c0*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F32*F35i*c5*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F32*F34i*c4*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F32*F33i*c3*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F32*F32i*c2*f3*u2**2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F32*F31i*c1*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F32*F30i*c0*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F31*F35i*c5*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F31*F34i*c4*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F31*F33i*c3*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F31*F32i*c2*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F31*F31i*c1*f3*u1**2*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F31*F30i*c0*f3*u0
     & *u1*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F30*F35i*c5*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F30*F34i*c4*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F30*F33i*c3*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F30*F32i*c2*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F30*F31i*c1*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**2*PC_uh*qq*svp*svm*F30*F30i*c0*f3*u0**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F45*F45r*c5*f4*u5**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F45*F44r*c4*f4*u4*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F45*F43r*c3*f4*u3*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F45*F42r*c2*f4*u2*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F45*F41r*c1*f4*u1*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F45*F40r*c0*f4*u0*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F44*F45r*c5*f4*u4*u5*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F44*F44r*c4*f4*
     & u4**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F44*F43r*c3*f4*u3*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F44*F42r*c2*f4*u2*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F44*F41r*c1*f4*u1*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F44*F40r*c0*f4*u0*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F43*F45r*c5*f4*u3*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F43*F44r*c4*f4*u3*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F43*F43r*c3*f4*u3**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F43*F42r*c2*f4*u2*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F43*F41r*c1*f4*u1*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F43*F40r*c0*f4*u0*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F42*F45r*c5*f4*u2*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F42*F44r*c4*f4*u2*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F42*F43r*c3*f4*u2*u3*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F42*F42r*c2*f4*
     & u2**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F42*F41r*c1*f4*u1*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F42*F40r*c0*f4*u0*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F41*F45r*c5*f4*u1*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F41*F44r*c4*f4*u1*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F41*F43r*c3*f4*u1*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F41*F42r*c2*f4*u1*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F41*F41r*c1*f4*u1**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F41*F40r*c0*f4*u0*u1*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F40*F45r*c5*f4*u0*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F40*F44r*c4*f4*u0*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F40*F43r*c3*f4*u0*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F40*F42r*c2*f4*u0*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F40*F41r*c1*f4*u0*u1*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q**2*PC_uh*i_*svp*svm*F40*F40r*c0*f4*
     & u0**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F35*F35r*c5*f3*u5**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F35*F34r*c4*f3*u4*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F35*F33r*c3*f3*u3*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F35*F32r*c2*f3*u2*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F35*F31r*c1*f3*u1*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F35*F30r*c0*f3*u0*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F34*F35r*c5*f3*u4*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F34*F34r*c4*f3*u4**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F34*F33r*c3*f3*u3*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F34*F32r*c2*f3*u2*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F34*F31r*c1*f3*u1*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F34*F30r*c0*f3*u0*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F33*F35r*c5*f3*u3*u5*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F33*F34r*c4*f3
     & *u3*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F33*F33r*c3*f3*u3**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F33*F32r*c2*f3*u2*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F33*F31r*c1*f3*u1*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F33*F30r*c0*f3*u0*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F32*F35r*c5*f3*u2*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F32*F34r*c4*f3*u2*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F32*F33r*c3*f3*u2*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F32*F32r*c2*f3*u2**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F32*F31r*c1*f3*u1*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F32*F30r*c0*f3*u0*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F31*F35r*c5*f3*u1*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F31*F34r*c4*f3*u1*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F31*F33r*c3*f3*u1*u3*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F31*F32r*c2*f3
     & *u1*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F31*F31r*c1*f3*u1**2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F31*F30r*c0*f3*u0*u1*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F30*F35r*c5*f3*u0*u5*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F30*F34r*c4*f3*u0*u4*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F30*F33r*c3*f3*u0*u3*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F30*F32r*c2*f3*u0*u2*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F30*F31r*c1*f3*u0*u1*w1*w2
     &  + 8.D0*PC_q**2*PC_uh*i_*qq*svp*svm*F30*F30r*c0*f3*u0**2*w1*w2
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F45*F45i*c5*f4*u5**2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F45*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F45*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F45*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F45*F43i*c3*f4*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*svp*svm*F45*F43i*c3*f4*u3*u5*
     & w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F45*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F45*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F45*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F45*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F45*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F45*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F44*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F44*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F44*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F44*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F44*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F44*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F44*F42i*c2*f4*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*svp*svm*F44*F42i*c2*f4*u2*u4*
     & w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F44*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F44*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F44*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F44*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F43*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F43*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F43*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F43*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F43*F43i*c3*f4*u3**2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F43*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F43*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F43*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F43*F41i*c1*f4*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*svp*svm*F43*F41i*c1*f4*u1*u3*
     & w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F43*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F43*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F42*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F42*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F42*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F42*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F42*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F42*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F42*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F42*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F42*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F42*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F42*F40i*c0*f4*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*svp*svm*F42*F40i*c0*f4*u0*u2*
     & w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F41*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F41*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F41*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F41*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F41*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F41*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F41*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F41*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F41*F41i*c1*f4*u1**2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F41*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F41*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F41*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F40*F45i*c5*f4*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*svp*svm*F40*F45i*c5*f4*u0*u5*
     & w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F40*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F40*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F40*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F40*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F40*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F40*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F40*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F40*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*svp*svm*F40*F40i*c0*f4*u0**2*w2
     &  + 4.D0*PC_q**2*q_uh*svp*svm*F40*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F45*F35i*c5*f3*u5**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F45*F34i*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F45*F33i*c3*f3*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*ssm*svp*F45*F32i*c2*f3*u2*u5*
     & w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F45*F31i*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F45*F30i*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F44*F35i*c5*f3*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F44*F34i*c4*f3*u4**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F44*F33i*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F44*F32i*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F44*F31i*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F44*F30i*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F43*F35i*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F43*F34i*c4*f3*u3*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F43*F33i*c3*f3*u3**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F43*F32i*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F43*F31i*c1*f3*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*ssm*svp*F43*F30i*c0*f3*u0*u3*
     & w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F42*F35i*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F42*F34i*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F42*F33i*c3*f3*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F42*F32i*c2*f3*u2**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F42*F31i*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F42*F30i*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F41*F35i*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F41*F34i*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F41*F33i*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F41*F32i*c2*f3*u1*u2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F41*F31i*c1*f3*u1**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F41*F30i*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F40*F35i*c5*f3*u0*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*ssm*svp*F40*F34i*c4*f3*u0*u4*
     & w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F40*F33i*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F40*F32i*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F40*F31i*c1*f3*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F40*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F35*F45i*c5*f4*u5**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F35*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F35*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F35*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F35*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F35*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F34*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F34*F44i*c4*f4*u4**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F34*F43i*c3*f4*u3*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*ssm*svp*F34*F42i*c2*f4*u2*u4*
     & w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F34*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F34*F40i*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F33*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F33*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F33*F43i*c3*f4*u3**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F33*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F33*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F33*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F32*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F32*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F32*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F32*F42i*c2*f4*u2**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F32*F41i*c1*f4*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*ssm*svp*F32*F40i*c0*f4*u0*u2*
     & w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F31*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F31*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F31*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F31*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F31*F41i*c1*f4*u1**2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F31*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F30*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F30*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F30*F43i*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F30*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F30*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*ssm*svp*F30*F40i*c0*f4*u0**2*w1
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F45*F35i*c5*f3*u5**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*ssp*svm*F45*F34i*c4*f3*u4*u5*
     & w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F45*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F45*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F45*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F45*F30i*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F44*F35i*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F44*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F44*F33i*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F44*F32i*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F44*F31i*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F44*F30i*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F43*F35i*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F43*F34i*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F43*F33i*c3*f3*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*ssp*svm*F43*F32i*c2*f3*u2*u3*
     & w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F43*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F43*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F42*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F42*F34i*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F42*F33i*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F42*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F42*F31i*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F42*F30i*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F41*F35i*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F41*F34i*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F41*F33i*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F41*F32i*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F41*F31i*c1*f3*u1**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*ssp*svm*F41*F30i*c0*f3*u0*u1*
     & w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F40*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F40*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F40*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F40*F32i*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F40*F31i*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F40*F30i*c0*f3*u0**2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F35*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F35*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F35*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F35*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F35*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F35*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F34*F45i*c5*f4*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*ssp*svm*F34*F44i*c4*f4*u4**2*
     & w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F34*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F34*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F34*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F34*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F33*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F33*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F33*F43i*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F33*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F33*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F33*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F32*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F32*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F32*F43i*c3*f4*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*ssp*svm*F32*F42i*c2*f4*u2**2*
     & w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F32*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F32*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F31*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F31*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F31*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F31*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F31*F41i*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F31*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F30*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F30*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F30*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F30*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*ssp*svm*F30*F41i*c1*f4*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*ssp*svm*F30*F40i*c0*f4*u0**2*
     & w2
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F35i*c5*f3*u5**2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F35i*c5*f3*u5**2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F34i*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F34i*c4*f3*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F33i*c3*f3*u3*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F32i*c2*f3*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F31i*c1*f3*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F30i*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F35*F30i*c0*f3*u0*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F35i*c5*f3*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F35i*c5*f3*u4*
     & u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F34i*c4*f3*u4**2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F33i*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F33i*c3*f3*u3*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F32i*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F32i*c2*f3*u2*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F31i*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F31i*c1*f3*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F30i*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F34*F30i*c0*f3*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F35i*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F35i*c5*f3*u3*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F34i*c4*f3*u3*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F34i*c4*f3*u3*
     & u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F33i*c3*f3*u3**2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F33i*c3*f3*u3**2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F32i*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F32i*c2*f3*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F31i*c1*f3*u1*u3*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F33*F30i*c0*f3*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F35i*c5*f3*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F34i*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F34i*c4*f3*u2*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F33i*c3*f3*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F33i*c3*f3*u2*
     & u3*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F32i*c2*f3*u2**2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F31i*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F31i*c1*f3*u1*u2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F30i*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F32*F30i*c0*f3*u0*u2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F35i*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F35i*c5*f3*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F34i*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F34i*c4*f3*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F33i*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F33i*c3*f3*u1*u3*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F32i*c2*f3*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F32i*c2*f3*u1*
     & u2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F31i*c1*f3*u1**2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F31i*c1*f3*u1**2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F30i*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F31*F30i*c0*f3*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F35i*c5*f3*u0*u5*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F34i*c4*f3*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F33i*c3*f3*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F32i*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F32i*c2*f3*u0*u2*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F31i*c1*f3*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F31i*c1*f3*u0*
     & u1*w1
     &  + 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F30i*c0*f3*u0**2*w2
     &  - 4.D0*PC_q**2*q_uh*qq*svp*svm*F30*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F40r*c0*f4*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*svp*svm*F45*F40r*c0*f4*u0*
     & u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F44*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F45r*c5*f4*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F45r*c5*f4*u3*
     & u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F43*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F44r*c4*f4*u2*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F44r*c4*f4*u2*
     & u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F42*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F43r*c3*f4*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F43r*c3*f4*u1*
     & u3*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F42r*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F41*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F42r*c2*f4*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F42r*c2*f4*u0*
     & u2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*svp*svm*F40*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F45*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F45*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F45*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F45*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F45*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F45*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F44*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F44*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F44*F33r*c3*f3*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F44*F32r*c2*f3*u2*
     & u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F44*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F44*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F43*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F43*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F43*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F43*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F43*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F43*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F42*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F42*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F42*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F42*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F42*F31r*c1*f3*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F42*F30r*c0*f3*u0*
     & u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F41*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F41*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F41*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F41*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F41*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F41*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F40*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F40*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F40*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F40*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F40*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F40*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F35*F45r*c5*f4*u5**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F35*F44r*c4*f4*u4*
     & u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F35*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F35*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F35*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F35*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F34*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F34*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F34*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F34*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F34*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F34*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F33*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F33*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F33*F43r*c3*f4*u3**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F33*F42r*c2*f4*u2*
     & u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F33*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F33*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F32*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F32*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F32*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F32*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F32*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F32*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F31*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F31*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F31*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F31*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F31*F41r*c1*f4*u1**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F31*F40r*c0*f4*u0*
     & u1*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F30*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F30*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F30*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F30*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F30*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*i_*ssm*svp*F30*F40r*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F45*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F45*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F45*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F45*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F45*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F45*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F44*F35r*c5*f3*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F44*F34r*c4*f3*
     & u4**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F44*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F44*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F44*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F44*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F43*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F43*F34r*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F43*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F43*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F43*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F43*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F42*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F42*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F42*F33r*c3*f3*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F42*F32r*c2*f3*
     & u2**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F42*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F42*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F41*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F41*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F41*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F41*F32r*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F41*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F41*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F40*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F40*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F40*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F40*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F40*F31r*c1*f3*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F40*F30r*c0*f3*
     & u0**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F35*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F35*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F35*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F35*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F35*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F35*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F34*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F34*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F34*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F34*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F34*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F34*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F33*F45r*c5*f4*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F33*F44r*c4*f4*u3*
     & u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F33*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F33*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F33*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F33*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F32*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F32*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F32*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F32*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F32*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F32*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F31*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F31*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F31*F43r*c3*f4*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F31*F42r*c2*f4*u1*
     & u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F31*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F31*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F30*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F30*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F30*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F30*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F30*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*i_*ssp*svm*F30*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F33r*c3*f3*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F33r*c3*f3*
     & u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F35*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F35r*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F32r*c2*f3*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F32r*c2*f3*
     & u2*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F34*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F34r*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F31r*c1*f3*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F31r*c1*f3*
     & u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F33*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F33r*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F30r*c0*f3*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F32*F30r*c0*f3*
     & u0*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F32r*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F31*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F35r*c5*f3*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F35r*c5*f3*
     & u0*u5*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F31r*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**2*q_uh*i_*qq*svp*svm*F30*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F14r*c4*f1*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F45*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F45*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F13r*c3*f1*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F44*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F44*F10r*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F11r*c1*f1*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F43*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F43*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F42*F10r*c0*f1*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F41*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F41*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F14r*c4*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F40*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F40*F10r*c0*f1*u0**2*w1
     &  + 8.D0*PC_q**2*svm*dsvp*F25*F25r*c5*f2*u5**2
     &  + 8.D0*PC_q**2*svm*dsvp*F25*F24r*c4*f2*u4*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F25*F23r*c3*f2*u3*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F25*F22r*c2*f2*u2*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F25*F21r*c1*f2*u1*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F25*F20r*c0*f2*u0*u5
     &
      traza1 = traza1 + 8.D0*PC_q**2*svm*dsvp*F24*F25r*c5*f2*u4*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F24*F24r*c4*f2*u4**2
     &  + 8.D0*PC_q**2*svm*dsvp*F24*F23r*c3*f2*u3*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F24*F22r*c2*f2*u2*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F24*F21r*c1*f2*u1*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F24*F20r*c0*f2*u0*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F23*F25r*c5*f2*u3*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F23*F24r*c4*f2*u3*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F23*F23r*c3*f2*u3**2
     &  + 8.D0*PC_q**2*svm*dsvp*F23*F22r*c2*f2*u2*u3
     &  + 8.D0*PC_q**2*svm*dsvp*F23*F21r*c1*f2*u1*u3
     &  + 8.D0*PC_q**2*svm*dsvp*F23*F20r*c0*f2*u0*u3
     &  + 8.D0*PC_q**2*svm*dsvp*F22*F25r*c5*f2*u2*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F22*F24r*c4*f2*u2*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F22*F23r*c3*f2*u2*u3
     &
      traza1 = traza1 + 8.D0*PC_q**2*svm*dsvp*F22*F22r*c2*f2*u2**2
     &  + 8.D0*PC_q**2*svm*dsvp*F22*F21r*c1*f2*u1*u2
     &  + 8.D0*PC_q**2*svm*dsvp*F22*F20r*c0*f2*u0*u2
     &  + 8.D0*PC_q**2*svm*dsvp*F21*F25r*c5*f2*u1*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F21*F24r*c4*f2*u1*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F21*F23r*c3*f2*u1*u3
     &  + 8.D0*PC_q**2*svm*dsvp*F21*F22r*c2*f2*u1*u2
     &  + 8.D0*PC_q**2*svm*dsvp*F21*F21r*c1*f2*u1**2
     &  + 8.D0*PC_q**2*svm*dsvp*F21*F20r*c0*f2*u0*u1
     &  + 8.D0*PC_q**2*svm*dsvp*F20*F25r*c5*f2*u0*u5
     &  + 8.D0*PC_q**2*svm*dsvp*F20*F24r*c4*f2*u0*u4
     &  + 8.D0*PC_q**2*svm*dsvp*F20*F23r*c3*f2*u0*u3
     &  + 8.D0*PC_q**2*svm*dsvp*F20*F22r*c2*f2*u0*u2
     &  + 8.D0*PC_q**2*svm*dsvp*F20*F21r*c1*f2*u0*u1
     &  + 8.D0*PC_q**2*svm*dsvp*F20*F20r*c0*f2*u0**2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F15*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F45r*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F15*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F44r*c4*f4*u4**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F14*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F14*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F43r*c3*f4*u3**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F13*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F42r*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F13*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F41r*c1*f4*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F12*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F12*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F42r*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F11*F40r*c0*f4*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dsvp*F10*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**2*svm*dsvp*F10*F40r*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**2*svm*dssp*F45*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*svm*dssp*F45*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*svm*dssp*F45*F23r*c3*f2*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*svm*dssp*F45*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*svm*dssp*F45*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*svm*dssp*F45*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*svm*dssp*F44*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*svm*dssp*F44*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*svm*dssp*F44*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*svm*dssp*F44*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*svm*dssp*F44*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*svm*dssp*F44*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*svm*dssp*F43*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*svm*dssp*F43*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*svm*dssp*F43*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q**2*svm*dssp*F43*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*svm*dssp*F43*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*svm*dssp*F43*F20r*c0*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*svm*dssp*F42*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*svm*dssp*F42*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*svm*dssp*F42*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*svm*dssp*F42*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*svm*dssp*F42*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*svm*dssp*F42*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*svm*dssp*F41*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*svm*dssp*F41*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*svm*dssp*F41*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*svm*dssp*F41*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*svm*dssp*F41*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*svm*dssp*F41*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*svm*dssp*F40*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*svm*dssp*F40*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*svm*dssp*F40*F23r*c3*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*svm*dssp*F40*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*svm*dssp*F40*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*svm*dssp*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC_q**2*svm*dssp*F35*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F35*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F35*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F35*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F35*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F35*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F34*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F34*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F34*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F34*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F34*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F34*F10r*c0*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dssp*F33*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F33*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F33*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F33*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F33*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F33*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F32*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F32*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F32*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F32*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F32*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F32*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F31*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F31*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F31*F13r*c3*f1*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dssp*F31*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F31*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F31*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dssp*F30*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F30*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F30*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F30*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F30*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dssp*F30*F10r*c0*f1*u0**2*w2
     &  + 4.D0*PC_q**2*svm*dssp*F25*F45r*c5*f4*u5**2
     &  + 4.D0*PC_q**2*svm*dssp*F25*F44r*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*svm*dssp*F25*F43r*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*svm*dssp*F25*F42r*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*svm*dssp*F25*F41r*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*svm*dssp*F25*F40r*c0*f4*u0*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*svm*dssp*F24*F45r*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*svm*dssp*F24*F44r*c4*f4*u4**2
     &  + 4.D0*PC_q**2*svm*dssp*F24*F43r*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*svm*dssp*F24*F42r*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*svm*dssp*F24*F41r*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*svm*dssp*F24*F40r*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*svm*dssp*F23*F45r*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*svm*dssp*F23*F44r*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*svm*dssp*F23*F43r*c3*f4*u3**2
     &  + 4.D0*PC_q**2*svm*dssp*F23*F42r*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*svm*dssp*F23*F41r*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*svm*dssp*F23*F40r*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*svm*dssp*F22*F45r*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*svm*dssp*F22*F44r*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*svm*dssp*F22*F43r*c3*f4*u2*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*svm*dssp*F22*F42r*c2*f4*u2**2
     &  + 4.D0*PC_q**2*svm*dssp*F22*F41r*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*svm*dssp*F22*F40r*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*svm*dssp*F21*F45r*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*svm*dssp*F21*F44r*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*svm*dssp*F21*F43r*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*svm*dssp*F21*F42r*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*svm*dssp*F21*F41r*c1*f4*u1**2
     &  + 4.D0*PC_q**2*svm*dssp*F21*F40r*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*svm*dssp*F20*F45r*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*svm*dssp*F20*F44r*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*svm*dssp*F20*F43r*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*svm*dssp*F20*F42r*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*svm*dssp*F20*F41r*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*svm*dssp*F20*F40r*c0*f4*u0**2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dssp*F15*F35r*c5*f3*u5**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F15*F34r*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F15*F33r*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F15*F32r*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F15*F31r*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F15*F30r*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F14*F35r*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F14*F34r*c4*f3*u4**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F14*F33r*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F14*F32r*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F14*F31r*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F14*F30r*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F13*F35r*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F13*F34r*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F13*F33r*c3*f3*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dssp*F13*F32r*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F13*F31r*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F13*F30r*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F12*F35r*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F12*F34r*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F12*F33r*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F12*F32r*c2*f3*u2**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F12*F31r*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F12*F30r*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F11*F35r*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F11*F34r*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F11*F33r*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F11*F32r*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F11*F31r*c1*f3*u1**2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F11*F30r*c0*f3*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svm*dssp*F10*F35r*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*svm*dssp*F10*F34r*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*svm*dssp*F10*F33r*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*svm*dssp*F10*F32r*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*svm*dssp*F10*F31r*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*svm*dssp*F10*F30r*c0*f3*u0**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F11r*c1*f1*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F45*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F45*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F44*F10r*c0*f1*u0*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F43*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F43*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F14r*c4*f1*u2*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F42*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F42*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F13r*c3*f1*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F41*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F41*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F11r*c1*f1*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F40*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F40*F10r*c0*f1*u0**2*w1
     &  + 8.D0*PC_q**2*svp*dsvm*F25*F25r*c5*f2*u5**2
     &  + 8.D0*PC_q**2*svp*dsvm*F25*F24r*c4*f2*u4*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F25*F23r*c3*f2*u3*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F25*F22r*c2*f2*u2*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F25*F21r*c1*f2*u1*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F25*F20r*c0*f2*u0*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F24*F25r*c5*f2*u4*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F24*F24r*c4*f2*u4**2
     &  + 8.D0*PC_q**2*svp*dsvm*F24*F23r*c3*f2*u3*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F24*F22r*c2*f2*u2*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F24*F21r*c1*f2*u1*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F24*F20r*c0*f2*u0*u4
     &
      traza1 = traza1 + 8.D0*PC_q**2*svp*dsvm*F23*F25r*c5*f2*u3*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F23*F24r*c4*f2*u3*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F23*F23r*c3*f2*u3**2
     &  + 8.D0*PC_q**2*svp*dsvm*F23*F22r*c2*f2*u2*u3
     &  + 8.D0*PC_q**2*svp*dsvm*F23*F21r*c1*f2*u1*u3
     &  + 8.D0*PC_q**2*svp*dsvm*F23*F20r*c0*f2*u0*u3
     &  + 8.D0*PC_q**2*svp*dsvm*F22*F25r*c5*f2*u2*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F22*F24r*c4*f2*u2*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F22*F23r*c3*f2*u2*u3
     &  + 8.D0*PC_q**2*svp*dsvm*F22*F22r*c2*f2*u2**2
     &  + 8.D0*PC_q**2*svp*dsvm*F22*F21r*c1*f2*u1*u2
     &  + 8.D0*PC_q**2*svp*dsvm*F22*F20r*c0*f2*u0*u2
     &  + 8.D0*PC_q**2*svp*dsvm*F21*F25r*c5*f2*u1*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F21*F24r*c4*f2*u1*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F21*F23r*c3*f2*u1*u3
     &
      traza1 = traza1 + 8.D0*PC_q**2*svp*dsvm*F21*F22r*c2*f2*u1*u2
     &  + 8.D0*PC_q**2*svp*dsvm*F21*F21r*c1*f2*u1**2
     &  + 8.D0*PC_q**2*svp*dsvm*F21*F20r*c0*f2*u0*u1
     &  + 8.D0*PC_q**2*svp*dsvm*F20*F25r*c5*f2*u0*u5
     &  + 8.D0*PC_q**2*svp*dsvm*F20*F24r*c4*f2*u0*u4
     &  + 8.D0*PC_q**2*svp*dsvm*F20*F23r*c3*f2*u0*u3
     &  + 8.D0*PC_q**2*svp*dsvm*F20*F22r*c2*f2*u0*u2
     &  + 8.D0*PC_q**2*svp*dsvm*F20*F21r*c1*f2*u0*u1
     &  + 8.D0*PC_q**2*svp*dsvm*F20*F20r*c0*f2*u0**2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F45r*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F43r*c3*f4*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F15*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F15*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F41r*c1*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F14*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F14*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F43r*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F42r*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F13*F40r*c0*f4*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F12*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F12*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F44r*c4*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F11*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F42r*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F11*F40r*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F43r*c3*f4*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dsvm*F10*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**2*svp*dsvm*F10*F40r*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**2*svp*dssm*F45*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*svp*dssm*F45*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*svp*dssm*F45*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*svp*dssm*F45*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*svp*dssm*F45*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*svp*dssm*F45*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*svp*dssm*F44*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*svp*dssm*F44*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*svp*dssm*F44*F23r*c3*f2*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*svp*dssm*F44*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*svp*dssm*F44*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*svp*dssm*F44*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*svp*dssm*F43*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*svp*dssm*F43*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*svp*dssm*F43*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q**2*svp*dssm*F43*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*svp*dssm*F43*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*svp*dssm*F43*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*svp*dssm*F42*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*svp*dssm*F42*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*svp*dssm*F42*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*svp*dssm*F42*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*svp*dssm*F42*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*svp*dssm*F42*F20r*c0*f2*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*svp*dssm*F41*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*svp*dssm*F41*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*svp*dssm*F41*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*svp*dssm*F41*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*svp*dssm*F41*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*svp*dssm*F41*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*svp*dssm*F40*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*svp*dssm*F40*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*svp*dssm*F40*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*svp*dssm*F40*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*svp*dssm*F40*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*svp*dssm*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC_q**2*svp*dssm*F35*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F35*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F35*F13r*c3*f1*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dssm*F35*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F35*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F35*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F34*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F34*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F34*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F34*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F34*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F34*F10r*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F33*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F33*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F33*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F33*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F33*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F33*F10r*c0*f1*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dssm*F32*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F32*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F32*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F32*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F32*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F32*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F31*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F31*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F31*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F31*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F31*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F31*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dssm*F30*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F30*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F30*F13r*c3*f1*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dssm*F30*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F30*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dssm*F30*F10r*c0*f1*u0**2*w1
     &  + 4.D0*PC_q**2*svp*dssm*F25*F45r*c5*f4*u5**2
     &  + 4.D0*PC_q**2*svp*dssm*F25*F44r*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*svp*dssm*F25*F43r*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*svp*dssm*F25*F42r*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*svp*dssm*F25*F41r*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*svp*dssm*F25*F40r*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*svp*dssm*F24*F45r*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*svp*dssm*F24*F44r*c4*f4*u4**2
     &  + 4.D0*PC_q**2*svp*dssm*F24*F43r*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*svp*dssm*F24*F42r*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*svp*dssm*F24*F41r*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*svp*dssm*F24*F40r*c0*f4*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*svp*dssm*F23*F45r*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*svp*dssm*F23*F44r*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*svp*dssm*F23*F43r*c3*f4*u3**2
     &  + 4.D0*PC_q**2*svp*dssm*F23*F42r*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*svp*dssm*F23*F41r*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*svp*dssm*F23*F40r*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*svp*dssm*F22*F45r*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*svp*dssm*F22*F44r*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*svp*dssm*F22*F43r*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*svp*dssm*F22*F42r*c2*f4*u2**2
     &  + 4.D0*PC_q**2*svp*dssm*F22*F41r*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*svp*dssm*F22*F40r*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*svp*dssm*F21*F45r*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*svp*dssm*F21*F44r*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*svp*dssm*F21*F43r*c3*f4*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*svp*dssm*F21*F42r*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*svp*dssm*F21*F41r*c1*f4*u1**2
     &  + 4.D0*PC_q**2*svp*dssm*F21*F40r*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*svp*dssm*F20*F45r*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*svp*dssm*F20*F44r*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*svp*dssm*F20*F43r*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*svp*dssm*F20*F42r*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*svp*dssm*F20*F41r*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*svp*dssm*F20*F40r*c0*f4*u0**2
     &  - 4.D0*PC_q**2*svp*dssm*F15*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F15*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F15*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F15*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F15*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F15*F30r*c0*f3*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dssm*F14*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F14*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F14*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F14*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F14*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F14*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F13*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F13*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F13*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F13*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F13*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F13*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F12*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F12*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F12*F33r*c3*f3*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*svp*dssm*F12*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F12*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F12*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F11*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F11*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F11*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F11*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F11*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F11*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dssm*F10*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*svp*dssm*F10*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*svp*dssm*F10*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*svp*dssm*F10*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*svp*dssm*F10*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*svp*dssm*F10*F30r*c0*f3*u0**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dsvp*F45*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F45*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F45*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F45*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F45*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F45*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F44*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F44*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F44*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F44*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F44*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F44*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F43*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F43*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F43*F23r*c3*f2*u3**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dsvp*F43*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F43*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F43*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F42*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F42*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F42*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F42*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F42*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F42*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F41*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F41*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F41*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F41*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F41*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F41*F20r*c0*f2*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dsvp*F40*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F40*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F40*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F40*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F40*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*ssm*dsvp*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC_q**2*ssm*dsvp*F35*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F35*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F35*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F35*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F35*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F35*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F34*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F34*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F34*F13r*c3*f1*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*dsvp*F34*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F34*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F34*F10r*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F33*F15r*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F33*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F33*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F33*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F33*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F33*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F32*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F32*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F32*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F32*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F32*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F32*F10r*c0*f1*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*dsvp*F31*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F31*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F31*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F31*F12r*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F31*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F31*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F30*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F30*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F30*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F30*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F30*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F30*F10r*c0*f1*u0**2*w1
     &  + 4.D0*PC_q**2*ssm*dsvp*F25*F45r*c5*f4*u5**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F25*F44r*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F25*F43r*c3*f4*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dsvp*F25*F42r*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F25*F41r*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F25*F40r*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F24*F45r*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F24*F44r*c4*f4*u4**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F24*F43r*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F24*F42r*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F24*F41r*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F24*F40r*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F23*F45r*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F23*F44r*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F23*F43r*c3*f4*u3**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F23*F42r*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F23*F41r*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F23*F40r*c0*f4*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dsvp*F22*F45r*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F22*F44r*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F22*F43r*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F22*F42r*c2*f4*u2**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F22*F41r*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F22*F40r*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F21*F45r*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F21*F44r*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F21*F43r*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*ssm*dsvp*F21*F42r*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F21*F41r*c1*f4*u1**2
     &  + 4.D0*PC_q**2*ssm*dsvp*F21*F40r*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*ssm*dsvp*F20*F45r*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*ssm*dsvp*F20*F44r*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*ssm*dsvp*F20*F43r*c3*f4*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dsvp*F20*F42r*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*ssm*dsvp*F20*F41r*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*ssm*dsvp*F20*F40r*c0*f4*u0**2
     &  - 4.D0*PC_q**2*ssm*dsvp*F15*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F15*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F15*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F15*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F15*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F15*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F14*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F14*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F14*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F14*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F14*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F14*F30r*c0*f3*u0*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*dsvp*F13*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F13*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F13*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F13*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F13*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F13*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F12*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F12*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F12*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F12*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F12*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F12*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F11*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F11*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F11*F33r*c3*f3*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssm*dsvp*F11*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F11*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F11*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F10*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F10*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F10*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F10*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F10*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*ssm*dsvp*F10*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q**2*ssm*dssp*F45*F45r*c5*f4*u5**2
     &  + 4.D0*PC_q**2*ssm*dssp*F45*F44r*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F45*F43r*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F45*F42r*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F45*F41r*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F45*F40r*c0*f4*u0*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F44*F45r*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F44*F44r*c4*f4*u4**2
     &  + 4.D0*PC_q**2*ssm*dssp*F44*F43r*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F44*F42r*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F44*F41r*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F44*F40r*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F43*F45r*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F43*F44r*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F43*F43r*c3*f4*u3**2
     &  + 4.D0*PC_q**2*ssm*dssp*F43*F42r*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F43*F41r*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F43*F40r*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F42*F45r*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F42*F44r*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F42*F43r*c3*f4*u2*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F42*F42r*c2*f4*u2**2
     &  + 4.D0*PC_q**2*ssm*dssp*F42*F41r*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F42*F40r*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F41*F45r*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F41*F44r*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F41*F43r*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F41*F42r*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F41*F41r*c1*f4*u1**2
     &  + 4.D0*PC_q**2*ssm*dssp*F41*F40r*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*ssm*dssp*F40*F45r*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F40*F44r*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F40*F43r*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F40*F42r*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F40*F41r*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*ssm*dssp*F40*F40r*c0*f4*u0**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F35*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*ssm*dssp*F35*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F35*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F35*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F35*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F35*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F34*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F34*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*ssm*dssp*F34*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F34*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F34*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F34*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F33*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F33*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F33*F23r*c3*f2*u3**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F33*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F33*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F33*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F32*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F32*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F32*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F32*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*ssm*dssp*F32*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F32*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F31*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F31*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F31*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F31*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F31*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*ssm*dssp*F31*F20r*c0*f2*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F30*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F30*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F30*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F30*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F30*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*ssm*dssp*F30*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q**2*ssm*dssp*F25*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*ssm*dssp*F25*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F25*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F25*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F25*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F25*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F24*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F24*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*ssm*dssp*F24*F33r*c3*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F24*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F24*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F24*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F23*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F23*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F23*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*ssm*dssp*F23*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F23*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F23*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F22*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F22*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F22*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F22*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*ssm*dssp*F22*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F22*F30r*c0*f3*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssm*dssp*F21*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F21*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F21*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F21*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F21*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*ssm*dssp*F21*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*ssm*dssp*F20*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*ssm*dssp*F20*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*ssm*dssp*F20*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*ssm*dssp*F20*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*ssm*dssp*F20*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*ssm*dssp*F20*F30r*c0*f3*u0**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F45*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F45*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F45*F23r*c3*f2*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dsvm*F45*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F45*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F45*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F44*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F44*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F44*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F44*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F44*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F44*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F43*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F43*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F43*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F43*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F43*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F43*F20r*c0*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dsvm*F42*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F42*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F42*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F42*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F42*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F42*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F41*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F41*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F41*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F41*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F41*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F41*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*ssp*dsvm*F40*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F40*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F40*F23r*c3*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dsvm*F40*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F40*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*ssp*dsvm*F40*F20r*c0*f2*u0**2
     &  - 4.D0*PC_q**2*ssp*dsvm*F35*F15r*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F35*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F35*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F35*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F35*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F35*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F34*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F34*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F34*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F34*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F34*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F34*F10r*c0*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*dsvm*F33*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F33*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F33*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F33*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F33*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F33*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F32*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F32*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F32*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F32*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F32*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F32*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F31*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F31*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F31*F13r*c3*f1*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*dsvm*F31*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F31*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F31*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F30*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F30*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F30*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F30*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F30*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F30*F10r*c0*f1*u0**2*w2
     &  + 4.D0*PC_q**2*ssp*dsvm*F25*F45r*c5*f4*u5**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F25*F44r*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F25*F43r*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F25*F42r*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F25*F41r*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F25*F40r*c0*f4*u0*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dsvm*F24*F45r*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F24*F44r*c4*f4*u4**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F24*F43r*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F24*F42r*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F24*F41r*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F24*F40r*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F23*F45r*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F23*F44r*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F23*F43r*c3*f4*u3**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F23*F42r*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F23*F41r*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F23*F40r*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F22*F45r*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F22*F44r*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F22*F43r*c3*f4*u2*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dsvm*F22*F42r*c2*f4*u2**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F22*F41r*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F22*F40r*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F21*F45r*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F21*F44r*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F21*F43r*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F21*F42r*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F21*F41r*c1*f4*u1**2
     &  + 4.D0*PC_q**2*ssp*dsvm*F21*F40r*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*ssp*dsvm*F20*F45r*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*ssp*dsvm*F20*F44r*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*ssp*dsvm*F20*F43r*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*ssp*dsvm*F20*F42r*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*ssp*dsvm*F20*F41r*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*ssp*dsvm*F20*F40r*c0*f4*u0**2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*dsvm*F15*F35r*c5*f3*u5**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F15*F34r*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F15*F33r*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F15*F32r*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F15*F31r*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F15*F30r*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F14*F35r*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F14*F34r*c4*f3*u4**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F14*F33r*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F14*F32r*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F14*F31r*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F14*F30r*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F13*F35r*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F13*F34r*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F13*F33r*c3*f3*u3**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*dsvm*F13*F32r*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F13*F31r*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F13*F30r*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F12*F35r*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F12*F34r*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F12*F33r*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F12*F32r*c2*f3*u2**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F12*F31r*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F12*F30r*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F11*F35r*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F11*F34r*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F11*F33r*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F11*F32r*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F11*F31r*c1*f3*u1**2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F11*F30r*c0*f3*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*ssp*dsvm*F10*F35r*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F10*F34r*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F10*F33r*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F10*F32r*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F10*F31r*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*ssp*dsvm*F10*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**2*ssp*dssm*F45*F45r*c5*f4*u5**2
     &  + 4.D0*PC_q**2*ssp*dssm*F45*F44r*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F45*F43r*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F45*F42r*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F45*F41r*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F45*F40r*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F44*F45r*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F44*F44r*c4*f4*u4**2
     &  + 4.D0*PC_q**2*ssp*dssm*F44*F43r*c3*f4*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F44*F42r*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F44*F41r*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F44*F40r*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F43*F45r*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F43*F44r*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F43*F43r*c3*f4*u3**2
     &  + 4.D0*PC_q**2*ssp*dssm*F43*F42r*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F43*F41r*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F43*F40r*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F42*F45r*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F42*F44r*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F42*F43r*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F42*F42r*c2*f4*u2**2
     &  + 4.D0*PC_q**2*ssp*dssm*F42*F41r*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F42*F40r*c0*f4*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F41*F45r*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F41*F44r*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F41*F43r*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F41*F42r*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F41*F41r*c1*f4*u1**2
     &  + 4.D0*PC_q**2*ssp*dssm*F41*F40r*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*ssp*dssm*F40*F45r*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F40*F44r*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F40*F43r*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F40*F42r*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F40*F41r*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*ssp*dssm*F40*F40r*c0*f4*u0**2
     &  + 4.D0*PC_q**2*ssp*dssm*F35*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*ssp*dssm*F35*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F35*F23r*c3*f2*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F35*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F35*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F35*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F34*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F34*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*ssp*dssm*F34*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F34*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F34*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F34*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F33*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F33*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F33*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q**2*ssp*dssm*F33*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F33*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F33*F20r*c0*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F32*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F32*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F32*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F32*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*ssp*dssm*F32*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F32*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F31*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F31*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F31*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F31*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F31*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*ssp*dssm*F31*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*ssp*dssm*F30*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F30*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F30*F23r*c3*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F30*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F30*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*ssp*dssm*F30*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q**2*ssp*dssm*F25*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*ssp*dssm*F25*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F25*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F25*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F25*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F25*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F24*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F24*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*ssp*dssm*F24*F33r*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F24*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F24*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F24*F30r*c0*f3*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F23*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F23*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F23*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*ssp*dssm*F23*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F23*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F23*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F22*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F22*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F22*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F22*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*ssp*dssm*F22*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F22*F30r*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F21*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F21*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F21*F33r*c3*f3*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*ssp*dssm*F21*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F21*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*ssp*dssm*F21*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*ssp*dssm*F20*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*ssp*dssm*F20*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*ssp*dssm*F20*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*ssp*dssm*F20*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*ssp*dssm*F20*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*ssp*dssm*F20*F30r*c0*f3*u0**2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F45*F45r*c5*f4*u5**2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F45*F44r*c4*f4*u4*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F45*F43r*c3*f4*u3*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F45*F42r*c2*f4*u2*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F45*F41r*c1*f4*u1*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F45*F40r*c0*f4*u0*u5
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svm*dsvp*F44*F45r*c5*f4*u4*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F44*F44r*c4*f4*u4**2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F44*F43r*c3*f4*u3*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F44*F42r*c2*f4*u2*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F44*F41r*c1*f4*u1*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F44*F40r*c0*f4*u0*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F43*F45r*c5*f4*u3*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F43*F44r*c4*f4*u3*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F43*F43r*c3*f4*u3**2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F43*F42r*c2*f4*u2*u3
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F43*F41r*c1*f4*u1*u3
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F43*F40r*c0*f4*u0*u3
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F42*F45r*c5*f4*u2*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F42*F44r*c4*f4*u2*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F42*F43r*c3*f4*u2*u3
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svm*dsvp*F42*F42r*c2*f4*u2**2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F42*F41r*c1*f4*u1*u2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F42*F40r*c0*f4*u0*u2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F41*F45r*c5*f4*u1*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F41*F44r*c4*f4*u1*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F41*F43r*c3*f4*u1*u3
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F41*F42r*c2*f4*u1*u2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F41*F41r*c1*f4*u1**2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F41*F40r*c0*f4*u0*u1
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F40*F45r*c5*f4*u0*u5
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F40*F44r*c4*f4*u0*u4
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F40*F43r*c3*f4*u0*u3
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F40*F42r*c2*f4*u0*u2
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F40*F41r*c1*f4*u0*u1
     &  - 4.D0*PC_q**2*qq*svm*dsvp*F40*F40r*c0*f4*u0**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svm*dsvp*F35*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F35*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F35*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F35*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F35*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F35*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F34*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F34*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F34*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F34*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F34*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F34*F20r*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F33*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F33*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F33*F23r*c3*f2*u3**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svm*dsvp*F33*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F33*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F33*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F32*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F32*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F32*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F32*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F32*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F32*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F31*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F31*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F31*F23r*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F31*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F31*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F31*F20r*c0*f2*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svm*dsvp*F30*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F30*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F30*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F30*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F30*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F30*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F25*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F25*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F25*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F25*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F25*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F25*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F24*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F24*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F24*F33r*c3*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svm*dsvp*F24*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F24*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F24*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F23*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F23*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F23*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F23*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F23*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F23*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F22*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F22*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F22*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F22*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F22*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F22*F30r*c0*f3*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svm*dsvp*F21*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F21*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F21*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F21*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F21*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F21*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F20*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F20*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F20*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F20*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F20*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*svm*dsvp*F20*F30r*c0*f3*u0**2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F45*F45r*c5*f4*u5**2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F45*F44r*c4*f4*u4*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F45*F43r*c3*f4*u3*u5
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svp*dsvm*F45*F42r*c2*f4*u2*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F45*F41r*c1*f4*u1*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F45*F40r*c0*f4*u0*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F44*F45r*c5*f4*u4*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F44*F44r*c4*f4*u4**2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F44*F43r*c3*f4*u3*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F44*F42r*c2*f4*u2*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F44*F41r*c1*f4*u1*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F44*F40r*c0*f4*u0*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F43*F45r*c5*f4*u3*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F43*F44r*c4*f4*u3*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F43*F43r*c3*f4*u3**2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F43*F42r*c2*f4*u2*u3
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F43*F41r*c1*f4*u1*u3
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F43*F40r*c0*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svp*dsvm*F42*F45r*c5*f4*u2*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F42*F44r*c4*f4*u2*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F42*F43r*c3*f4*u2*u3
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F42*F42r*c2*f4*u2**2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F42*F41r*c1*f4*u1*u2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F42*F40r*c0*f4*u0*u2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F41*F45r*c5*f4*u1*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F41*F44r*c4*f4*u1*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F41*F43r*c3*f4*u1*u3
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F41*F42r*c2*f4*u1*u2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F41*F41r*c1*f4*u1**2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F41*F40r*c0*f4*u0*u1
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F40*F45r*c5*f4*u0*u5
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F40*F44r*c4*f4*u0*u4
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F40*F43r*c3*f4*u0*u3
     &
      traza1 = traza1 - 4.D0*PC_q**2*qq*svp*dsvm*F40*F42r*c2*f4*u0*u2
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F40*F41r*c1*f4*u0*u1
     &  - 4.D0*PC_q**2*qq*svp*dsvm*F40*F40r*c0*f4*u0**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F35*F25r*c5*f2*u5**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F35*F24r*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F35*F23r*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F35*F22r*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F35*F21r*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F35*F20r*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F34*F25r*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F34*F24r*c4*f2*u4**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F34*F23r*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F34*F22r*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F34*F21r*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F34*F20r*c0*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*dsvm*F33*F25r*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F33*F24r*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F33*F23r*c3*f2*u3**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F33*F22r*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F33*F21r*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F33*F20r*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F32*F25r*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F32*F24r*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F32*F23r*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F32*F22r*c2*f2*u2**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F32*F21r*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F32*F20r*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F31*F25r*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F31*F24r*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F31*F23r*c3*f2*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*dsvm*F31*F22r*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F31*F21r*c1*f2*u1**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F31*F20r*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F30*F25r*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F30*F24r*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F30*F23r*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F30*F22r*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F30*F21r*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F30*F20r*c0*f2*u0**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F25*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F25*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F25*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F25*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F25*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F25*F30r*c0*f3*u0*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*dsvm*F24*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F24*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F24*F33r*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F24*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F24*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F24*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F23*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F23*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F23*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F23*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F23*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F23*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F22*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F22*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F22*F33r*c3*f3*u2*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*svp*dsvm*F22*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F22*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F22*F30r*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F21*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F21*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F21*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F21*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F21*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F21*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F20*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F20*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F20*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F20*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F20*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*svp*dsvm*F20*F30r*c0*f3*u0**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssm*dssp*F35*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F35*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F35*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F35*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F35*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F35*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F34*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F34*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F34*F33r*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F34*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F34*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F34*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F33*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F33*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F33*F33r*c3*f3*u3**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssm*dssp*F33*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F33*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F33*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F32*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F32*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F32*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F32*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F32*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F32*F30r*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F31*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F31*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F31*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F31*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F31*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F31*F30r*c0*f3*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssm*dssp*F30*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F30*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F30*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F30*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F30*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*ssm*dssp*F30*F30r*c0*f3*u0**2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F35*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F35*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F35*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F35*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F35*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F35*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F34*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F34*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F34*F33r*c3*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssp*dssm*F34*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F34*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F34*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F33*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F33*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F33*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F33*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F33*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F33*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F32*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F32*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F32*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F32*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F32*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F32*F30r*c0*f3*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq*ssp*dssm*F31*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F31*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F31*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F31*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F31*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F31*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F30*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F30*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F30*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F30*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F30*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*qq*ssp*dssm*F30*F30r*c0*f3*u0**2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F35*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F35*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F35*F33r*c3*f3*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svm*dsvp*F35*F32r*c2*f3*u2*
     & u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F35*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F35*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F34*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F34*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F34*F33r*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F34*F32r*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F34*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F34*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F33*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F33*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F33*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F33*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F33*F31r*c1*f3*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svm*dsvp*F33*F30r*c0*f3*u0*
     & u3
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F32*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F32*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F32*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F32*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F32*F31r*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F32*F30r*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F31*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F31*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F31*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F31*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F31*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F31*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F30*F35r*c5*f3*u0*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svm*dsvp*F30*F34r*c4*f3*u0*
     & u4
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F30*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F30*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F30*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*qq**2*svm*dsvp*F30*F30r*c0*f3*u0**2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F35*F35r*c5*f3*u5**2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F35*F34r*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F35*F33r*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F35*F32r*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F35*F31r*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F35*F30r*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F34*F35r*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F34*F34r*c4*f3*u4**2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F34*F33r*c3*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svp*dsvm*F34*F32r*c2*f3*u2*
     & u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F34*F31r*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F34*F30r*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F33*F35r*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F33*F34r*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F33*F33r*c3*f3*u3**2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F33*F32r*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F33*F31r*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F33*F30r*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F32*F35r*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F32*F34r*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F32*F33r*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F32*F32r*c2*f3*u2**2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F32*F31r*c1*f3*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*qq**2*svp*dsvm*F32*F30r*c0*f3*u0*
     & u2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F31*F35r*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F31*F34r*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F31*F33r*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F31*F32r*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F31*F31r*c1*f3*u1**2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F31*F30r*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F30*F35r*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F30*F34r*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F30*F33r*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F30*F32r*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F30*F31r*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*qq**2*svp*dsvm*F30*F30r*c0*f3*u0**2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F45*F45r*c5*f4*u5**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svm*dsvp*F45*F44r*c4*f4*u4*u5*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F45*F43r*c3*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F45*F42r*c2*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F45*F41r*c1*f4*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F45*F40r*c0*f4*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F44*F45r*c5*f4*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F44*F44r*c4*f4*u4**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F44*F43r*c3*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F44*F42r*c2*f4*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F44*F41r*c1*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F44*F40r*c0*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F43*F45r*c5*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F43*F44r*c4*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F43*F43r*c3*f4*u3**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svm*dsvp*F43*F42r*c2*f4*u2*u3*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F43*F41r*c1*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F43*F40r*c0*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F42*F45r*c5*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F42*F44r*c4*f4*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F42*F43r*c3*f4*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F42*F42r*c2*f4*u2**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F42*F41r*c1*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F42*F40r*c0*f4*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F41*F45r*c5*f4*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F41*F44r*c4*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F41*F43r*c3*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F41*F42r*c2*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F41*F41r*c1*f4*u1**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svm*dsvp*F41*F40r*c0*f4*u0*u1*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F40*F45r*c5*f4*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F40*F44r*c4*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F40*F43r*c3*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F40*F42r*c2*f4*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F40*F41r*c1*f4*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*svm*dsvp*F40*F40r*c0*f4*u0**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F35*F25r*c5*f2*u5**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F35*F24r*c4*f2*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F35*F23r*c3*f2*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F35*F22r*c2*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F35*F21r*c1*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F35*F20r*c0*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F34*F25r*c5*f2*u4*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svm*dsvp*F34*F24r*c4*f2*u4**2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F34*F23r*c3*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F34*F22r*c2*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F34*F21r*c1*f2*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F34*F20r*c0*f2*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F33*F25r*c5*f2*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F33*F24r*c4*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F33*F23r*c3*f2*u3**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F33*F22r*c2*f2*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F33*F21r*c1*f2*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F33*F20r*c0*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F32*F25r*c5*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F32*F24r*c4*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F32*F23r*c3*f2*u2*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svm*dsvp*F32*F22r*c2*f2*u2**2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F32*F21r*c1*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F32*F20r*c0*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F31*F25r*c5*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F31*F24r*c4*f2*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F31*F23r*c3*f2*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F31*F22r*c2*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F31*F21r*c1*f2*u1**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F31*F20r*c0*f2*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F30*F25r*c5*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F30*F24r*c4*f2*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F30*F23r*c3*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F30*F22r*c2*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F30*F21r*c1*f2*u0*u1*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svm*dsvp*F30*F20r*c0*f2*u0**2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F25*F35r*c5*f3*u5**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F25*F34r*c4*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F25*F33r*c3*f3*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F25*F32r*c2*f3*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F25*F31r*c1*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F25*F30r*c0*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F24*F35r*c5*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F24*F34r*c4*f3*u4**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F24*F33r*c3*f3*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F24*F32r*c2*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F24*F31r*c1*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F24*F30r*c0*f3*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F23*F35r*c5*f3*u3*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svm*dsvp*F23*F34r*c4*f3*u3*u4*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F23*F33r*c3*f3*u3**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F23*F32r*c2*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F23*F31r*c1*f3*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F23*F30r*c0*f3*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F22*F35r*c5*f3*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F22*F34r*c4*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F22*F33r*c3*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F22*F32r*c2*f3*u2**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F22*F31r*c1*f3*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F22*F30r*c0*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F21*F35r*c5*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F21*F34r*c4*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F21*F33r*c3*f3*u1*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svm*dsvp*F21*F32r*c2*f3*u1*u2*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F21*F31r*c1*f3*u1**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F21*F30r*c0*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F20*F35r*c5*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F20*F34r*c4*f3*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F20*F33r*c3*f3*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F20*F32r*c2*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F20*F31r*c1*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svm*dsvp*F20*F30r*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F45*F45r*c5*f4*u5**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F45*F44r*c4*f4*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F45*F43r*c3*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F45*F42r*c2*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F45*F41r*c1*f4*u1*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svp*dsvm*F45*F40r*c0*f4*u0*u5*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F44*F45r*c5*f4*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F44*F44r*c4*f4*u4**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F44*F43r*c3*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F44*F42r*c2*f4*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F44*F41r*c1*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F44*F40r*c0*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F43*F45r*c5*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F43*F44r*c4*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F43*F43r*c3*f4*u3**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F43*F42r*c2*f4*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F43*F41r*c1*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F43*F40r*c0*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F42*F45r*c5*f4*u2*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svp*dsvm*F42*F44r*c4*f4*u2*u4*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F42*F43r*c3*f4*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F42*F42r*c2*f4*u2**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F42*F41r*c1*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F42*F40r*c0*f4*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F41*F45r*c5*f4*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F41*F44r*c4*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F41*F43r*c3*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F41*F42r*c2*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F41*F41r*c1*f4*u1**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F41*F40r*c0*f4*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F40*F45r*c5*f4*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F40*F44r*c4*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F40*F43r*c3*f4*u0*u3*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*svp*dsvm*F40*F42r*c2*f4*u0*u2*
     & w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F40*F41r*c1*f4*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*svp*dsvm*F40*F40r*c0*f4*u0**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F35*F25r*c5*f2*u5**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F35*F24r*c4*f2*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F35*F23r*c3*f2*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F35*F22r*c2*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F35*F21r*c1*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F35*F20r*c0*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F34*F25r*c5*f2*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F34*F24r*c4*f2*u4**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F34*F23r*c3*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F34*F22r*c2*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F34*F21r*c1*f2*u1*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*dsvm*F34*F20r*c0*f2*u0*u4*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F33*F25r*c5*f2*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F33*F24r*c4*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F33*F23r*c3*f2*u3**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F33*F22r*c2*f2*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F33*F21r*c1*f2*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F33*F20r*c0*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F32*F25r*c5*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F32*F24r*c4*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F32*F23r*c3*f2*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F32*F22r*c2*f2*u2**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F32*F21r*c1*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F32*F20r*c0*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F31*F25r*c5*f2*u1*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*dsvm*F31*F24r*c4*f2*u1*u4*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F31*F23r*c3*f2*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F31*F22r*c2*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F31*F21r*c1*f2*u1**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F31*F20r*c0*f2*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F30*F25r*c5*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F30*F24r*c4*f2*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F30*F23r*c3*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F30*F22r*c2*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F30*F21r*c1*f2*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F30*F20r*c0*f2*u0**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F25*F35r*c5*f3*u5**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F25*F34r*c4*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F25*F33r*c3*f3*u3*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*dsvm*F25*F32r*c2*f3*u2*u5*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F25*F31r*c1*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F25*F30r*c0*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F24*F35r*c5*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F24*F34r*c4*f3*u4**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F24*F33r*c3*f3*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F24*F32r*c2*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F24*F31r*c1*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F24*F30r*c0*f3*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F23*F35r*c5*f3*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F23*F34r*c4*f3*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F23*F33r*c3*f3*u3**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F23*F32r*c2*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F23*F31r*c1*f3*u1*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*dsvm*F23*F30r*c0*f3*u0*u3*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F22*F35r*c5*f3*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F22*F34r*c4*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F22*F33r*c3*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F22*F32r*c2*f3*u2**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F22*F31r*c1*f3*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F22*F30r*c0*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F21*F35r*c5*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F21*F34r*c4*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F21*F33r*c3*f3*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F21*F32r*c2*f3*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F21*F31r*c1*f3*u1**2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F21*F30r*c0*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F20*F35r*c5*f3*u0*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*PC2*svp*dsvm*F20*F34r*c4*f3*u0*u4*
     & w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F20*F33r*c3*f3*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F20*F32r*c2*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F20*F31r*c1*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*PC2*svp*dsvm*F20*F30r*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F35*F35r*c5*f3*u5**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F35*F34r*c4*f3*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F35*F33r*c3*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F35*F32r*c2*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F35*F31r*c1*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F35*F30r*c0*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F34*F35r*c5*f3*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F34*F34r*c4*f3*u4**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F34*F33r*c3*f3*u3*u4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F34*F32r*c2*f3*u2*
     & u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F34*F31r*c1*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F34*F30r*c0*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F33*F35r*c5*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F33*F34r*c4*f3*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F33*F33r*c3*f3*u3**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F33*F32r*c2*f3*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F33*F31r*c1*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F33*F30r*c0*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F32*F35r*c5*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F32*F34r*c4*f3*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F32*F33r*c3*f3*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F32*F32r*c2*f3*u2**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F32*F31r*c1*f3*u1*u2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F32*F30r*c0*f3*u0*
     & u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F31*F35r*c5*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F31*F34r*c4*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F31*F33r*c3*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F31*F32r*c2*f3*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F31*F31r*c1*f3*u1**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F31*F30r*c0*f3*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F30*F35r*c5*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F30*F34r*c4*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F30*F33r*c3*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F30*F32r*c2*f3*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F30*F31r*c1*f3*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svm*dsvp*F30*F30r*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F35*F35r*c5*f3*u5**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F35*F34r*c4*f3*u4*
     & u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F35*F33r*c3*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F35*F32r*c2*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F35*F31r*c1*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F35*F30r*c0*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F34*F35r*c5*f3*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F34*F34r*c4*f3*u4**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F34*F33r*c3*f3*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F34*F32r*c2*f3*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F34*F31r*c1*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F34*F30r*c0*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F33*F35r*c5*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F33*F34r*c4*f3*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F33*F33r*c3*f3*u3**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F33*F32r*c2*f3*u2*
     & u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F33*F31r*c1*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F33*F30r*c0*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F32*F35r*c5*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F32*F34r*c4*f3*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F32*F33r*c3*f3*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F32*F32r*c2*f3*u2**2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F32*F31r*c1*f3*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F32*F30r*c0*f3*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F31*F35r*c5*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F31*F34r*c4*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F31*F33r*c3*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F31*F32r*c2*f3*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F31*F31r*c1*f3*u1**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F31*F30r*c0*f3*u0*
     & u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F30*F35r*c5*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F30*F34r*c4*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F30*F33r*c3*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F30*F32r*c2*f3*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F30*F31r*c1*f3*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*PC2*qq*svp*dsvm*F30*F30r*c0*f3*u0**2*w1*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F15i*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F12i*c2*f1*u2*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F45*F12i*c2*f1*u2*u5*
     & w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F45*F10i*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F14i*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F13i*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F11i*c1*f1*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F44*F11i*c1*f1*u1*u4*
     & w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F44*F10i*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F15i*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F13i*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F43*F10i*c0*f1*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F43*F10i*c0*f1*u0*u3*
     & w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F13i*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F12i*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F42*F10i*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F15i*c5*f1*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F41*F15i*c5*f1*u1*u5*
     & w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F13i*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F11i*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F41*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F14i*c4*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F40*F14i*c4*f1*u0*u4*
     & w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F10i*c0*f1*u0**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F40*F10i*c0*f1*u0**2*w1
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F25*F25i*c5*f2*u5**2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F25*F24i*c4*f2*u4*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F25*F23i*c3*f2*u3*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F25*F22i*c2*f2*u2*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F25*F21i*c1*f2*u1*u5
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svm*dsvp*F25*F20i*c0*f2*u0*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F24*F25i*c5*f2*u4*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F24*F24i*c4*f2*u4**2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F24*F23i*c3*f2*u3*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F24*F22i*c2*f2*u2*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F24*F21i*c1*f2*u1*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F24*F20i*c0*f2*u0*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F23*F25i*c5*f2*u3*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F23*F24i*c4*f2*u3*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F23*F23i*c3*f2*u3**2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F23*F22i*c2*f2*u2*u3
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F23*F21i*c1*f2*u1*u3
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F23*F20i*c0*f2*u0*u3
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F22*F25i*c5*f2*u2*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F22*F24i*c4*f2*u2*u4
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svm*dsvp*F22*F23i*c3*f2*u2*u3
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F22*F22i*c2*f2*u2**2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F22*F21i*c1*f2*u1*u2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F22*F20i*c0*f2*u0*u2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F21*F25i*c5*f2*u1*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F21*F24i*c4*f2*u1*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F21*F23i*c3*f2*u1*u3
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F21*F22i*c2*f2*u1*u2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F21*F21i*c1*f2*u1**2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F21*F20i*c0*f2*u0*u1
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F20*F25i*c5*f2*u0*u5
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F20*F24i*c4*f2*u0*u4
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F20*F23i*c3*f2*u0*u3
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F20*F22i*c2*f2*u0*u2
     &  + 8.D0*PC_q**2*i_*svm*dsvp*F20*F21i*c1*f2*u0*u1
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svm*dsvp*F20*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F15*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F45i*c5*f4*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F14*F44i*c4*f4*u4**2*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F14*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F44i*c4*f4*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F13*F43i*c3*f4*u3**2*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F13*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F43i*c3*f4*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F12*F42i*c2*f4*u2**2*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F12*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F42i*c2*f4*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F11*F41i*c1*f4*u1**2*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F11*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F41i*c1*f4*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dsvp*F10*F40i*c0*f4*u0**2*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dsvp*F10*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**2*i_*svm*dssp*F45*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F45*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F45*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F45*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F45*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F45*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F44*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F44*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F44*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F44*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F44*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F44*F20i*c0*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svm*dssp*F43*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F43*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F43*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F43*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F43*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F43*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F42*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F42*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F42*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F42*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F42*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F42*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F41*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F41*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F41*F23i*c3*f2*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svm*dssp*F41*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F41*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F41*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*svm*dssp*F40*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F40*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F40*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F40*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F40*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*svm*dssp*F40*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F35*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F35*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F35*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F35*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F35*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F35*F10i*c0*f1*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dssp*F34*F15i*c5*f1*u4*u5*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F34*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F34*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F34*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F34*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F34*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F33*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F33*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F33*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F33*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F33*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F33*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F32*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F32*F14i*c4*f1*u2*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dssp*F32*F13i*c3*f1*u2*u3*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F32*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F32*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F32*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F31*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F31*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F31*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F31*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F31*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F31*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F30*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F30*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F30*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F30*F12i*c2*f1*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dssp*F30*F11i*c1*f1*u0*u1*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F30*F10i*c0*f1*u0**2*w2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F25*F45i*c5*f4*u5**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F25*F44i*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F25*F43i*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F25*F42i*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F25*F41i*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F25*F40i*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F24*F45i*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F24*F44i*c4*f4*u4**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F24*F43i*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F24*F42i*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F24*F41i*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F24*F40i*c0*f4*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svm*dssp*F23*F45i*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F23*F44i*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F23*F43i*c3*f4*u3**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F23*F42i*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F23*F41i*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F23*F40i*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F22*F45i*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F22*F44i*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F22*F43i*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F22*F42i*c2*f4*u2**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F22*F41i*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F22*F40i*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F21*F45i*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F21*F44i*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F21*F43i*c3*f4*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svm*dssp*F21*F42i*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F21*F41i*c1*f4*u1**2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F21*F40i*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*svm*dssp*F20*F45i*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*svm*dssp*F20*F44i*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*svm*dssp*F20*F43i*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*svm*dssp*F20*F42i*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*svm*dssp*F20*F41i*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*svm*dssp*F20*F40i*c0*f4*u0**2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F15*F35i*c5*f3*u5**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F15*F34i*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F15*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F15*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F15*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F15*F30i*c0*f3*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dssp*F14*F35i*c5*f3*u4*u5*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F14*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F14*F33i*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F14*F32i*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F14*F31i*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F14*F30i*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F13*F35i*c5*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F13*F34i*c4*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F13*F33i*c3*f3*u3**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F13*F32i*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F13*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F13*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F12*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F12*F34i*c4*f3*u2*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dssp*F12*F33i*c3*f3*u2*u3*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F12*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F12*F31i*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F12*F30i*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F11*F35i*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F11*F34i*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F11*F33i*c3*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F11*F32i*c2*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F11*F31i*c1*f3*u1**2*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F11*F30i*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F10*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F10*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F10*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F10*F32i*c2*f3*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svm*dssp*F10*F31i*c1*f3*u0*u1*
     & w2
     &  - 4.D0*PC_q**2*i_*svm*dssp*F10*F30i*c0*f3*u0**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F15i*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F45*F10i*c0*f1*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F44*F15i*c5*f1*u4*u5*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F14i*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F13i*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F11i*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F44*F10i*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F15i*c5*f1*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F43*F14i*c4*f1*u3*u4*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F13i*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F43*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F14i*c4*f1*u2*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F42*F13i*c3*f1*u2*u3*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F12i*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F42*F10i*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F15i*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F13i*c3*f1*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F41*F12i*c2*f1*u1*u2*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F11i*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F41*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F12i*c2*f1*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F40*F11i*c1*f1*u0*u1*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F10i*c0*f1*u0**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F40*F10i*c0*f1*u0**2*w1
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F25*F25i*c5*f2*u5**2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F25*F24i*c4*f2*u4*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F25*F23i*c3*f2*u3*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F25*F22i*c2*f2*u2*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F25*F21i*c1*f2*u1*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F25*F20i*c0*f2*u0*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F24*F25i*c5*f2*u4*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F24*F24i*c4*f2*u4**2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F24*F23i*c3*f2*u3*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F24*F22i*c2*f2*u2*u4
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svp*dsvm*F24*F21i*c1*f2*u1*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F24*F20i*c0*f2*u0*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F23*F25i*c5*f2*u3*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F23*F24i*c4*f2*u3*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F23*F23i*c3*f2*u3**2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F23*F22i*c2*f2*u2*u3
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F23*F21i*c1*f2*u1*u3
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F23*F20i*c0*f2*u0*u3
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F22*F25i*c5*f2*u2*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F22*F24i*c4*f2*u2*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F22*F23i*c3*f2*u2*u3
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F22*F22i*c2*f2*u2**2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F22*F21i*c1*f2*u1*u2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F22*F20i*c0*f2*u0*u2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F21*F25i*c5*f2*u1*u5
     &
      traza1 = traza1 + 8.D0*PC_q**2*i_*svp*dsvm*F21*F24i*c4*f2*u1*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F21*F23i*c3*f2*u1*u3
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F21*F22i*c2*f2*u1*u2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F21*F21i*c1*f2*u1**2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F21*F20i*c0*f2*u0*u1
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F20*F25i*c5*f2*u0*u5
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F20*F24i*c4*f2*u0*u4
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F20*F23i*c3*f2*u0*u3
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F20*F22i*c2*f2*u0*u2
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F20*F21i*c1*f2*u0*u1
     &  + 8.D0*PC_q**2*i_*svp*dsvm*F20*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F44i*c4*f4*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F15*F43i*c3*f4*u3*u5*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F15*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F44i*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F43i*c3*f4*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F14*F42i*c2*f4*u2*u4*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F14*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F43i*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F42i*c2*f4*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F13*F41i*c1*f4*u1*u3*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F13*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F42i*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F41i*c1*f4*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F12*F40i*c0*f4*u0*u2*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F12*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F41i*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F11*F40i*c0*f4*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dsvm*F10*F45i*c5*f4*u0*u5*
     & w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F40i*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**2*i_*svp*dsvm*F10*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**2*i_*svp*dssm*F45*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F45*F24i*c4*f2*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svp*dssm*F45*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F45*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F45*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F45*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F44*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F44*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F44*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F44*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F44*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F44*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F43*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F43*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F43*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F43*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F43*F21i*c1*f2*u1*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svp*dssm*F43*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F42*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F42*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F42*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F42*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F42*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F42*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F41*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F41*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F41*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F41*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F41*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F41*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*svp*dssm*F40*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F40*F24i*c4*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svp*dssm*F40*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F40*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F40*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*svp*dssm*F40*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q**2*i_*svp*dssm*F35*F15i*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F35*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F35*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F35*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F35*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F35*F10i*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F34*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F34*F14i*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F34*F13i*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F34*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F34*F11i*c1*f1*u1*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dssm*F34*F10i*c0*f1*u0*u4*
     & w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F33*F15i*c5*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F33*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F33*F13i*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F33*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F33*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F33*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F32*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F32*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F32*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F32*F12i*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F32*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F32*F10i*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F31*F15i*c5*f1*u1*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dssm*F31*F14i*c4*f1*u1*u4*
     & w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F31*F13i*c3*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F31*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F31*F11i*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F31*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F30*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F30*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F30*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F30*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F30*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F30*F10i*c0*f1*u0**2*w1
     &  + 4.D0*PC_q**2*i_*svp*dssm*F25*F45i*c5*f4*u5**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F25*F44i*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F25*F43i*c3*f4*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svp*dssm*F25*F42i*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F25*F41i*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F25*F40i*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F24*F45i*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F24*F44i*c4*f4*u4**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F24*F43i*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F24*F42i*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F24*F41i*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F24*F40i*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F23*F45i*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F23*F44i*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F23*F43i*c3*f4*u3**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F23*F42i*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F23*F41i*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F23*F40i*c0*f4*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svp*dssm*F22*F45i*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F22*F44i*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F22*F43i*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F22*F42i*c2*f4*u2**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F22*F41i*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F22*F40i*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F21*F45i*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F21*F44i*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F21*F43i*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*svp*dssm*F21*F42i*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F21*F41i*c1*f4*u1**2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F21*F40i*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*svp*dssm*F20*F45i*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*svp*dssm*F20*F44i*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*svp*dssm*F20*F43i*c3*f4*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*svp*dssm*F20*F42i*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*svp*dssm*F20*F41i*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*svp*dssm*F20*F40i*c0*f4*u0**2
     &  - 4.D0*PC_q**2*i_*svp*dssm*F15*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F15*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F15*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F15*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F15*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F15*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F14*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F14*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F14*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F14*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F14*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F14*F30i*c0*f3*u0*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dssm*F13*F35i*c5*f3*u3*u5*
     & w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F13*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F13*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F13*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F13*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F13*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F12*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F12*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F12*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F12*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F12*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F12*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F11*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F11*F34i*c4*f3*u1*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*svp*dssm*F11*F33i*c3*f3*u1*u3*
     & w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F11*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F11*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F11*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F10*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F10*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F10*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F10*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F10*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*svp*dssm*F10*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F45*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F45*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F45*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F45*F22i*c2*f2*u2*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dsvp*F45*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F45*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F44*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F44*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F44*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F44*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F44*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F44*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F43*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F43*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F43*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F43*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F43*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F43*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F42*F25i*c5*f2*u2*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dsvp*F42*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F42*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F42*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F42*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F42*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F41*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F41*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F41*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F41*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F41*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F41*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F40*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F40*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F40*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F40*F22i*c2*f2*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dsvp*F40*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F40*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F35*F15i*c5*f1*u5**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F35*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F35*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F35*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F35*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F35*F10i*c0*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F34*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F34*F14i*c4*f1*u4**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F34*F13i*c3*f1*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F34*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F34*F11i*c1*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F34*F10i*c0*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F33*F15i*c5*f1*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*dsvp*F33*F14i*c4*f1*u3*u4*
     & w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F33*F13i*c3*f1*u3**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F33*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F33*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F33*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F32*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F32*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F32*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F32*F12i*c2*f1*u2**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F32*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F32*F10i*c0*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F31*F15i*c5*f1*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F31*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F31*F13i*c3*f1*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*dsvp*F31*F12i*c2*f1*u1*u2*
     & w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F31*F11i*c1*f1*u1**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F31*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F30*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F30*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F30*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F30*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F30*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F30*F10i*c0*f1*u0**2*w1
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F25*F45i*c5*f4*u5**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F25*F44i*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F25*F43i*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F25*F42i*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F25*F41i*c1*f4*u1*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dsvp*F25*F40i*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F24*F45i*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F24*F44i*c4*f4*u4**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F24*F43i*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F24*F42i*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F24*F41i*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F24*F40i*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F23*F45i*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F23*F44i*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F23*F43i*c3*f4*u3**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F23*F42i*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F23*F41i*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F23*F40i*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F22*F45i*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F22*F44i*c4*f4*u2*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dsvp*F22*F43i*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F22*F42i*c2*f4*u2**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F22*F41i*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F22*F40i*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F21*F45i*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F21*F44i*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F21*F43i*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F21*F42i*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F21*F41i*c1*f4*u1**2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F21*F40i*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F20*F45i*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F20*F44i*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F20*F43i*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F20*F42i*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dsvp*F20*F41i*c1*f4*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dsvp*F20*F40i*c0*f4*u0**2
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F15*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F15*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F15*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F15*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F15*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F15*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F14*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F14*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F14*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F14*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F14*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F14*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F13*F35i*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F13*F34i*c4*f3*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*dsvp*F13*F33i*c3*f3*u3**2*
     & w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F13*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F13*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F13*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F12*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F12*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F12*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F12*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F12*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F12*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F11*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F11*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F11*F33i*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F11*F32i*c2*f3*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssm*dsvp*F11*F31i*c1*f3*u1**2*
     & w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F11*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F10*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F10*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F10*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F10*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F10*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**2*i_*ssm*dsvp*F10*F30i*c0*f3*u0**2*w1
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F45*F45i*c5*f4*u5**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F45*F44i*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F45*F43i*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F45*F42i*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F45*F41i*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F45*F40i*c0*f4*u0*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F44*F45i*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F44*F44i*c4*f4*u4**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F44*F43i*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F44*F42i*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F44*F41i*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F44*F40i*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F43*F45i*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F43*F44i*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F43*F43i*c3*f4*u3**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F43*F42i*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F43*F41i*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F43*F40i*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F42*F45i*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F42*F44i*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F42*F43i*c3*f4*u2*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F42*F42i*c2*f4*u2**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F42*F41i*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F42*F40i*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F41*F45i*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F41*F44i*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F41*F43i*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F41*F42i*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F41*F41i*c1*f4*u1**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F41*F40i*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F40*F45i*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F40*F44i*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F40*F43i*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F40*F42i*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F40*F41i*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F40*F40i*c0*f4*u0**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F35*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F35*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F35*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F35*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F35*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F35*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F34*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F34*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F34*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F34*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F34*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F34*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F33*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F33*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F33*F23i*c3*f2*u3**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F33*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F33*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F33*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F32*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F32*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F32*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F32*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F32*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F32*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F31*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F31*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F31*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F31*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F31*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F31*F20i*c0*f2*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F30*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F30*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F30*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F30*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F30*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F30*F20i*c0*f2*u0**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F25*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F25*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F25*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F25*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F25*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F25*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F24*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F24*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F24*F33i*c3*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F24*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F24*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F24*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F23*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F23*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F23*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F23*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F23*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F23*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F22*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F22*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F22*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F22*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F22*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F22*F30i*c0*f3*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssm*dssp*F21*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F21*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F21*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F21*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F21*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F21*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F20*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F20*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F20*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F20*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F20*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*ssm*dssp*F20*F30i*c0*f3*u0**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F45*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F45*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F45*F23i*c3*f2*u3*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dsvm*F45*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F45*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F45*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F44*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F44*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F44*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F44*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F44*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F44*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F43*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F43*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F43*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F43*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F43*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F43*F20i*c0*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dsvm*F42*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F42*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F42*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F42*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F42*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F42*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F41*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F41*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F41*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F41*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F41*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F41*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F40*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F40*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F40*F23i*c3*f2*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dsvm*F40*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F40*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F40*F20i*c0*f2*u0**2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F35*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F35*F14i*c4*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F35*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F35*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F35*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F35*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F34*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F34*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F34*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F34*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F34*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F34*F10i*c0*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*dsvm*F33*F15i*c5*f1*u3*u5*
     & w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F33*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F33*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F33*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F33*F11i*c1*f1*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F33*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F32*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F32*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F32*F13i*c3*f1*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F32*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F32*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F32*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F31*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F31*F14i*c4*f1*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*dsvm*F31*F13i*c3*f1*u1*u3*
     & w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F31*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F31*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F31*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F30*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F30*F14i*c4*f1*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F30*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F30*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F30*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F30*F10i*c0*f1*u0**2*w2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F25*F45i*c5*f4*u5**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F25*F44i*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F25*F43i*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F25*F42i*c2*f4*u2*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dsvm*F25*F41i*c1*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F25*F40i*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F24*F45i*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F24*F44i*c4*f4*u4**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F24*F43i*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F24*F42i*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F24*F41i*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F24*F40i*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F23*F45i*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F23*F44i*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F23*F43i*c3*f4*u3**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F23*F42i*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F23*F41i*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F23*F40i*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F22*F45i*c5*f4*u2*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dsvm*F22*F44i*c4*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F22*F43i*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F22*F42i*c2*f4*u2**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F22*F41i*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F22*F40i*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F21*F45i*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F21*F44i*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F21*F43i*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F21*F42i*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F21*F41i*c1*f4*u1**2
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F21*F40i*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F20*F45i*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F20*F44i*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F20*F43i*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F20*F42i*c2*f4*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dsvm*F20*F41i*c1*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dsvm*F20*F40i*c0*f4*u0**2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F15*F35i*c5*f3*u5**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F15*F34i*c4*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F15*F33i*c3*f3*u3*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F15*F32i*c2*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F15*F31i*c1*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F15*F30i*c0*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F14*F35i*c5*f3*u4*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F14*F34i*c4*f3*u4**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F14*F33i*c3*f3*u3*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F14*F32i*c2*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F14*F31i*c1*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F14*F30i*c0*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F13*F35i*c5*f3*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*dsvm*F13*F34i*c4*f3*u3*u4*
     & w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F13*F33i*c3*f3*u3**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F13*F32i*c2*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F13*F31i*c1*f3*u1*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F13*F30i*c0*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F12*F35i*c5*f3*u2*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F12*F34i*c4*f3*u2*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F12*F33i*c3*f3*u2*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F12*F32i*c2*f3*u2**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F12*F31i*c1*f3*u1*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F12*F30i*c0*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F11*F35i*c5*f3*u1*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F11*F34i*c4*f3*u1*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F11*F33i*c3*f3*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*ssp*dsvm*F11*F32i*c2*f3*u1*u2*
     & w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F11*F31i*c1*f3*u1**2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F11*F30i*c0*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F10*F35i*c5*f3*u0*u5*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F10*F34i*c4*f3*u0*u4*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F10*F33i*c3*f3*u0*u3*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F10*F32i*c2*f3*u0*u2*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F10*F31i*c1*f3*u0*u1*w2
     &  - 4.D0*PC_q**2*i_*ssp*dsvm*F10*F30i*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F45*F45i*c5*f4*u5**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F45*F44i*c4*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F45*F43i*c3*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F45*F42i*c2*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F45*F41i*c1*f4*u1*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F45*F40i*c0*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F44*F45i*c5*f4*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F44*F44i*c4*f4*u4**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F44*F43i*c3*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F44*F42i*c2*f4*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F44*F41i*c1*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F44*F40i*c0*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F43*F45i*c5*f4*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F43*F44i*c4*f4*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F43*F43i*c3*f4*u3**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F43*F42i*c2*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F43*F41i*c1*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F43*F40i*c0*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F42*F45i*c5*f4*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F42*F44i*c4*f4*u2*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F42*F43i*c3*f4*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F42*F42i*c2*f4*u2**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F42*F41i*c1*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F42*F40i*c0*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F41*F45i*c5*f4*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F41*F44i*c4*f4*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F41*F43i*c3*f4*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F41*F42i*c2*f4*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F41*F41i*c1*f4*u1**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F41*F40i*c0*f4*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F40*F45i*c5*f4*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F40*F44i*c4*f4*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F40*F43i*c3*f4*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F40*F42i*c2*f4*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F40*F41i*c1*f4*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F40*F40i*c0*f4*u0**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F35*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F35*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F35*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F35*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F35*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F35*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F34*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F34*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F34*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F34*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F34*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F34*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F33*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F33*F24i*c4*f2*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F33*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F33*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F33*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F33*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F32*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F32*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F32*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F32*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F32*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F32*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F31*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F31*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F31*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F31*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F31*F21i*c1*f2*u1**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F31*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F30*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F30*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F30*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F30*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F30*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F30*F20i*c0*f2*u0**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F25*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F25*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F25*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F25*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F25*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F25*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F24*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F24*F34i*c4*f3*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F24*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F24*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F24*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F24*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F23*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F23*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F23*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F23*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F23*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F23*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F22*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F22*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F22*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F22*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F22*F31i*c1*f3*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*ssp*dssm*F22*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F21*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F21*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F21*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F21*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F21*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F21*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F20*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F20*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F20*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F20*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F20*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*ssp*dssm*F20*F30i*c0*f3*u0**2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F45*F45i*c5*f4*u5**2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F45*F44i*c4*f4*u4*u5
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svm*dsvp*F45*F43i*c3*f4*u3*
     & u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F45*F42i*c2*f4*u2*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F45*F41i*c1*f4*u1*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F45*F40i*c0*f4*u0*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F44*F45i*c5*f4*u4*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F44*F44i*c4*f4*u4**2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F44*F43i*c3*f4*u3*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F44*F42i*c2*f4*u2*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F44*F41i*c1*f4*u1*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F44*F40i*c0*f4*u0*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F43*F45i*c5*f4*u3*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F43*F44i*c4*f4*u3*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F43*F43i*c3*f4*u3**2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F43*F42i*c2*f4*u2*u3
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svm*dsvp*F43*F41i*c1*f4*u1*
     & u3
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F43*F40i*c0*f4*u0*u3
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F42*F45i*c5*f4*u2*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F42*F44i*c4*f4*u2*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F42*F43i*c3*f4*u2*u3
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F42*F42i*c2*f4*u2**2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F42*F41i*c1*f4*u1*u2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F42*F40i*c0*f4*u0*u2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F41*F45i*c5*f4*u1*u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F41*F44i*c4*f4*u1*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F41*F43i*c3*f4*u1*u3
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F41*F42i*c2*f4*u1*u2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F41*F41i*c1*f4*u1**2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F41*F40i*c0*f4*u0*u1
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svm*dsvp*F40*F45i*c5*f4*u0*
     & u5
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F40*F44i*c4*f4*u0*u4
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F40*F43i*c3*f4*u0*u3
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F40*F42i*c2*f4*u0*u2
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F40*F41i*c1*f4*u0*u1
     &  - 4.D0*PC_q**2*i_*qq*svm*dsvp*F40*F40i*c0*f4*u0**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F35*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F35*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F35*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F35*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F35*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F35*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F34*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F34*F24i*c4*f2*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svm*dsvp*F34*F23i*c3*f2*u3*
     & u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F34*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F34*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F34*F20i*c0*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F33*F25i*c5*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F33*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F33*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F33*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F33*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F33*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F32*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F32*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F32*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F32*F22i*c2*f2*u2**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svm*dsvp*F32*F21i*c1*f2*u1*
     & u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F32*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F31*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F31*F24i*c4*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F31*F23i*c3*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F31*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F31*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F31*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F30*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F30*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F30*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F30*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F30*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F30*F20i*c0*f2*u0**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svm*dsvp*F25*F35i*c5*f3*
     & u5**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F25*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F25*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F25*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F25*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F25*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F24*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F24*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F24*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F24*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F24*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F24*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F23*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F23*F34i*c4*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svm*dsvp*F23*F33i*c3*f3*
     & u3**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F23*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F23*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F23*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F22*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F22*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F22*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F22*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F22*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F22*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F21*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F21*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F21*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F21*F32i*c2*f3*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svm*dsvp*F21*F31i*c1*f3*
     & u1**2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F21*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F20*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F20*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F20*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F20*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F20*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svm*dsvp*F20*F30i*c0*f3*u0**2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F45*F45i*c5*f4*u5**2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F45*F44i*c4*f4*u4*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F45*F43i*c3*f4*u3*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F45*F42i*c2*f4*u2*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F45*F41i*c1*f4*u1*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F45*F40i*c0*f4*u0*u5
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svp*dsvm*F44*F45i*c5*f4*u4*
     & u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F44*F44i*c4*f4*u4**2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F44*F43i*c3*f4*u3*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F44*F42i*c2*f4*u2*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F44*F41i*c1*f4*u1*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F44*F40i*c0*f4*u0*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F43*F45i*c5*f4*u3*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F43*F44i*c4*f4*u3*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F43*F43i*c3*f4*u3**2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F43*F42i*c2*f4*u2*u3
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F43*F41i*c1*f4*u1*u3
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F43*F40i*c0*f4*u0*u3
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F42*F45i*c5*f4*u2*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F42*F44i*c4*f4*u2*u4
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svp*dsvm*F42*F43i*c3*f4*u2*
     & u3
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F42*F42i*c2*f4*u2**2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F42*F41i*c1*f4*u1*u2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F42*F40i*c0*f4*u0*u2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F41*F45i*c5*f4*u1*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F41*F44i*c4*f4*u1*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F41*F43i*c3*f4*u1*u3
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F41*F42i*c2*f4*u1*u2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F41*F41i*c1*f4*u1**2
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F41*F40i*c0*f4*u0*u1
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F40*F45i*c5*f4*u0*u5
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F40*F44i*c4*f4*u0*u4
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F40*F43i*c3*f4*u0*u3
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F40*F42i*c2*f4*u0*u2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*qq*svp*dsvm*F40*F41i*c1*f4*u0*
     & u1
     &  - 4.D0*PC_q**2*i_*qq*svp*dsvm*F40*F40i*c0*f4*u0**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F35*F25i*c5*f2*u5**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F35*F24i*c4*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F35*F23i*c3*f2*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F35*F22i*c2*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F35*F21i*c1*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F35*F20i*c0*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F34*F25i*c5*f2*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F34*F24i*c4*f2*u4**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F34*F23i*c3*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F34*F22i*c2*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F34*F21i*c1*f2*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F34*F20i*c0*f2*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*dsvm*F33*F25i*c5*f2*u3*
     & u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F33*F24i*c4*f2*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F33*F23i*c3*f2*u3**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F33*F22i*c2*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F33*F21i*c1*f2*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F33*F20i*c0*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F32*F25i*c5*f2*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F32*F24i*c4*f2*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F32*F23i*c3*f2*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F32*F22i*c2*f2*u2**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F32*F21i*c1*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F32*F20i*c0*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F31*F25i*c5*f2*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F31*F24i*c4*f2*u1*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*dsvm*F31*F23i*c3*f2*u1*
     & u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F31*F22i*c2*f2*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F31*F21i*c1*f2*u1**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F31*F20i*c0*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F30*F25i*c5*f2*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F30*F24i*c4*f2*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F30*F23i*c3*f2*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F30*F22i*c2*f2*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F30*F21i*c1*f2*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F30*F20i*c0*f2*u0**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F25*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F25*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F25*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F25*F32i*c2*f3*u2*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*dsvm*F25*F31i*c1*f3*u1*
     & u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F25*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F24*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F24*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F24*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F24*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F24*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F24*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F23*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F23*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F23*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F23*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F23*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F23*F30i*c0*f3*u0*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*dsvm*F22*F35i*c5*f3*u2*
     & u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F22*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F22*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F22*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F22*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F22*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F21*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F21*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F21*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F21*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F21*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F21*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F20*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F20*F34i*c4*f3*u0*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*svp*dsvm*F20*F33i*c3*f3*u0*
     & u3
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F20*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F20*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*svp*dsvm*F20*F30i*c0*f3*u0**2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F35*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F35*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F35*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F35*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F35*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F35*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F34*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F34*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F34*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F34*F32i*c2*f3*u2*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssm*dssp*F34*F31i*c1*f3*u1*
     & u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F34*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F33*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F33*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F33*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F33*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F33*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F33*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F32*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F32*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F32*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F32*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F32*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F32*F30i*c0*f3*u0*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssm*dssp*F31*F35i*c5*f3*u1*
     & u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F31*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F31*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F31*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F31*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F31*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F30*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F30*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F30*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F30*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F30*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*ssm*dssp*F30*F30i*c0*f3*u0**2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F35*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F35*F34i*c4*f3*u4*u5
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssp*dssm*F35*F33i*c3*f3*u3*
     & u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F35*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F35*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F35*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F34*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F34*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F34*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F34*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F34*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F34*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F33*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F33*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F33*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F33*F32i*c2*f3*u2*u3
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssp*dssm*F33*F31i*c1*f3*u1*
     & u3
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F33*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F32*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F32*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F32*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F32*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F32*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F32*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F31*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F31*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F31*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F31*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F31*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F31*F30i*c0*f3*u0*u1
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq*ssp*dssm*F30*F35i*c5*f3*u0*
     & u5
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F30*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F30*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F30*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F30*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq*ssp*dssm*F30*F30i*c0*f3*u0**2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F35*F35i*c5*f3*u5**2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F35*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F35*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F35*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F35*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F35*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F34*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F34*F34i*c4*f3*u4**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F34*F33i*c3*f3*
     & u3*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F34*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F34*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F34*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F33*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F33*F34i*c4*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F33*F33i*c3*f3*u3**2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F33*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F33*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F33*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F32*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F32*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F32*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F32*F32i*c2*f3*u2**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F32*F31i*c1*f3*
     & u1*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F32*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F31*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F31*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F31*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F31*F32i*c2*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F31*F31i*c1*f3*u1**2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F31*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F30*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F30*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F30*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F30*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F30*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq**2*svm*dsvp*F30*F30i*c0*f3*u0**2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F35*F35i*c5*f3*
     & u5**2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F35*F34i*c4*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F35*F33i*c3*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F35*F32i*c2*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F35*F31i*c1*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F35*F30i*c0*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F34*F35i*c5*f3*u4*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F34*F34i*c4*f3*u4**2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F34*F33i*c3*f3*u3*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F34*F32i*c2*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F34*F31i*c1*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F34*F30i*c0*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F33*F35i*c5*f3*u3*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F33*F34i*c4*f3*u3*u4
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F33*F33i*c3*f3*
     & u3**2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F33*F32i*c2*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F33*F31i*c1*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F33*F30i*c0*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F32*F35i*c5*f3*u2*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F32*F34i*c4*f3*u2*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F32*F33i*c3*f3*u2*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F32*F32i*c2*f3*u2**2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F32*F31i*c1*f3*u1*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F32*F30i*c0*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F31*F35i*c5*f3*u1*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F31*F34i*c4*f3*u1*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F31*F33i*c3*f3*u1*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F31*F32i*c2*f3*u1*u2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F31*F31i*c1*f3*
     & u1**2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F31*F30i*c0*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F30*F35i*c5*f3*u0*u5
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F30*F34i*c4*f3*u0*u4
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F30*F33i*c3*f3*u0*u3
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F30*F32i*c2*f3*u0*u2
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F30*F31i*c1*f3*u0*u1
     &  + 4.D0*PC_q**2*i_*qq**2*svp*dsvm*F30*F30i*c0*f3*u0**2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F45*F45i*c5*f4*u5**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F45*F44i*c4*f4*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F45*F43i*c3*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F45*F42i*c2*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F45*F41i*c1*f4*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F45*F40i*c0*f4*u0*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F44*F45i*c5*f4*u4*
     & u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F44*F44i*c4*f4*u4**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F44*F43i*c3*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F44*F42i*c2*f4*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F44*F41i*c1*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F44*F40i*c0*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F43*F45i*c5*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F43*F44i*c4*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F43*F43i*c3*f4*u3**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F43*F42i*c2*f4*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F43*F41i*c1*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F43*F40i*c0*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F42*F45i*c5*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F42*F44i*c4*f4*u2*u4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F42*F43i*c3*f4*u2*
     & u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F42*F42i*c2*f4*u2**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F42*F41i*c1*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F42*F40i*c0*f4*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F41*F45i*c5*f4*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F41*F44i*c4*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F41*F43i*c3*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F41*F42i*c2*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F41*F41i*c1*f4*u1**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F41*F40i*c0*f4*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F40*F45i*c5*f4*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F40*F44i*c4*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F40*F43i*c3*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F40*F42i*c2*f4*u0*u2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F40*F41i*c1*f4*u0*
     & u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svm*dsvp*F40*F40i*c0*f4*u0**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F35*F25i*c5*f2*u5**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F35*F24i*c4*f2*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F35*F23i*c3*f2*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F35*F22i*c2*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F35*F21i*c1*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F35*F20i*c0*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F34*F25i*c5*f2*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F34*F24i*c4*f2*u4**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F34*F23i*c3*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F34*F22i*c2*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F34*F21i*c1*f2*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F34*F20i*c0*f2*u0*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F33*F25i*c5*f2*u3*
     & u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F33*F24i*c4*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F33*F23i*c3*f2*u3**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F33*F22i*c2*f2*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F33*F21i*c1*f2*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F33*F20i*c0*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F32*F25i*c5*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F32*F24i*c4*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F32*F23i*c3*f2*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F32*F22i*c2*f2*u2**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F32*F21i*c1*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F32*F20i*c0*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F31*F25i*c5*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F31*F24i*c4*f2*u1*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F31*F23i*c3*f2*u1*
     & u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F31*F22i*c2*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F31*F21i*c1*f2*u1**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F31*F20i*c0*f2*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F30*F25i*c5*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F30*F24i*c4*f2*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F30*F23i*c3*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F30*F22i*c2*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F30*F21i*c1*f2*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F30*F20i*c0*f2*u0**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F25*F35i*c5*f3*u5**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F25*F34i*c4*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F25*F33i*c3*f3*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F25*F32i*c2*f3*u2*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F25*F31i*c1*f3*u1*
     & u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F25*F30i*c0*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F24*F35i*c5*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F24*F34i*c4*f3*u4**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F24*F33i*c3*f3*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F24*F32i*c2*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F24*F31i*c1*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F24*F30i*c0*f3*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F23*F35i*c5*f3*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F23*F34i*c4*f3*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F23*F33i*c3*f3*u3**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F23*F32i*c2*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F23*F31i*c1*f3*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F23*F30i*c0*f3*u0*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F22*F35i*c5*f3*u2*
     & u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F22*F34i*c4*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F22*F33i*c3*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F22*F32i*c2*f3*u2**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F22*F31i*c1*f3*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F22*F30i*c0*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F21*F35i*c5*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F21*F34i*c4*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F21*F33i*c3*f3*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F21*F32i*c2*f3*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F21*F31i*c1*f3*u1**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F21*F30i*c0*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F20*F35i*c5*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F20*F34i*c4*f3*u0*u4*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F20*F33i*c3*f3*u0*
     & u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F20*F32i*c2*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F20*F31i*c1*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svm*dsvp*F20*F30i*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F45*F45i*c5*f4*u5**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F45*F44i*c4*f4*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F45*F43i*c3*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F45*F42i*c2*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F45*F41i*c1*f4*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F45*F40i*c0*f4*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F44*F45i*c5*f4*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F44*F44i*c4*f4*u4**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F44*F43i*c3*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F44*F42i*c2*f4*u2*u4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F44*F41i*c1*f4*u1*
     & u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F44*F40i*c0*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F43*F45i*c5*f4*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F43*F44i*c4*f4*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F43*F43i*c3*f4*u3**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F43*F42i*c2*f4*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F43*F41i*c1*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F43*F40i*c0*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F42*F45i*c5*f4*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F42*F44i*c4*f4*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F42*F43i*c3*f4*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F42*F42i*c2*f4*u2**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F42*F41i*c1*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F42*F40i*c0*f4*u0*u2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F41*F45i*c5*f4*u1*
     & u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F41*F44i*c4*f4*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F41*F43i*c3*f4*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F41*F42i*c2*f4*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F41*F41i*c1*f4*u1**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F41*F40i*c0*f4*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F40*F45i*c5*f4*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F40*F44i*c4*f4*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F40*F43i*c3*f4*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F40*F42i*c2*f4*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F40*F41i*c1*f4*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*svp*dsvm*F40*F40i*c0*f4*u0**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F35*F25i*c5*f2*u5**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F35*F24i*c4*f2*u4*u5*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F35*F23i*c3*f2*u3*
     & u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F35*F22i*c2*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F35*F21i*c1*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F35*F20i*c0*f2*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F34*F25i*c5*f2*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F34*F24i*c4*f2*u4**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F34*F23i*c3*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F34*F22i*c2*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F34*F21i*c1*f2*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F34*F20i*c0*f2*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F33*F25i*c5*f2*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F33*F24i*c4*f2*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F33*F23i*c3*f2*u3**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F33*F22i*c2*f2*u2*u3*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F33*F21i*c1*f2*u1*
     & u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F33*F20i*c0*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F32*F25i*c5*f2*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F32*F24i*c4*f2*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F32*F23i*c3*f2*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F32*F22i*c2*f2*u2**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F32*F21i*c1*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F32*F20i*c0*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F31*F25i*c5*f2*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F31*F24i*c4*f2*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F31*F23i*c3*f2*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F31*F22i*c2*f2*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F31*F21i*c1*f2*u1**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F31*F20i*c0*f2*u0*u1*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F30*F25i*c5*f2*u0*
     & u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F30*F24i*c4*f2*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F30*F23i*c3*f2*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F30*F22i*c2*f2*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F30*F21i*c1*f2*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F30*F20i*c0*f2*u0**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F25*F35i*c5*f3*u5**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F25*F34i*c4*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F25*F33i*c3*f3*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F25*F32i*c2*f3*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F25*F31i*c1*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F25*F30i*c0*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F24*F35i*c5*f3*u4*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F24*F34i*c4*f3*u4**2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F24*F33i*c3*f3*u3*
     & u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F24*F32i*c2*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F24*F31i*c1*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F24*F30i*c0*f3*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F23*F35i*c5*f3*u3*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F23*F34i*c4*f3*u3*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F23*F33i*c3*f3*u3**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F23*F32i*c2*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F23*F31i*c1*f3*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F23*F30i*c0*f3*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F22*F35i*c5*f3*u2*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F22*F34i*c4*f3*u2*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F22*F33i*c3*f3*u2*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F22*F32i*c2*f3*u2**2*w1*w2
     &
      traza1 = traza1 - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F22*F31i*c1*f3*u1*
     & u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F22*F30i*c0*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F21*F35i*c5*f3*u1*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F21*F34i*c4*f3*u1*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F21*F33i*c3*f3*u1*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F21*F32i*c2*f3*u1*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F21*F31i*c1*f3*u1**2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F21*F30i*c0*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F20*F35i*c5*f3*u0*u5*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F20*F34i*c4*f3*u0*u4*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F20*F33i*c3*f3*u0*u3*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F20*F32i*c2*f3*u0*u2*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F20*F31i*c1*f3*u0*u1*w1*w2
     &  - 4.D0*PC_q**2*i_*PC2*svp*dsvm*F20*F30i*c0*f3*u0**2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F35*F35i*c5*f3*
     & u5**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F35*F34i*c4*f3*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F35*F33i*c3*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F35*F32i*c2*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F35*F31i*c1*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F35*F30i*c0*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F34*F35i*c5*f3*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F34*F34i*c4*f3*u4**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F34*F33i*c3*f3*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F34*F32i*c2*f3*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F34*F31i*c1*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F34*F30i*c0*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F33*F35i*c5*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F33*F34i*c4*f3*u3*u4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F33*F33i*c3*f3*
     & u3**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F33*F32i*c2*f3*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F33*F31i*c1*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F33*F30i*c0*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F32*F35i*c5*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F32*F34i*c4*f3*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F32*F33i*c3*f3*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F32*F32i*c2*f3*u2**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F32*F31i*c1*f3*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F32*F30i*c0*f3*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F31*F35i*c5*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F31*F34i*c4*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F31*F33i*c3*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F31*F32i*c2*f3*u1*u2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F31*F31i*c1*f3*
     & u1**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F31*F30i*c0*f3*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F30*F35i*c5*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F30*F34i*c4*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F30*F33i*c3*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F30*F32i*c2*f3*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F30*F31i*c1*f3*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svm*dsvp*F30*F30i*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F35*F35i*c5*f3*u5**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F35*F34i*c4*f3*u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F35*F33i*c3*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F35*F32i*c2*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F35*F31i*c1*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F35*F30i*c0*f3*u0*u5*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F34*F35i*c5*f3*
     & u4*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F34*F34i*c4*f3*u4**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F34*F33i*c3*f3*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F34*F32i*c2*f3*u2*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F34*F31i*c1*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F34*F30i*c0*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F33*F35i*c5*f3*u3*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F33*F34i*c4*f3*u3*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F33*F33i*c3*f3*u3**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F33*F32i*c2*f3*u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F33*F31i*c1*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F33*F30i*c0*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F32*F35i*c5*f3*u2*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F32*F34i*c4*f3*u2*u4*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F32*F33i*c3*f3*
     & u2*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F32*F32i*c2*f3*u2**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F32*F31i*c1*f3*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F32*F30i*c0*f3*u0*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F31*F35i*c5*f3*u1*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F31*F34i*c4*f3*u1*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F31*F33i*c3*f3*u1*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F31*F32i*c2*f3*u1*u2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F31*F31i*c1*f3*u1**2*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F31*F30i*c0*f3*u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F30*F35i*c5*f3*u0*u5*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F30*F34i*c4*f3*u0*u4*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F30*F33i*c3*f3*u0*u3*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F30*F32i*c2*f3*u0*u2*w1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F30*F31i*c1*f3*
     & u0*u1*w1*w2
     &  + 4.D0*PC_q**2*i_*PC2*qq*svp*dsvm*F30*F30i*c0*f3*u0**2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F35*F35i*c5*f3*u5**2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F35*F34i*c4*f3*u4*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F35*F33i*c3*f3*u3*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F35*F32i*c2*f3*u2*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F35*F31i*c1*f3*u1*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F35*F30i*c0*f3*u0*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F34*F35i*c5*f3*u4*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F34*F34i*c4*f3*u4**2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F34*F33i*c3*f3*u3*u4*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F34*F32i*c2*f3*u2*u4*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F34*F31i*c1*f3*u1*u4*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F34*F30i*c0*f3*u0*u4*w1*w2
     &
      traza1 = traza1 + 16.D0*PC_q**3*q_uh*svp*svm*F33*F35i*c5*f3*u3*u5
     & *w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F33*F34i*c4*f3*u3*u4*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F33*F33i*c3*f3*u3**2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F33*F32i*c2*f3*u2*u3*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F33*F31i*c1*f3*u1*u3*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F33*F30i*c0*f3*u0*u3*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F32*F35i*c5*f3*u2*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F32*F34i*c4*f3*u2*u4*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F32*F33i*c3*f3*u2*u3*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F32*F32i*c2*f3*u2**2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F32*F31i*c1*f3*u1*u2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F32*F30i*c0*f3*u0*u2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F31*F35i*c5*f3*u1*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F31*F34i*c4*f3*u1*u4*w1*w2
     &
      traza1 = traza1 + 16.D0*PC_q**3*q_uh*svp*svm*F31*F33i*c3*f3*u1*u3
     & *w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F31*F32i*c2*f3*u1*u2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F31*F31i*c1*f3*u1**2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F31*F30i*c0*f3*u0*u1*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F30*F35i*c5*f3*u0*u5*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F30*F34i*c4*f3*u0*u4*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F30*F33i*c3*f3*u0*u3*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F30*F32i*c2*f3*u0*u2*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F30*F31i*c1*f3*u0*u1*w1*w2
     &  + 16.D0*PC_q**3*q_uh*svp*svm*F30*F30i*c0*f3*u0**2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F35*F35r*c5*f3*u5**2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F35*F34r*c4*f3*u4*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F35*F33r*c3*f3*u3*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F35*F32r*c2*f3*u2*u5*w1*w2
     &
      traza1 = traza1 - 16.D0*PC_q**3*q_uh*i_*svp*svm*F35*F31r*c1*f3*u1
     & *u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F35*F30r*c0*f3*u0*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F34*F35r*c5*f3*u4*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F34*F34r*c4*f3*u4**2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F34*F33r*c3*f3*u3*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F34*F32r*c2*f3*u2*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F34*F31r*c1*f3*u1*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F34*F30r*c0*f3*u0*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F33*F35r*c5*f3*u3*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F33*F34r*c4*f3*u3*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F33*F33r*c3*f3*u3**2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F33*F32r*c2*f3*u2*u3*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F33*F31r*c1*f3*u1*u3*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F33*F30r*c0*f3*u0*u3*w1*w2
     &
      traza1 = traza1 - 16.D0*PC_q**3*q_uh*i_*svp*svm*F32*F35r*c5*f3*u2
     & *u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F32*F34r*c4*f3*u2*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F32*F33r*c3*f3*u2*u3*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F32*F32r*c2*f3*u2**2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F32*F31r*c1*f3*u1*u2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F32*F30r*c0*f3*u0*u2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F31*F35r*c5*f3*u1*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F31*F34r*c4*f3*u1*u4*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F31*F33r*c3*f3*u1*u3*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F31*F32r*c2*f3*u1*u2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F31*F31r*c1*f3*u1**2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F31*F30r*c0*f3*u0*u1*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F30*F35r*c5*f3*u0*u5*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F30*F34r*c4*f3*u0*u4*w1*w2
     &
      traza1 = traza1 - 16.D0*PC_q**3*q_uh*i_*svp*svm*F30*F33r*c3*f3*u0
     & *u3*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F30*F32r*c2*f3*u0*u2*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F30*F31r*c1*f3*u0*u1*w1*w2
     &  - 16.D0*PC_q**3*q_uh*i_*svp*svm*F30*F30r*c0*f3*u0**2*w1*w2
     &  + 4.D0*PC_q**3*svm*dsvp*F45*F45r*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F45*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F45*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F45*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F45*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F45*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F45*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F45*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F45*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F45*F41r*c1*f4*u1*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dsvp*F45*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F45*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F44*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F44*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F44*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F44*F44r*c4*f4*u4**2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F44*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F44*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F44*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F44*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F44*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F44*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F44*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F44*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F43*F45r*c5*f4*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svm*dsvp*F43*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F43*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F43*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F43*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F43*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F43*F42r*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F43*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F43*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F43*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F43*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F43*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F42*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F42*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F42*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F42*F44r*c4*f4*u2*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dsvp*F42*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F42*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F42*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F42*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F42*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F42*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F42*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F42*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F41*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F41*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F41*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F41*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F41*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F41*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F41*F42r*c2*f4*u1*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svm*dsvp*F41*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F41*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F41*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F41*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F41*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F40*F45r*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F40*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F40*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F40*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F40*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F40*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F40*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F40*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*svm*dsvp*F40*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F40*F41r*c1*f4*u0*u1*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dsvp*F40*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**3*svm*dsvp*F40*F40r*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**3*svm*dssp*F45*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F45*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F45*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F45*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F45*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F45*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F44*F35r*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F44*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F44*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F44*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F44*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F44*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F43*F35r*c5*f3*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dssp*F43*F34r*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F43*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F43*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F43*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F43*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F42*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F42*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F42*F33r*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F42*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F42*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F42*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F41*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F41*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F41*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F41*F32r*c2*f3*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dssp*F41*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F41*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*svm*dssp*F40*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F40*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F40*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F40*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F40*F31r*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*svm*dssp*F40*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F35*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F35*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F35*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F35*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F35*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F35*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F34*F45r*c5*f4*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dssp*F34*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F34*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F34*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F34*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F34*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F33*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F33*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F33*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F33*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F33*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F33*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F32*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F32*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F32*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F32*F42r*c2*f4*u2**2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*svm*dssp*F32*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F32*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F31*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F31*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F31*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F31*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F31*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F31*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*svm*dssp*F30*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*svm*dssp*F30*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*svm*dssp*F30*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*svm*dssp*F30*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q**3*svm*dssp*F30*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*svm*dssp*F30*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC_q**3*svp*dsvm*F45*F45r*c5*f4*u5**2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dsvm*F45*F45r*c5*f4*u5**2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F45*F44r*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F45*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F45*F43r*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F45*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F45*F42r*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F45*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F45*F41r*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F45*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F45*F40r*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F45*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F44*F45r*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F44*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F44*F44r*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F44*F44r*c4*f4*u4**2*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svp*dsvm*F44*F43r*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F44*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F44*F42r*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F44*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F44*F41r*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F44*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F44*F40r*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F44*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F43*F45r*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F43*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F43*F44r*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F43*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F43*F43r*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F43*F43r*c3*f4*u3**2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F43*F42r*c2*f4*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dsvm*F43*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F43*F41r*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F43*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F43*F40r*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F43*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F42*F45r*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F42*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F42*F44r*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F42*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F42*F43r*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F42*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F42*F42r*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F42*F42r*c2*f4*u2**2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F42*F41r*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F42*F41r*c1*f4*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_q**3*svp*dsvm*F42*F40r*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F42*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F41*F45r*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F41*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F41*F44r*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F41*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F41*F43r*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F41*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F41*F42r*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F41*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F41*F41r*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F41*F41r*c1*f4*u1**2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F41*F40r*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F41*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F40*F45r*c5*f4*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dsvm*F40*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F40*F44r*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F40*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F40*F43r*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F40*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F40*F42r*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F40*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F40*F41r*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F40*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*svp*dsvm*F40*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**3*svp*dsvm*F40*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F45*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F45*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F45*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F45*F32r*c2*f3*u2*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dssm*F45*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F45*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F44*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F44*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F44*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F44*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F44*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F44*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F43*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F43*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F43*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F43*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F43*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F43*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F42*F35r*c5*f3*u2*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dssm*F42*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F42*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F42*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F42*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F42*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F41*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F41*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F41*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F41*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F41*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F41*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*svp*dssm*F40*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F40*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F40*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F40*F32r*c2*f3*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dssm*F40*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*svp*dssm*F40*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F35*F45r*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F35*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F35*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F35*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F35*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F35*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F34*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F34*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F34*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F34*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F34*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F34*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F33*F45r*c5*f4*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dssm*F33*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F33*F43r*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F33*F42r*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F33*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F33*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F32*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F32*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F32*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F32*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F32*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F32*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F31*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F31*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F31*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F31*F42r*c2*f4*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*svp*dssm*F31*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F31*F40r*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*svp*dssm*F30*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*svp*dssm*F30*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*svp*dssm*F30*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*svp*dssm*F30*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*svp*dssm*F30*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*svp*dssm*F30*F40r*c0*f4*u0**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F45*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F45*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F45*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F45*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F45*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F45*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F44*F35r*c5*f3*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*dsvp*F44*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F44*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F44*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F44*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F44*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F43*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F43*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F43*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F43*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F43*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F43*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F42*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F42*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F42*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F42*F32r*c2*f3*u2**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*dsvp*F42*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F42*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F41*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F41*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F41*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F41*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F41*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F41*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F40*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F40*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F40*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F40*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F40*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F40*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F35*F45r*c5*f4*u5**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*dsvp*F35*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F35*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F35*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F35*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F35*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F34*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F34*F44r*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F34*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F34*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F34*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F34*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F33*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F33*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F33*F43r*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F33*F42r*c2*f4*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*dsvp*F33*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F33*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F32*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F32*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F32*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F32*F42r*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F32*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F32*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F31*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F31*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F31*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F31*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F31*F41r*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F31*F40r*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F30*F45r*c5*f4*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*ssm*dsvp*F30*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F30*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F30*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F30*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*ssm*dsvp*F30*F40r*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**3*ssp*dsvm*F45*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F45*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F45*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F45*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F45*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F45*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F44*F35r*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F44*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F44*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F44*F32r*c2*f3*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*dsvm*F44*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F44*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F43*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F43*F34r*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F43*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F43*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F43*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F43*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F42*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F42*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F42*F33r*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F42*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F42*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F42*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F41*F35r*c5*f3*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*dsvm*F41*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F41*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F41*F32r*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F41*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F41*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F40*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F40*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F40*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F40*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F40*F31r*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F40*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F35*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F35*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F35*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F35*F42r*c2*f4*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*dsvm*F35*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F35*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F34*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F34*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F34*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F34*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F34*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F34*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F33*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F33*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F33*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F33*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F33*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F33*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F32*F45r*c5*f4*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*dsvm*F32*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F32*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F32*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F32*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F32*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F31*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F31*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F31*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F31*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F31*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F31*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F30*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F30*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F30*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F30*F42r*c2*f4*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*ssp*dsvm*F30*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*ssp*dsvm*F30*F40r*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F35*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F35*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F35*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F35*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F35*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F35*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F35*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F35*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F35*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F35*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F35*F30r*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F35*F30r*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F34*F35r*c5*f3*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svm*dsvp*F34*F35r*c5*f3*u4*u5*
     & w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F34*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F34*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F34*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F34*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F34*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F34*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F34*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F34*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F34*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F34*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F33*F35r*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F33*F35r*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F33*F34r*c4*f3*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svm*dsvp*F33*F34r*c4*f3*u3*u4*
     & w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F33*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F33*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F33*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F33*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F33*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F33*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F33*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F33*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F32*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F32*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F32*F34r*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F32*F34r*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F32*F33r*c3*f3*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svm*dsvp*F32*F33r*c3*f3*u2*u3*
     & w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F32*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F32*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F32*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F32*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F32*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F32*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F31*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F31*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F31*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F31*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F31*F33r*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F31*F33r*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F31*F32r*c2*f3*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svm*dsvp*F31*F32r*c2*f3*u1*u2*
     & w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F31*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F31*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F31*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F31*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F30*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F30*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F30*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F30*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F30*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F30*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F30*F32r*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F30*F32r*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F30*F31r*c1*f3*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svm*dsvp*F30*F31r*c1*f3*u0*u1*
     & w1
     &  - 4.D0*PC_q**3*qq*svm*dsvp*F30*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*qq*svm*dsvp*F30*F30r*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F35*F35r*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F35*F35r*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F35*F34r*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F35*F34r*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F35*F33r*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F35*F33r*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F35*F32r*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F35*F32r*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F35*F31r*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F35*F31r*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F35*F30r*c0*f3*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*dsvm*F35*F30r*c0*f3*u0*u5*
     & w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F34*F35r*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F34*F35r*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F34*F34r*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F34*F34r*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F34*F33r*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F34*F33r*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F34*F32r*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F34*F32r*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F34*F31r*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F34*F31r*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F34*F30r*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F34*F30r*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F33*F35r*c5*f3*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*dsvm*F33*F35r*c5*f3*u3*u5*
     & w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F33*F34r*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F33*F34r*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F33*F33r*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F33*F33r*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F33*F32r*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F33*F32r*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F33*F31r*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F33*F31r*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F33*F30r*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F33*F30r*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F32*F35r*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F32*F35r*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F32*F34r*c4*f3*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*dsvm*F32*F34r*c4*f3*u2*u4*
     & w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F32*F33r*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F32*F33r*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F32*F32r*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F32*F32r*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F32*F31r*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F32*F31r*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F32*F30r*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F32*F30r*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F31*F35r*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F31*F35r*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F31*F34r*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F31*F34r*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F31*F33r*c3*f3*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*dsvm*F31*F33r*c3*f3*u1*u3*
     & w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F31*F32r*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F31*F32r*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F31*F31r*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F31*F31r*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F31*F30r*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F31*F30r*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F30*F35r*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F30*F35r*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F30*F34r*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F30*F34r*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F30*F33r*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F30*F33r*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F30*F32r*c2*f3*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*qq*svp*dsvm*F30*F32r*c2*f3*u0*u2*
     & w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F30*F31r*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F30*F31r*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*qq*svp*dsvm*F30*F30r*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*qq*svp*dsvm*F30*F30r*c0*f3*u0**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F45*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F45*F45i*c5*f4*u5**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F45*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F45*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F45*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F45*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F45*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F45*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F45*F41i*c1*f4*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svm*dsvp*F45*F41i*c1*f4*u1*u5*
     & w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F45*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F45*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F44*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F44*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F44*F44i*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F44*F44i*c4*f4*u4**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F44*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F44*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F44*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F44*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F44*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F44*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F44*F40i*c0*f4*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svm*dsvp*F44*F40i*c0*f4*u0*u4*
     & w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F43*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F43*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F43*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F43*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F43*F43i*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F43*F43i*c3*f4*u3**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F43*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F43*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F43*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F43*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F43*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F43*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F42*F45i*c5*f4*u2*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svm*dsvp*F42*F45i*c5*f4*u2*u5*
     & w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F42*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F42*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F42*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F42*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F42*F42i*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F42*F42i*c2*f4*u2**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F42*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F42*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F42*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F42*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F41*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F41*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F41*F44i*c4*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svm*dsvp*F41*F44i*c4*f4*u1*u4*
     & w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F41*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F41*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F41*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F41*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F41*F41i*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F41*F41i*c1*f4*u1**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F41*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F41*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F40*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F40*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F40*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F40*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F40*F43i*c3*f4*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svm*dsvp*F40*F43i*c3*f4*u0*u3*
     & w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F40*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F40*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F40*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F40*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*i_*svm*dsvp*F40*F40i*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**3*i_*svm*dsvp*F40*F40i*c0*f4*u0**2*w1
     &  + 4.D0*PC_q**3*i_*svm*dssp*F45*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F45*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F45*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F45*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F45*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F45*F30i*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F44*F35i*c5*f3*u4*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svm*dssp*F44*F34i*c4*f3*u4**2*
     & w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F44*F33i*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F44*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F44*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F44*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F43*F35i*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F43*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F43*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F43*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F43*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F43*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F42*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F42*F34i*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F42*F33i*c3*f3*u2*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svm*dssp*F42*F32i*c2*f3*u2**2*
     & w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F42*F31i*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F42*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F41*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F41*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F41*F33i*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F41*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F41*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F41*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F40*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F40*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F40*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F40*F32i*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F40*F31i*c1*f3*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svm*dssp*F40*F30i*c0*f3*u0**2*
     & w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F35*F45i*c5*f4*u5**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F35*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F35*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F35*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F35*F41i*c1*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F35*F40i*c0*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F34*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F34*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F34*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F34*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F34*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F34*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F33*F45i*c5*f4*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svm*dssp*F33*F44i*c4*f4*u3*u4*
     & w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F33*F43i*c3*f4*u3**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F33*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F33*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F33*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F32*F45i*c5*f4*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F32*F44i*c4*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F32*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F32*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F32*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F32*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F31*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F31*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F31*F43i*c3*f4*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*svm*dssp*F31*F42i*c2*f4*u1*u2*
     & w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F31*F41i*c1*f4*u1**2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F31*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F30*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F30*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F30*F43i*c3*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F30*F42i*c2*f4*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F30*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*svm*dssp*F30*F40i*c0*f4*u0**2*w2
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F45*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F45*F45i*c5*f4*u5**2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F45*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F45*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F45*F43i*c3*f4*u3*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dsvm*F45*F43i*c3*f4*u3*u5*
     & w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F45*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F45*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F45*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F45*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F45*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F45*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F44*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F44*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F44*F44i*c4*f4*u4**2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F44*F44i*c4*f4*u4**2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F44*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F44*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F44*F42i*c2*f4*u2*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dsvm*F44*F42i*c2*f4*u2*u4*
     & w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F44*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F44*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F44*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F44*F40i*c0*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F43*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F43*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F43*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F43*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F43*F43i*c3*f4*u3**2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F43*F43i*c3*f4*u3**2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F43*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F43*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F43*F41i*c1*f4*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dsvm*F43*F41i*c1*f4*u1*u3*
     & w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F43*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F43*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F42*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F42*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F42*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F42*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F42*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F42*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F42*F42i*c2*f4*u2**2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F42*F42i*c2*f4*u2**2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F42*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F42*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F42*F40i*c0*f4*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dsvm*F42*F40i*c0*f4*u0*u2*
     & w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F41*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F41*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F41*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F41*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F41*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F41*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F41*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F41*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F41*F41i*c1*f4*u1**2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F41*F41i*c1*f4*u1**2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F41*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F41*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F40*F45i*c5*f4*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dsvm*F40*F45i*c5*f4*u0*u5*
     & w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F40*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F40*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F40*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F40*F43i*c3*f4*u0*u3*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F40*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F40*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F40*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F40*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*PC_q**3*i_*svp*dsvm*F40*F40i*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**3*i_*svp*dsvm*F40*F40i*c0*f4*u0**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F45*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F45*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F45*F33i*c3*f3*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dssm*F45*F32i*c2*f3*u2*u5*
     & w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F45*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F45*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F44*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F44*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F44*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F44*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F44*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F44*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F43*F35i*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F43*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F43*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F43*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F43*F31i*c1*f3*u1*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dssm*F43*F30i*c0*f3*u0*u3*
     & w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F42*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F42*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F42*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F42*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F42*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F42*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F41*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F41*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F41*F33i*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F41*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F41*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F41*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F40*F35i*c5*f3*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dssm*F40*F34i*c4*f3*u0*u4*
     & w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F40*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F40*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F40*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F40*F30i*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F35*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F35*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F35*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F35*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F35*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F35*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F34*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F34*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F34*F43i*c3*f4*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dssm*F34*F42i*c2*f4*u2*u4*
     & w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F34*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F34*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F33*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F33*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F33*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F33*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F33*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F33*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F32*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F32*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F32*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F32*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F32*F41i*c1*f4*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*svp*dssm*F32*F40i*c0*f4*u0*u2*
     & w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F31*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F31*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F31*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F31*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F31*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F31*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F30*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F30*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F30*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F30*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F30*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*svp*dssm*F30*F40i*c0*f4*u0**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F45*F35i*c5*f3*u5**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*dsvp*F45*F34i*c4*f3*u4*u5*
     & w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F45*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F45*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F45*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F45*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F44*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F44*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F44*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F44*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F44*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F44*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F43*F35i*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F43*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F43*F33i*c3*f3*u3**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*dsvp*F43*F32i*c2*f3*u2*u3*
     & w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F43*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F43*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F42*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F42*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F42*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F42*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F42*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F42*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F41*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F41*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F41*F33i*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F41*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F41*F31i*c1*f3*u1**2*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*dsvp*F41*F30i*c0*f3*u0*u1*
     & w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F40*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F40*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F40*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F40*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F40*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F40*F30i*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F35*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F35*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F35*F43i*c3*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F35*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F35*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F35*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F34*F45i*c5*f4*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*dsvp*F34*F44i*c4*f4*u4**2*
     & w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F34*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F34*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F34*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F34*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F33*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F33*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F33*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F33*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F33*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F33*F40i*c0*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F32*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F32*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F32*F43i*c3*f4*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*dsvp*F32*F42i*c2*f4*u2**2*
     & w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F32*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F32*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F31*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F31*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F31*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F31*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F31*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F31*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F30*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F30*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F30*F43i*c3*f4*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F30*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*ssm*dsvp*F30*F41i*c1*f4*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*PC_q**3*i_*ssm*dsvp*F30*F40i*c0*f4*u0**2*
     & w1
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F45*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F45*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F45*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F45*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F45*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F45*F30i*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F44*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F44*F34i*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F44*F33i*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F44*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F44*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F44*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F43*F35i*c5*f3*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*dsvm*F43*F34i*c4*f3*u3*u4*
     & w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F43*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F43*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F43*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F43*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F42*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F42*F34i*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F42*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F42*F32i*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F42*F31i*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F42*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F41*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F41*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F41*F33i*c3*f3*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*dsvm*F41*F32i*c2*f3*u1*u2*
     & w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F41*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F41*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F40*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F40*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F40*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F40*F32i*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F40*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F40*F30i*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F35*F45i*c5*f4*u5**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F35*F44i*c4*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F35*F43i*c3*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F35*F42i*c2*f4*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F35*F41i*c1*f4*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*dsvm*F35*F40i*c0*f4*u0*u5*
     & w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F34*F45i*c5*f4*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F34*F44i*c4*f4*u4**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F34*F43i*c3*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F34*F42i*c2*f4*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F34*F41i*c1*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F34*F40i*c0*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F33*F45i*c5*f4*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F33*F44i*c4*f4*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F33*F43i*c3*f4*u3**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F33*F42i*c2*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F33*F41i*c1*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F33*F40i*c0*f4*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F32*F45i*c5*f4*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*dsvm*F32*F44i*c4*f4*u2*u4*
     & w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F32*F43i*c3*f4*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F32*F42i*c2*f4*u2**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F32*F41i*c1*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F32*F40i*c0*f4*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F31*F45i*c5*f4*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F31*F44i*c4*f4*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F31*F43i*c3*f4*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F31*F42i*c2*f4*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F31*F41i*c1*f4*u1**2*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F31*F40i*c0*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F30*F45i*c5*f4*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F30*F44i*c4*f4*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F30*F43i*c3*f4*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*ssp*dsvm*F30*F42i*c2*f4*u0*u2*
     & w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F30*F41i*c1*f4*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*ssp*dsvm*F30*F40i*c0*f4*u0**2*w2
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F31i*c1*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F31i*c1*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F30i*c0*f3*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svm*dsvp*F35*F30i*c0*f3*u0*
     & u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F34i*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F33i*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F30i*c0*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F34*F30i*c0*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F35i*c5*f3*u3*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F35i*c5*f3*u3*
     & u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F33*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F35i*c5*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F35i*c5*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F34i*c4*f3*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F34i*c4*f3*u2*
     & u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F32i*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F31i*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F32*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F34i*c4*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F34i*c4*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F33i*c3*f3*u1*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F33i*c3*f3*u1*
     & u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F31*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F33i*c3*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F33i*c3*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F32i*c2*f3*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F32i*c2*f3*u0*
     & u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F30i*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svm*dsvp*F30*F30i*c0*f3*u0**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F35i*c5*f3*u5**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F35i*c5*f3*u5**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F34i*c4*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F34i*c4*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F33i*c3*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F33i*c3*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F32i*c2*f3*u2*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F32i*c2*f3*u2*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F31i*c1*f3*u1*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F31i*c1*f3*u1*
     & u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F30i*c0*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F35*F30i*c0*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F35i*c5*f3*u4*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F35i*c5*f3*u4*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F34i*c4*f3*u4**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F34i*c4*f3*u4**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F33i*c3*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F33i*c3*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F32i*c2*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F32i*c2*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F31i*c1*f3*u1*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F31i*c1*f3*u1*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F30i*c0*f3*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*dsvm*F34*F30i*c0*f3*u0*
     & u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F35i*c5*f3*u3*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F35i*c5*f3*u3*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F34i*c4*f3*u3*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F34i*c4*f3*u3*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F33i*c3*f3*u3**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F33i*c3*f3*u3**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F32i*c2*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F32i*c2*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F31i*c1*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F31i*c1*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F30i*c0*f3*u0*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F33*F30i*c0*f3*u0*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F35i*c5*f3*u2*u5*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F35i*c5*f3*u2*
     & u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F34i*c4*f3*u2*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F34i*c4*f3*u2*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F33i*c3*f3*u2*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F33i*c3*f3*u2*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F32i*c2*f3*u2**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F32i*c2*f3*u2**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F31i*c1*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F31i*c1*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F30i*c0*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F32*F30i*c0*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F35i*c5*f3*u1*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F35i*c5*f3*u1*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F34i*c4*f3*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F34i*c4*f3*u1*
     & u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F33i*c3*f3*u1*u3*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F33i*c3*f3*u1*u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F32i*c2*f3*u1*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F32i*c2*f3*u1*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F31i*c1*f3*u1**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F31i*c1*f3*u1**2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F30i*c0*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F31*F30i*c0*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F35i*c5*f3*u0*u5*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F35i*c5*f3*u0*u5*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F34i*c4*f3*u0*u4*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F34i*c4*f3*u0*u4*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F33i*c3*f3*u0*u3*w2
     &
      traza1 = traza1 + 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F33i*c3*f3*u0*
     & u3*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F32i*c2*f3*u0*u2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F32i*c2*f3*u0*u2*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F31i*c1*f3*u0*u1*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F31i*c1*f3*u0*u1*w1
     &  - 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F30i*c0*f3*u0**2*w2
     &  + 4.D0*PC_q**3*i_*qq*svp*dsvm*F30*F30i*c0*f3*u0**2*w1
     &  - 8.D0*PC_q**4*svm*dsvp*F35*F35r*c5*f3*u5**2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F35*F34r*c4*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F35*F33r*c3*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F35*F32r*c2*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F35*F31r*c1*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F35*F30r*c0*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F34*F35r*c5*f3*u4*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svm*dsvp*F34*F34r*c4*f3*u4**2*w1*
     & w2
     &  - 8.D0*PC_q**4*svm*dsvp*F34*F33r*c3*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F34*F32r*c2*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F34*F31r*c1*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F34*F30r*c0*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F33*F35r*c5*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F33*F34r*c4*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F33*F33r*c3*f3*u3**2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F33*F32r*c2*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F33*F31r*c1*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F33*F30r*c0*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F32*F35r*c5*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F32*F34r*c4*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F32*F33r*c3*f3*u2*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svm*dsvp*F32*F32r*c2*f3*u2**2*w1*
     & w2
     &  - 8.D0*PC_q**4*svm*dsvp*F32*F31r*c1*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F32*F30r*c0*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F31*F35r*c5*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F31*F34r*c4*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F31*F33r*c3*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F31*F32r*c2*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F31*F31r*c1*f3*u1**2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F31*F30r*c0*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F30*F35r*c5*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F30*F34r*c4*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F30*F33r*c3*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F30*F32r*c2*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*svm*dsvp*F30*F31r*c1*f3*u0*u1*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svm*dsvp*F30*F30r*c0*f3*u0**2*w1*
     & w2
     &  - 8.D0*PC_q**4*svp*dsvm*F35*F35r*c5*f3*u5**2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F35*F34r*c4*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F35*F33r*c3*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F35*F32r*c2*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F35*F31r*c1*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F35*F30r*c0*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F34*F35r*c5*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F34*F34r*c4*f3*u4**2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F34*F33r*c3*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F34*F32r*c2*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F34*F31r*c1*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F34*F30r*c0*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F33*F35r*c5*f3*u3*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svp*dsvm*F33*F34r*c4*f3*u3*u4*w1*
     & w2
     &  - 8.D0*PC_q**4*svp*dsvm*F33*F33r*c3*f3*u3**2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F33*F32r*c2*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F33*F31r*c1*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F33*F30r*c0*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F32*F35r*c5*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F32*F34r*c4*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F32*F33r*c3*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F32*F32r*c2*f3*u2**2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F32*F31r*c1*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F32*F30r*c0*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F31*F35r*c5*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F31*F34r*c4*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F31*F33r*c3*f3*u1*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*svp*dsvm*F31*F32r*c2*f3*u1*u2*w1*
     & w2
     &  - 8.D0*PC_q**4*svp*dsvm*F31*F31r*c1*f3*u1**2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F31*F30r*c0*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F30*F35r*c5*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F30*F34r*c4*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F30*F33r*c3*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F30*F32r*c2*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F30*F31r*c1*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*svp*dsvm*F30*F30r*c0*f3*u0**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F35*F35i*c5*f3*u5**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F35*F34i*c4*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F35*F33i*c3*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F35*F32i*c2*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F35*F31i*c1*f3*u1*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svm*dsvp*F35*F30i*c0*f3*u0*u5*
     & w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F34*F35i*c5*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F34*F34i*c4*f3*u4**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F34*F33i*c3*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F34*F32i*c2*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F34*F31i*c1*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F34*F30i*c0*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F33*F35i*c5*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F33*F34i*c4*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F33*F33i*c3*f3*u3**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F33*F32i*c2*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F33*F31i*c1*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F33*F30i*c0*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F32*F35i*c5*f3*u2*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svm*dsvp*F32*F34i*c4*f3*u2*u4*
     & w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F32*F33i*c3*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F32*F32i*c2*f3*u2**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F32*F31i*c1*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F32*F30i*c0*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F31*F35i*c5*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F31*F34i*c4*f3*u1*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F31*F33i*c3*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F31*F32i*c2*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F31*F31i*c1*f3*u1**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F31*F30i*c0*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F30*F35i*c5*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F30*F34i*c4*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F30*F33i*c3*f3*u0*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svm*dsvp*F30*F32i*c2*f3*u0*u2*
     & w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F30*F31i*c1*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*i_*svm*dsvp*F30*F30i*c0*f3*u0**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F35*F35i*c5*f3*u5**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F35*F34i*c4*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F35*F33i*c3*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F35*F32i*c2*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F35*F31i*c1*f3*u1*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F35*F30i*c0*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F34*F35i*c5*f3*u4*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F34*F34i*c4*f3*u4**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F34*F33i*c3*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F34*F32i*c2*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F34*F31i*c1*f3*u1*u4*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svp*dsvm*F34*F30i*c0*f3*u0*u4*
     & w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F33*F35i*c5*f3*u3*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F33*F34i*c4*f3*u3*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F33*F33i*c3*f3*u3**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F33*F32i*c2*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F33*F31i*c1*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F33*F30i*c0*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F32*F35i*c5*f3*u2*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F32*F34i*c4*f3*u2*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F32*F33i*c3*f3*u2*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F32*F32i*c2*f3*u2**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F32*F31i*c1*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F32*F30i*c0*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F31*F35i*c5*f3*u1*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_q**4*i_*svp*dsvm*F31*F34i*c4*f3*u1*u4*
     & w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F31*F33i*c3*f3*u1*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F31*F32i*c2*f3*u1*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F31*F31i*c1*f3*u1**2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F31*F30i*c0*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F30*F35i*c5*f3*u0*u5*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F30*F34i*c4*f3*u0*u4*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F30*F33i*c3*f3*u0*u3*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F30*F32i*c2*f3*u0*u2*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F30*F31i*c1*f3*u0*u1*w1*w2
     &  - 8.D0*PC_q**4*i_*svp*dsvm*F30*F30i*c0*f3*u0**2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F15*F15i*c5*f1*u5**2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F15*F14i*c4*f1*u4*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F15*F13i*c3*f1*u3*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*svp*svm*F15*F12i*c2*f1*u2*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F15*F11i*c1*f1*u1*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F15*F10i*c0*f1*u0*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F14*F15i*c5*f1*u4*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F14*F14i*c4*f1*u4**2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F14*F13i*c3*f1*u3*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F14*F12i*c2*f1*u2*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F14*F11i*c1*f1*u1*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F14*F10i*c0*f1*u0*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F13*F15i*c5*f1*u3*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F13*F14i*c4*f1*u3*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F13*F13i*c3*f1*u3**2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F13*F12i*c2*f1*u2*u3*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F13*F11i*c1*f1*u1*u3*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F13*F10i*c0*f1*u0*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*svp*svm*F12*F15i*c5*f1*u2*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F12*F14i*c4*f1*u2*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F12*F13i*c3*f1*u2*u3*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F12*F12i*c2*f1*u2**2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F12*F11i*c1*f1*u1*u2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F12*F10i*c0*f1*u0*u2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F11*F15i*c5*f1*u1*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F11*F14i*c4*f1*u1*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F11*F13i*c3*f1*u1*u3*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F11*F12i*c2*f1*u1*u2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F11*F11i*c1*f1*u1**2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F11*F10i*c0*f1*u0*u1*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F10*F15i*c5*f1*u0*u5*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F10*F14i*c4*f1*u0*u4*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F10*F13i*c3*f1*u0*u3*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*svp*svm*F10*F12i*c2*f1*u0*u2*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F10*F11i*c1*f1*u0*u1*w1*w2
     &  - 8.D0*PC_uh*svp*svm*F10*F10i*c0*f1*u0**2*w1*w2
     &  + 4.D0*PC_uh*ssm*svp*F25*F15i*c5*f1*u5**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F25*F14i*c4*f1*u4*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F25*F13i*c3*f1*u3*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F25*F12i*c2*f1*u2*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F25*F11i*c1*f1*u1*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F25*F10i*c0*f1*u0*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F24*F15i*c5*f1*u4*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F24*F14i*c4*f1*u4**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F24*F13i*c3*f1*u3*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F24*F12i*c2*f1*u2*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F24*F11i*c1*f1*u1*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F24*F10i*c0*f1*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*ssm*svp*F23*F15i*c5*f1*u3*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F23*F14i*c4*f1*u3*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F23*F13i*c3*f1*u3**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F23*F12i*c2*f1*u2*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F23*F11i*c1*f1*u1*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F23*F10i*c0*f1*u0*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F22*F15i*c5*f1*u2*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F22*F14i*c4*f1*u2*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F22*F13i*c3*f1*u2*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F22*F12i*c2*f1*u2**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F22*F11i*c1*f1*u1*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F22*F10i*c0*f1*u0*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F21*F15i*c5*f1*u1*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F21*F14i*c4*f1*u1*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F21*F13i*c3*f1*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*ssm*svp*F21*F12i*c2*f1*u1*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F21*F11i*c1*f1*u1**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F21*F10i*c0*f1*u0*u1*w1
     &  + 4.D0*PC_uh*ssm*svp*F20*F15i*c5*f1*u0*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F20*F14i*c4*f1*u0*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F20*F13i*c3*f1*u0*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F20*F12i*c2*f1*u0*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F20*F11i*c1*f1*u0*u1*w1
     &  + 4.D0*PC_uh*ssm*svp*F20*F10i*c0*f1*u0**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F15*F25i*c5*f2*u5**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F15*F24i*c4*f2*u4*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F15*F23i*c3*f2*u3*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F15*F22i*c2*f2*u2*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F15*F21i*c1*f2*u1*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F15*F20i*c0*f2*u0*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*ssm*svp*F14*F25i*c5*f2*u4*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F14*F24i*c4*f2*u4**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F14*F23i*c3*f2*u3*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F14*F22i*c2*f2*u2*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F14*F21i*c1*f2*u1*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F14*F20i*c0*f2*u0*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*F25i*c5*f2*u3*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*F24i*c4*f2*u3*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*F23i*c3*f2*u3**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*F22i*c2*f2*u2*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*F21i*c1*f2*u1*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F13*F20i*c0*f2*u0*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F12*F25i*c5*f2*u2*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F12*F24i*c4*f2*u2*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F12*F23i*c3*f2*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*ssm*svp*F12*F22i*c2*f2*u2**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F12*F21i*c1*f2*u1*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F12*F20i*c0*f2*u0*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*F25i*c5*f2*u1*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*F24i*c4*f2*u1*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*F23i*c3*f2*u1*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*F22i*c2*f2*u1*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*F21i*c1*f2*u1**2*w1
     &  + 4.D0*PC_uh*ssm*svp*F11*F20i*c0*f2*u0*u1*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*F25i*c5*f2*u0*u5*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*F24i*c4*f2*u0*u4*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*F23i*c3*f2*u0*u3*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*F22i*c2*f2*u0*u2*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*F21i*c1*f2*u0*u1*w1
     &  + 4.D0*PC_uh*ssm*svp*F10*F20i*c0*f2*u0**2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*ssp*svm*F25*F15i*c5*f1*u5**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F25*F14i*c4*f1*u4*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F25*F13i*c3*f1*u3*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F25*F12i*c2*f1*u2*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F25*F11i*c1*f1*u1*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F25*F10i*c0*f1*u0*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F24*F15i*c5*f1*u4*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F24*F14i*c4*f1*u4**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F24*F13i*c3*f1*u3*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F24*F12i*c2*f1*u2*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F24*F11i*c1*f1*u1*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F24*F10i*c0*f1*u0*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F23*F15i*c5*f1*u3*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F23*F14i*c4*f1*u3*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F23*F13i*c3*f1*u3**2*w2
     &
      traza1 = traza1 + 4.D0*PC_uh*ssp*svm*F23*F12i*c2*f1*u2*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F23*F11i*c1*f1*u1*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F23*F10i*c0*f1*u0*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F22*F15i*c5*f1*u2*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F22*F14i*c4*f1*u2*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F22*F13i*c3*f1*u2*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F22*F12i*c2*f1*u2**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F22*F11i*c1*f1*u1*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F22*F10i*c0*f1*u0*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F21*F15i*c5*f1*u1*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F21*F14i*c4*f1*u1*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F21*F13i*c3*f1*u1*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F21*F12i*c2*f1*u1*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F21*F11i*c1*f1*u1**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F21*F10i*c0*f1*u0*u1*w2
     &
      traza1 = traza1 + 4.D0*PC_uh*ssp*svm*F20*F15i*c5*f1*u0*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F20*F14i*c4*f1*u0*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F20*F13i*c3*f1*u0*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F20*F12i*c2*f1*u0*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F20*F11i*c1*f1*u0*u1*w2
     &  + 4.D0*PC_uh*ssp*svm*F20*F10i*c0*f1*u0**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F15*F25i*c5*f2*u5**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F15*F24i*c4*f2*u4*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F15*F23i*c3*f2*u3*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F15*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F15*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F15*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F14*F25i*c5*f2*u4*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F14*F24i*c4*f2*u4**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F14*F23i*c3*f2*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*PC_uh*ssp*svm*F14*F22i*c2*f2*u2*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F14*F21i*c1*f2*u1*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F14*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*F24i*c4*f2*u3*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*F23i*c3*f2*u3**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*F22i*c2*f2*u2*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*F21i*c1*f2*u1*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F13*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F12*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F12*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F12*F23i*c3*f2*u2*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F12*F22i*c2*f2*u2**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F12*F21i*c1*f2*u1*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F12*F20i*c0*f2*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*PC_uh*ssp*svm*F11*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F11*F24i*c4*f2*u1*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F11*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F11*F22i*c2*f2*u1*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F11*F21i*c1*f2*u1**2*w2
     &  + 4.D0*PC_uh*ssp*svm*F11*F20i*c0*f2*u0*u1*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*F25i*c5*f2*u0*u5*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*F21i*c1*f2*u0*u1*w2
     &  + 4.D0*PC_uh*ssp*svm*F10*F20i*c0*f2*u0**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F15i*c5*f1*u5**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F15i*c5*f1*u5**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F14i*c4*f1*u4*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F45*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F13i*c3*f1*u3*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F12i*c2*f1*u2*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F11i*c1*f1*u1*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F10i*c0*f1*u0*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F45*F10i*c0*f1*u0*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F15i*c5*f1*u4*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F14i*c4*f1*u4**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F14i*c4*f1*u4**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F13i*c3*f1*u3*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F13i*c3*f1*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F44*F12i*c2*f1*u2*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F11i*c1*f1*u1*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F11i*c1*f1*u1*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F10i*c0*f1*u0*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F44*F10i*c0*f1*u0*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F15i*c5*f1*u3*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F15i*c5*f1*u3*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F14i*c4*f1*u3*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F13i*c3*f1*u3**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F13i*c3*f1*u3**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F12i*c2*f1*u2*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F11i*c1*f1*u1*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F43*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F10i*c0*f1*u0*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F43*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F15i*c5*f1*u2*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F14i*c4*f1*u2*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F13i*c3*f1*u2*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F13i*c3*f1*u2*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F12i*c2*f1*u2**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F12i*c2*f1*u2**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F11i*c1*f1*u1*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F10i*c0*f1*u0*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F42*F10i*c0*f1*u0*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F41*F15i*c5*f1*u1*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F15i*c5*f1*u1*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F14i*c4*f1*u1*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F13i*c3*f1*u1*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F13i*c3*f1*u1*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F12i*c2*f1*u1*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F11i*c1*f1*u1**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F11i*c1*f1*u1**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F10i*c0*f1*u0*u1*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F41*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F15i*c5*f1*u0*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F14i*c4*f1*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F40*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F13i*c3*f1*u0*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F12i*c2*f1*u0*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F11i*c1*f1*u0*u1*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F10i*c0*f1*u0**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F40*F10i*c0*f1*u0**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F45i*c5*f4*u5**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F45i*c5*f4*u5**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F44i*c4*f4*u4*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F43i*c3*f4*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F15*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F42i*c2*f4*u2*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F41i*c1*f4*u1*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F15*F40i*c0*f4*u0*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F45i*c5*f4*u4*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F44i*c4*f4*u4**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F44i*c4*f4*u4**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F43i*c3*f4*u3*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F42i*c2*f4*u2*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F41i*c1*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F14*F41i*c1*f4*u1*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F14*F40i*c0*f4*u0*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F45i*c5*f4*u3*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F44i*c4*f4*u3*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F43i*c3*f4*u3**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F43i*c3*f4*u3**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F42i*c2*f4*u2*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F41i*c1*f4*u1*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F13*F40i*c0*f4*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F12*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F45i*c5*f4*u2*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F44i*c4*f4*u2*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F43i*c3*f4*u2*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F42i*c2*f4*u2**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F42i*c2*f4*u2**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F41i*c1*f4*u1*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F12*F40i*c0*f4*u0*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F45i*c5*f4*u1*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F44i*c4*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F11*F44i*c4*f4*u1*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F43i*c3*f4*u1*u3*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F42i*c2*f4*u1*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F41i*c1*f4*u1**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F41i*c1*f4*u1**2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F11*F40i*c0*f4*u0*u1*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F45i*c5*f4*u0*u5*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F44i*c4*f4*u0*u4*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F43i*c3*f4*u0*u3*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*qq*svp*svm*F10*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F42i*c2*f4*u0*u2*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F41i*c1*f4*u0*u1*w1
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F40i*c0*f4*u0**2*w2
     &  - 4.D0*PC_uh*qq*svp*svm*F10*F40i*c0*f4*u0**2*w1
     &  + 8.D0*PC_uh*PC2*svp*svm*F25*F25i*c5*f2*u5**2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F25*F24i*c4*f2*u4*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F25*F23i*c3*f2*u3*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F25*F22i*c2*f2*u2*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F25*F21i*c1*f2*u1*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F25*F20i*c0*f2*u0*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F24*F25i*c5*f2*u4*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F24*F24i*c4*f2*u4**2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F24*F23i*c3*f2*u3*u4*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*PC2*svp*svm*F24*F22i*c2*f2*u2*u4*w1*
     & w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F24*F21i*c1*f2*u1*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F24*F20i*c0*f2*u0*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F23*F25i*c5*f2*u3*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F23*F24i*c4*f2*u3*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F23*F23i*c3*f2*u3**2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F23*F22i*c2*f2*u2*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F23*F21i*c1*f2*u1*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F23*F20i*c0*f2*u0*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F22*F25i*c5*f2*u2*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F22*F24i*c4*f2*u2*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F22*F23i*c3*f2*u2*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F22*F22i*c2*f2*u2**2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F22*F21i*c1*f2*u1*u2*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*PC2*svp*svm*F22*F20i*c0*f2*u0*u2*w1*
     & w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F21*F25i*c5*f2*u1*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F21*F24i*c4*f2*u1*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F21*F23i*c3*f2*u1*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F21*F22i*c2*f2*u1*u2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F21*F21i*c1*f2*u1**2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F21*F20i*c0*f2*u0*u1*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F20*F25i*c5*f2*u0*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F20*F24i*c4*f2*u0*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F20*F23i*c3*f2*u0*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F20*F22i*c2*f2*u0*u2*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F20*F21i*c1*f2*u0*u1*w1*w2
     &  + 8.D0*PC_uh*PC2*svp*svm*F20*F20i*c0*f2*u0**2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F45*F45i*c5*f4*u5**2*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*PC2*qq*svp*svm*F45*F44i*c4*f4*u4*u5*
     & w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F45*F43i*c3*f4*u3*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F45*F42i*c2*f4*u2*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F45*F41i*c1*f4*u1*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F45*F40i*c0*f4*u0*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F44*F45i*c5*f4*u4*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F44*F44i*c4*f4*u4**2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F44*F43i*c3*f4*u3*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F44*F42i*c2*f4*u2*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F44*F41i*c1*f4*u1*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F44*F40i*c0*f4*u0*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F43*F45i*c5*f4*u3*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F43*F44i*c4*f4*u3*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F43*F43i*c3*f4*u3**2*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*PC2*qq*svp*svm*F43*F42i*c2*f4*u2*u3*
     & w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F43*F41i*c1*f4*u1*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F43*F40i*c0*f4*u0*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F42*F45i*c5*f4*u2*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F42*F44i*c4*f4*u2*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F42*F43i*c3*f4*u2*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F42*F42i*c2*f4*u2**2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F42*F41i*c1*f4*u1*u2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F42*F40i*c0*f4*u0*u2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F41*F45i*c5*f4*u1*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F41*F44i*c4*f4*u1*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F41*F43i*c3*f4*u1*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F41*F42i*c2*f4*u1*u2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F41*F41i*c1*f4*u1**2*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*PC2*qq*svp*svm*F41*F40i*c0*f4*u0*u1*
     & w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F40*F45i*c5*f4*u0*u5*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F40*F44i*c4*f4*u0*u4*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F40*F43i*c3*f4*u0*u3*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F40*F42i*c2*f4*u0*u2*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F40*F41i*c1*f4*u0*u1*w1*w2
     &  + 8.D0*PC_uh*PC2*qq*svp*svm*F40*F40i*c0*f4*u0**2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F15*F15r*c5*f1*u5**2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F15*F14r*c4*f1*u4*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F15*F13r*c3*f1*u3*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F15*F12r*c2*f1*u2*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F15*F11r*c1*f1*u1*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F15*F10r*c0*f1*u0*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F14*F15r*c5*f1*u4*u5*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*i_*svp*svm*F14*F14r*c4*f1*u4**2*w1*
     & w2
     &  + 8.D0*PC_uh*i_*svp*svm*F14*F13r*c3*f1*u3*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F14*F12r*c2*f1*u2*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F14*F11r*c1*f1*u1*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F14*F10r*c0*f1*u0*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F13*F15r*c5*f1*u3*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F13*F14r*c4*f1*u3*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F13*F13r*c3*f1*u3**2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F13*F12r*c2*f1*u2*u3*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F13*F11r*c1*f1*u1*u3*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F13*F10r*c0*f1*u0*u3*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F12*F15r*c5*f1*u2*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F12*F14r*c4*f1*u2*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F12*F13r*c3*f1*u2*u3*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*i_*svp*svm*F12*F12r*c2*f1*u2**2*w1*
     & w2
     &  + 8.D0*PC_uh*i_*svp*svm*F12*F11r*c1*f1*u1*u2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F12*F10r*c0*f1*u0*u2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F11*F15r*c5*f1*u1*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F11*F14r*c4*f1*u1*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F11*F13r*c3*f1*u1*u3*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F11*F12r*c2*f1*u1*u2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F11*F11r*c1*f1*u1**2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F11*F10r*c0*f1*u0*u1*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F10*F15r*c5*f1*u0*u5*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F10*F14r*c4*f1*u0*u4*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F10*F13r*c3*f1*u0*u3*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F10*F12r*c2*f1*u0*u2*w1*w2
     &  + 8.D0*PC_uh*i_*svp*svm*F10*F11r*c1*f1*u0*u1*w1*w2
     &
      traza1 = traza1 + 8.D0*PC_uh*i_*svp*svm*F10*F10r*c0*f1*u0**2*w1*
     & w2
     &  - 4.D0*PC_uh*i_*ssm*svp*F25*F15r*c5*f1*u5**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F25*F14r*c4*f1*u4*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F25*F13r*c3*f1*u3*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F25*F12r*c2*f1*u2*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F25*F11r*c1*f1*u1*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F25*F10r*c0*f1*u0*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F24*F15r*c5*f1*u4*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F24*F14r*c4*f1*u4**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F24*F13r*c3*f1*u3*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F24*F12r*c2*f1*u2*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F24*F11r*c1*f1*u1*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F24*F10r*c0*f1*u0*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F23*F15r*c5*f1*u3*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssm*svp*F23*F14r*c4*f1*u3*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F23*F13r*c3*f1*u3**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F23*F12r*c2*f1*u2*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F23*F11r*c1*f1*u1*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F23*F10r*c0*f1*u0*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F22*F15r*c5*f1*u2*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F22*F14r*c4*f1*u2*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F22*F13r*c3*f1*u2*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F22*F12r*c2*f1*u2**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F22*F11r*c1*f1*u1*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F22*F10r*c0*f1*u0*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F21*F15r*c5*f1*u1*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F21*F14r*c4*f1*u1*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F21*F13r*c3*f1*u1*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F21*F12r*c2*f1*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssm*svp*F21*F11r*c1*f1*u1**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F21*F10r*c0*f1*u0*u1*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F20*F15r*c5*f1*u0*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F20*F14r*c4*f1*u0*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F20*F13r*c3*f1*u0*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F20*F12r*c2*f1*u0*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F20*F11r*c1*f1*u0*u1*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F20*F10r*c0*f1*u0**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F15*F25r*c5*f2*u5**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F15*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F15*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F15*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F15*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F15*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F14*F25r*c5*f2*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssm*svp*F14*F24r*c4*f2*u4**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F14*F23r*c3*f2*u3*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F14*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F14*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F14*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F13*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F13*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F13*F23r*c3*f2*u3**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F13*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F13*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F13*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F12*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F12*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F12*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F12*F22r*c2*f2*u2**2*w1
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssm*svp*F12*F21r*c1*f2*u1*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F12*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F11*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F11*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F11*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F11*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F11*F21r*c1*f2*u1**2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F11*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F10*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F10*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F10*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F10*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F10*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*PC_uh*i_*ssm*svp*F10*F20r*c0*f2*u0**2*w1
     &  - 4.D0*PC_uh*i_*ssp*svm*F25*F15r*c5*f1*u5**2*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssp*svm*F25*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F25*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F25*F12r*c2*f1*u2*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F25*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F25*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F24*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F24*F14r*c4*f1*u4**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F24*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F24*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F24*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F24*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F23*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F23*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F23*F13r*c3*f1*u3**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F23*F12r*c2*f1*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssp*svm*F23*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F23*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F22*F15r*c5*f1*u2*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F22*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F22*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F22*F12r*c2*f1*u2**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F22*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F22*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F21*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F21*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F21*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F21*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F21*F11r*c1*f1*u1**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F21*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F20*F15r*c5*f1*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssp*svm*F20*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F20*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F20*F12r*c2*f1*u0*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F20*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F20*F10r*c0*f1*u0**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F15*F25r*c5*f2*u5**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F15*F24r*c4*f2*u4*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F15*F23r*c3*f2*u3*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F15*F22r*c2*f2*u2*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F15*F21r*c1*f2*u1*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F15*F20r*c0*f2*u0*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F14*F25r*c5*f2*u4*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F14*F24r*c4*f2*u4**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F14*F23r*c3*f2*u3*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F14*F22r*c2*f2*u2*u4*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssp*svm*F14*F21r*c1*f2*u1*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F14*F20r*c0*f2*u0*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F13*F25r*c5*f2*u3*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F13*F24r*c4*f2*u3*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F13*F23r*c3*f2*u3**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F13*F22r*c2*f2*u2*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F13*F21r*c1*f2*u1*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F13*F20r*c0*f2*u0*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F12*F25r*c5*f2*u2*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F12*F24r*c4*f2*u2*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F12*F23r*c3*f2*u2*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F12*F22r*c2*f2*u2**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F12*F21r*c1*f2*u1*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F12*F20r*c0*f2*u0*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F11*F25r*c5*f2*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*PC_uh*i_*ssp*svm*F11*F24r*c4*f2*u1*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F11*F23r*c3*f2*u1*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F11*F22r*c2*f2*u1*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F11*F21r*c1*f2*u1**2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F11*F20r*c0*f2*u0*u1*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F10*F25r*c5*f2*u0*u5*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F10*F24r*c4*f2*u0*u4*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F10*F23r*c3*f2*u0*u3*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F10*F22r*c2*f2*u0*u2*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F10*F21r*c1*f2*u0*u1*w2
     &  - 4.D0*PC_uh*i_*ssp*svm*F10*F20r*c0*f2*u0**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F15r*c5*f1*u5**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F15r*c5*f1*u5**2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F14r*c4*f1*u4*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F14r*c4*f1*u4*u5*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F45*F13r*c3*f1*u3*u5*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F13r*c3*f1*u3*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F12r*c2*f1*u2*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F12r*c2*f1*u2*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F11r*c1*f1*u1*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F11r*c1*f1*u1*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F10r*c0*f1*u0*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F45*F10r*c0*f1*u0*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F15r*c5*f1*u4*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F15r*c5*f1*u4*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F14r*c4*f1*u4**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F14r*c4*f1*u4**2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F13r*c3*f1*u3*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F13r*c3*f1*u3*u4*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F44*F12r*c2*f1*u2*u4*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F12r*c2*f1*u2*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F11r*c1*f1*u1*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F11r*c1*f1*u1*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F10r*c0*f1*u0*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F44*F10r*c0*f1*u0*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F15r*c5*f1*u3*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F15r*c5*f1*u3*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F14r*c4*f1*u3*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F14r*c4*f1*u3*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F13r*c3*f1*u3**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F13r*c3*f1*u3**2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F12r*c2*f1*u2*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F12r*c2*f1*u2*u3*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F43*F11r*c1*f1*u1*u3*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F11r*c1*f1*u1*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F10r*c0*f1*u0*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F43*F10r*c0*f1*u0*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F15r*c5*f1*u2*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F15r*c5*f1*u2*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F14r*c4*f1*u2*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F14r*c4*f1*u2*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F13r*c3*f1*u2*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F13r*c3*f1*u2*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F12r*c2*f1*u2**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F12r*c2*f1*u2**2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F11r*c1*f1*u1*u2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F11r*c1*f1*u1*u2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F42*F10r*c0*f1*u0*u2*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F42*F10r*c0*f1*u0*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F15r*c5*f1*u1*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F15r*c5*f1*u1*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F14r*c4*f1*u1*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F14r*c4*f1*u1*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F13r*c3*f1*u1*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F13r*c3*f1*u1*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F12r*c2*f1*u1*u2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F12r*c2*f1*u1*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F11r*c1*f1*u1**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F11r*c1*f1*u1**2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F10r*c0*f1*u0*u1*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F41*F10r*c0*f1*u0*u1*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F40*F15r*c5*f1*u0*u5*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F15r*c5*f1*u0*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F14r*c4*f1*u0*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F14r*c4*f1*u0*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F13r*c3*f1*u0*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F13r*c3*f1*u0*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F12r*c2*f1*u0*u2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F12r*c2*f1*u0*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F11r*c1*f1*u0*u1*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F11r*c1*f1*u0*u1*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F10r*c0*f1*u0**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F40*F10r*c0*f1*u0**2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F45r*c5*f4*u5**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F45r*c5*f4*u5**2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F15*F44r*c4*f4*u4*u5*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F44r*c4*f4*u4*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F43r*c3*f4*u3*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F42r*c2*f4*u2*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F41r*c1*f4*u1*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F15*F40r*c0*f4*u0*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F45r*c5*f4*u4*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F44r*c4*f4*u4**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F44r*c4*f4*u4**2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F14*F43r*c3*f4*u3*u4*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F43r*c3*f4*u3*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F42r*c2*f4*u2*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F41r*c1*f4*u1*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F14*F40r*c0*f4*u0*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F45r*c5*f4*u3*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F44r*c4*f4*u3*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F43r*c3*f4*u3**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F43r*c3*f4*u3**2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F13*F42r*c2*f4*u2*u3*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F42r*c2*f4*u2*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F41r*c1*f4*u1*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F13*F40r*c0*f4*u0*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F45r*c5*f4*u2*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F44r*c4*f4*u2*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F43r*c3*f4*u2*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F42r*c2*f4*u2**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F42r*c2*f4*u2**2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F12*F41r*c1*f4*u1*u2*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F41r*c1*f4*u1*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F12*F40r*c0*f4*u0*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F45r*c5*f4*u1*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F44r*c4*f4*u1*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F43r*c3*f4*u1*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F42r*c2*f4*u1*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F41r*c1*f4*u1**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F41r*c1*f4*u1**2*w1
     &
      traza1 = traza1 + 4.D0*PC_uh*i_*qq*svp*svm*F11*F40r*c0*f4*u0*u1*
     & w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F11*F40r*c0*f4*u0*u1*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F45r*c5*f4*u0*u5*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F44r*c4*f4*u0*u4*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F43r*c3*f4*u0*u3*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F42r*c2*f4*u0*u2*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F41r*c1*f4*u0*u1*w1
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F40r*c0*f4*u0**2*w2
     &  + 4.D0*PC_uh*i_*qq*svp*svm*F10*F40r*c0*f4*u0**2*w1
     &
      traza1 = traza1 - 8.D0*PC_uh*i_*PC2*svp*svm*F25*F25r*c5*f2*u5**2*
     & w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F25*F24r*c4*f2*u4*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F25*F23r*c3*f2*u3*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F25*F22r*c2*f2*u2*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F25*F21r*c1*f2*u1*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F25*F20r*c0*f2*u0*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F24*F25r*c5*f2*u4*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F24*F24r*c4*f2*u4**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F24*F23r*c3*f2*u3*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F24*F22r*c2*f2*u2*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F24*F21r*c1*f2*u1*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F24*F20r*c0*f2*u0*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F23*F25r*c5*f2*u3*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F23*F24r*c4*f2*u3*u4*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*i_*PC2*svp*svm*F23*F23r*c3*f2*u3**2*
     & w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F23*F22r*c2*f2*u2*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F23*F21r*c1*f2*u1*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F23*F20r*c0*f2*u0*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F22*F25r*c5*f2*u2*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F22*F24r*c4*f2*u2*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F22*F23r*c3*f2*u2*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F22*F22r*c2*f2*u2**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F22*F21r*c1*f2*u1*u2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F22*F20r*c0*f2*u0*u2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F21*F25r*c5*f2*u1*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F21*F24r*c4*f2*u1*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F21*F23r*c3*f2*u1*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F21*F22r*c2*f2*u1*u2*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*i_*PC2*svp*svm*F21*F21r*c1*f2*u1**2*
     & w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F21*F20r*c0*f2*u0*u1*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F20*F25r*c5*f2*u0*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F20*F24r*c4*f2*u0*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F20*F23r*c3*f2*u0*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F20*F22r*c2*f2*u0*u2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F20*F21r*c1*f2*u0*u1*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*svp*svm*F20*F20r*c0*f2*u0**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F45*F45r*c5*f4*u5**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F45*F44r*c4*f4*u4*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F45*F43r*c3*f4*u3*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F45*F42r*c2*f4*u2*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F45*F41r*c1*f4*u1*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F45*F40r*c0*f4*u0*u5*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F44*F45r*c5*f4*u4*
     & u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F44*F44r*c4*f4*u4**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F44*F43r*c3*f4*u3*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F44*F42r*c2*f4*u2*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F44*F41r*c1*f4*u1*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F44*F40r*c0*f4*u0*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F43*F45r*c5*f4*u3*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F43*F44r*c4*f4*u3*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F43*F43r*c3*f4*u3**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F43*F42r*c2*f4*u2*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F43*F41r*c1*f4*u1*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F43*F40r*c0*f4*u0*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F42*F45r*c5*f4*u2*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F42*F44r*c4*f4*u2*u4*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F42*F43r*c3*f4*u2*
     & u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F42*F42r*c2*f4*u2**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F42*F41r*c1*f4*u1*u2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F42*F40r*c0*f4*u0*u2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F41*F45r*c5*f4*u1*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F41*F44r*c4*f4*u1*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F41*F43r*c3*f4*u1*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F41*F42r*c2*f4*u1*u2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F41*F41r*c1*f4*u1**2*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F41*F40r*c0*f4*u0*u1*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F40*F45r*c5*f4*u0*u5*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F40*F44r*c4*f4*u0*u4*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F40*F43r*c3*f4*u0*u3*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F40*F42r*c2*f4*u0*u2*w1*w2
     &
      traza1 = traza1 - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F40*F41r*c1*f4*u0*
     & u1*w1*w2
     &  - 8.D0*PC_uh*i_*PC2*qq*svp*svm*F40*F40r*c0*f4*u0**2*w1*w2
     &  - 4.D0*q_uh*svp*svm*F15*F15i*c5*f1*u5**2*w2
     &  + 4.D0*q_uh*svp*svm*F15*F15i*c5*f1*u5**2*w1
     &  - 4.D0*q_uh*svp*svm*F15*F14i*c4*f1*u4*u5*w2
     &  + 4.D0*q_uh*svp*svm*F15*F14i*c4*f1*u4*u5*w1
     &  - 4.D0*q_uh*svp*svm*F15*F13i*c3*f1*u3*u5*w2
     &  + 4.D0*q_uh*svp*svm*F15*F13i*c3*f1*u3*u5*w1
     &  - 4.D0*q_uh*svp*svm*F15*F12i*c2*f1*u2*u5*w2
     &  + 4.D0*q_uh*svp*svm*F15*F12i*c2*f1*u2*u5*w1
     &  - 4.D0*q_uh*svp*svm*F15*F11i*c1*f1*u1*u5*w2
     &  + 4.D0*q_uh*svp*svm*F15*F11i*c1*f1*u1*u5*w1
     &  - 4.D0*q_uh*svp*svm*F15*F10i*c0*f1*u0*u5*w2
     &  + 4.D0*q_uh*svp*svm*F15*F10i*c0*f1*u0*u5*w1
     &
      traza1 = traza1 - 4.D0*q_uh*svp*svm*F14*F15i*c5*f1*u4*u5*w2
     &  + 4.D0*q_uh*svp*svm*F14*F15i*c5*f1*u4*u5*w1
     &  - 4.D0*q_uh*svp*svm*F14*F14i*c4*f1*u4**2*w2
     &  + 4.D0*q_uh*svp*svm*F14*F14i*c4*f1*u4**2*w1
     &  - 4.D0*q_uh*svp*svm*F14*F13i*c3*f1*u3*u4*w2
     &  + 4.D0*q_uh*svp*svm*F14*F13i*c3*f1*u3*u4*w1
     &  - 4.D0*q_uh*svp*svm*F14*F12i*c2*f1*u2*u4*w2
     &  + 4.D0*q_uh*svp*svm*F14*F12i*c2*f1*u2*u4*w1
     &  - 4.D0*q_uh*svp*svm*F14*F11i*c1*f1*u1*u4*w2
     &  + 4.D0*q_uh*svp*svm*F14*F11i*c1*f1*u1*u4*w1
     &  - 4.D0*q_uh*svp*svm*F14*F10i*c0*f1*u0*u4*w2
     &  + 4.D0*q_uh*svp*svm*F14*F10i*c0*f1*u0*u4*w1
     &  - 4.D0*q_uh*svp*svm*F13*F15i*c5*f1*u3*u5*w2
     &  + 4.D0*q_uh*svp*svm*F13*F15i*c5*f1*u3*u5*w1
     &  - 4.D0*q_uh*svp*svm*F13*F14i*c4*f1*u3*u4*w2
     &
      traza1 = traza1 + 4.D0*q_uh*svp*svm*F13*F14i*c4*f1*u3*u4*w1
     &  - 4.D0*q_uh*svp*svm*F13*F13i*c3*f1*u3**2*w2
     &  + 4.D0*q_uh*svp*svm*F13*F13i*c3*f1*u3**2*w1
     &  - 4.D0*q_uh*svp*svm*F13*F12i*c2*f1*u2*u3*w2
     &  + 4.D0*q_uh*svp*svm*F13*F12i*c2*f1*u2*u3*w1
     &  - 4.D0*q_uh*svp*svm*F13*F11i*c1*f1*u1*u3*w2
     &  + 4.D0*q_uh*svp*svm*F13*F11i*c1*f1*u1*u3*w1
     &  - 4.D0*q_uh*svp*svm*F13*F10i*c0*f1*u0*u3*w2
     &  + 4.D0*q_uh*svp*svm*F13*F10i*c0*f1*u0*u3*w1
     &  - 4.D0*q_uh*svp*svm*F12*F15i*c5*f1*u2*u5*w2
     &  + 4.D0*q_uh*svp*svm*F12*F15i*c5*f1*u2*u5*w1
     &  - 4.D0*q_uh*svp*svm*F12*F14i*c4*f1*u2*u4*w2
     &  + 4.D0*q_uh*svp*svm*F12*F14i*c4*f1*u2*u4*w1
     &  - 4.D0*q_uh*svp*svm*F12*F13i*c3*f1*u2*u3*w2
     &  + 4.D0*q_uh*svp*svm*F12*F13i*c3*f1*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*q_uh*svp*svm*F12*F12i*c2*f1*u2**2*w2
     &  + 4.D0*q_uh*svp*svm*F12*F12i*c2*f1*u2**2*w1
     &  - 4.D0*q_uh*svp*svm*F12*F11i*c1*f1*u1*u2*w2
     &  + 4.D0*q_uh*svp*svm*F12*F11i*c1*f1*u1*u2*w1
     &  - 4.D0*q_uh*svp*svm*F12*F10i*c0*f1*u0*u2*w2
     &  + 4.D0*q_uh*svp*svm*F12*F10i*c0*f1*u0*u2*w1
     &  - 4.D0*q_uh*svp*svm*F11*F15i*c5*f1*u1*u5*w2
     &  + 4.D0*q_uh*svp*svm*F11*F15i*c5*f1*u1*u5*w1
     &  - 4.D0*q_uh*svp*svm*F11*F14i*c4*f1*u1*u4*w2
     &  + 4.D0*q_uh*svp*svm*F11*F14i*c4*f1*u1*u4*w1
     &  - 4.D0*q_uh*svp*svm*F11*F13i*c3*f1*u1*u3*w2
     &  + 4.D0*q_uh*svp*svm*F11*F13i*c3*f1*u1*u3*w1
     &  - 4.D0*q_uh*svp*svm*F11*F12i*c2*f1*u1*u2*w2
     &  + 4.D0*q_uh*svp*svm*F11*F12i*c2*f1*u1*u2*w1
     &  - 4.D0*q_uh*svp*svm*F11*F11i*c1*f1*u1**2*w2
     &
      traza1 = traza1 + 4.D0*q_uh*svp*svm*F11*F11i*c1*f1*u1**2*w1
     &  - 4.D0*q_uh*svp*svm*F11*F10i*c0*f1*u0*u1*w2
     &  + 4.D0*q_uh*svp*svm*F11*F10i*c0*f1*u0*u1*w1
     &  - 4.D0*q_uh*svp*svm*F10*F15i*c5*f1*u0*u5*w2
     &  + 4.D0*q_uh*svp*svm*F10*F15i*c5*f1*u0*u5*w1
     &  - 4.D0*q_uh*svp*svm*F10*F14i*c4*f1*u0*u4*w2
     &  + 4.D0*q_uh*svp*svm*F10*F14i*c4*f1*u0*u4*w1
     &  - 4.D0*q_uh*svp*svm*F10*F13i*c3*f1*u0*u3*w2
     &  + 4.D0*q_uh*svp*svm*F10*F13i*c3*f1*u0*u3*w1
     &  - 4.D0*q_uh*svp*svm*F10*F12i*c2*f1*u0*u2*w2
     &  + 4.D0*q_uh*svp*svm*F10*F12i*c2*f1*u0*u2*w1
     &  - 4.D0*q_uh*svp*svm*F10*F11i*c1*f1*u0*u1*w2
     &  + 4.D0*q_uh*svp*svm*F10*F11i*c1*f1*u0*u1*w1
     &  - 4.D0*q_uh*svp*svm*F10*F10i*c0*f1*u0**2*w2
     &  + 4.D0*q_uh*svp*svm*F10*F10i*c0*f1*u0**2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*svp*svm*F25*F25i*c5*f2*u5**2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F25*F25i*c5*f2*u5**2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F25*F24i*c4*f2*u4*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F25*F24i*c4*f2*u4*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F25*F23i*c3*f2*u3*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F25*F23i*c3*f2*u3*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F25*F22i*c2*f2*u2*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F25*F22i*c2*f2*u2*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F25*F21i*c1*f2*u1*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F25*F21i*c1*f2*u1*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F25*F20i*c0*f2*u0*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F25*F20i*c0*f2*u0*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F24*F25i*c5*f2*u4*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F24*F25i*c5*f2*u4*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F24*F24i*c4*f2*u4**2*w2
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*svp*svm*F24*F24i*c4*f2*u4**2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F24*F23i*c3*f2*u3*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F24*F23i*c3*f2*u3*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F24*F22i*c2*f2*u2*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F24*F22i*c2*f2*u2*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F24*F21i*c1*f2*u1*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F24*F21i*c1*f2*u1*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F24*F20i*c0*f2*u0*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F24*F20i*c0*f2*u0*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F23*F25i*c5*f2*u3*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F23*F25i*c5*f2*u3*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F23*F24i*c4*f2*u3*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F23*F24i*c4*f2*u3*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F23*F23i*c3*f2*u3**2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F23*F23i*c3*f2*u3**2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*svp*svm*F23*F22i*c2*f2*u2*u3*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F23*F22i*c2*f2*u2*u3*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F23*F21i*c1*f2*u1*u3*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F23*F21i*c1*f2*u1*u3*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F23*F20i*c0*f2*u0*u3*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F23*F20i*c0*f2*u0*u3*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F22*F25i*c5*f2*u2*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F22*F25i*c5*f2*u2*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F22*F24i*c4*f2*u2*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F22*F24i*c4*f2*u2*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F22*F23i*c3*f2*u2*u3*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F22*F23i*c3*f2*u2*u3*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F22*F22i*c2*f2*u2**2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F22*F22i*c2*f2*u2**2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F22*F21i*c1*f2*u1*u2*w2
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*svp*svm*F22*F21i*c1*f2*u1*u2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F22*F20i*c0*f2*u0*u2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F22*F20i*c0*f2*u0*u2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F21*F25i*c5*f2*u1*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F21*F25i*c5*f2*u1*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F21*F24i*c4*f2*u1*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F21*F24i*c4*f2*u1*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F21*F23i*c3*f2*u1*u3*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F21*F23i*c3*f2*u1*u3*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F21*F22i*c2*f2*u1*u2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F21*F22i*c2*f2*u1*u2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F21*F21i*c1*f2*u1**2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F21*F21i*c1*f2*u1**2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F21*F20i*c0*f2*u0*u1*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F21*F20i*c0*f2*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*svp*svm*F20*F25i*c5*f2*u0*u5*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F20*F25i*c5*f2*u0*u5*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F20*F24i*c4*f2*u0*u4*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F20*F24i*c4*f2*u0*u4*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F20*F23i*c3*f2*u0*u3*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F20*F23i*c3*f2*u0*u3*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F20*F22i*c2*f2*u0*u2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F20*F22i*c2*f2*u0*u2*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F20*F21i*c1*f2*u0*u1*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F20*F21i*c1*f2*u0*u1*w1
     &  - 4.D0*q_uh*PC2*svp*svm*F20*F20i*c0*f2*u0**2*w2
     &  + 4.D0*q_uh*PC2*svp*svm*F20*F20i*c0*f2*u0**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F45*F25i*c5*f2*u5**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F45*F24i*c4*f2*u4*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F45*F23i*c3*f2*u3*u5*w1
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*ssm*svp*F45*F22i*c2*f2*u2*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F45*F21i*c1*f2*u1*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F45*F20i*c0*f2*u0*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F44*F25i*c5*f2*u4*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F44*F24i*c4*f2*u4**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F44*F23i*c3*f2*u3*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F44*F22i*c2*f2*u2*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F44*F21i*c1*f2*u1*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F44*F20i*c0*f2*u0*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F43*F25i*c5*f2*u3*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F43*F24i*c4*f2*u3*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F43*F23i*c3*f2*u3**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F43*F22i*c2*f2*u2*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F43*F21i*c1*f2*u1*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F43*F20i*c0*f2*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*ssm*svp*F42*F25i*c5*f2*u2*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F42*F24i*c4*f2*u2*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F42*F23i*c3*f2*u2*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F42*F22i*c2*f2*u2**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F42*F21i*c1*f2*u1*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F42*F20i*c0*f2*u0*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F41*F25i*c5*f2*u1*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F41*F24i*c4*f2*u1*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F41*F23i*c3*f2*u1*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F41*F22i*c2*f2*u1*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F41*F21i*c1*f2*u1**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F41*F20i*c0*f2*u0*u1*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F40*F25i*c5*f2*u0*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F40*F24i*c4*f2*u0*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F40*F23i*c3*f2*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*ssm*svp*F40*F22i*c2*f2*u0*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F40*F21i*c1*f2*u0*u1*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F40*F20i*c0*f2*u0**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F25*F45i*c5*f4*u5**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F25*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F25*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F25*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F25*F41i*c1*f4*u1*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F25*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F24*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F24*F44i*c4*f4*u4**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F24*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F24*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F24*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F24*F40i*c0*f4*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*ssm*svp*F23*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F23*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F23*F43i*c3*f4*u3**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F23*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F23*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F23*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F22*F45i*c5*f4*u2*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F22*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F22*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F22*F42i*c2*f4*u2**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F22*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F22*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F21*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F21*F44i*c4*f4*u1*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F21*F43i*c3*f4*u1*u3*w1
     &
      traza1 = traza1 + 4.D0*q_uh*PC2*ssm*svp*F21*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F21*F41i*c1*f4*u1**2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F21*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F20*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F20*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F20*F43i*c3*f4*u0*u3*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F20*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F20*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*q_uh*PC2*ssm*svp*F20*F40i*c0*f4*u0**2*w1
     &  - 4.D0*q_uh*PC2*ssp*svm*F45*F25i*c5*f2*u5**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F45*F24i*c4*f2*u4*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F45*F23i*c3*f2*u3*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F45*F22i*c2*f2*u2*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F45*F21i*c1*f2*u1*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F45*F20i*c0*f2*u0*u5*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*ssp*svm*F44*F25i*c5*f2*u4*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F44*F24i*c4*f2*u4**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F44*F23i*c3*f2*u3*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F44*F22i*c2*f2*u2*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F44*F21i*c1*f2*u1*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F44*F20i*c0*f2*u0*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F43*F25i*c5*f2*u3*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F43*F24i*c4*f2*u3*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F43*F23i*c3*f2*u3**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F43*F22i*c2*f2*u2*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F43*F21i*c1*f2*u1*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F43*F20i*c0*f2*u0*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F42*F25i*c5*f2*u2*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F42*F24i*c4*f2*u2*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F42*F23i*c3*f2*u2*u3*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*ssp*svm*F42*F22i*c2*f2*u2**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F42*F21i*c1*f2*u1*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F42*F20i*c0*f2*u0*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F41*F25i*c5*f2*u1*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F41*F24i*c4*f2*u1*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F41*F23i*c3*f2*u1*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F41*F22i*c2*f2*u1*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F41*F21i*c1*f2*u1**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F41*F20i*c0*f2*u0*u1*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F40*F25i*c5*f2*u0*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F40*F24i*c4*f2*u0*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F40*F23i*c3*f2*u0*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F40*F22i*c2*f2*u0*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F40*F21i*c1*f2*u0*u1*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F40*F20i*c0*f2*u0**2*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*ssp*svm*F25*F45i*c5*f4*u5**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F25*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F25*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F25*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F25*F41i*c1*f4*u1*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F25*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F24*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F24*F44i*c4*f4*u4**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F24*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F24*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F24*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F24*F40i*c0*f4*u0*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F23*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F23*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F23*F43i*c3*f4*u3**2*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*ssp*svm*F23*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F23*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F23*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F22*F45i*c5*f4*u2*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F22*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F22*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F22*F42i*c2*f4*u2**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F22*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F22*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F21*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F21*F44i*c4*f4*u1*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F21*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F21*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F21*F41i*c1*f4*u1**2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F21*F40i*c0*f4*u0*u1*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*ssp*svm*F20*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F20*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F20*F43i*c3*f4*u0*u3*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F20*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F20*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*q_uh*PC2*ssp*svm*F20*F40i*c0*f4*u0**2*w2
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F45*F45i*c5*f4*u5**2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F45*F45i*c5*f4*u5**2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F45*F44i*c4*f4*u4*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F45*F44i*c4*f4*u4*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F45*F43i*c3*f4*u3*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F45*F43i*c3*f4*u3*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F45*F42i*c2*f4*u2*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F45*F42i*c2*f4*u2*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F45*F41i*c1*f4*u1*u5*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*qq*svp*svm*F45*F41i*c1*f4*u1*u5*
     & w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F45*F40i*c0*f4*u0*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F45*F40i*c0*f4*u0*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F44*F45i*c5*f4*u4*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F44*F45i*c5*f4*u4*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F44*F44i*c4*f4*u4**2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F44*F44i*c4*f4*u4**2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F44*F43i*c3*f4*u3*u4*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F44*F43i*c3*f4*u3*u4*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F44*F42i*c2*f4*u2*u4*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F44*F42i*c2*f4*u2*u4*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F44*F41i*c1*f4*u1*u4*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F44*F41i*c1*f4*u1*u4*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F44*F40i*c0*f4*u0*u4*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*qq*svp*svm*F44*F40i*c0*f4*u0*u4*
     & w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F43*F45i*c5*f4*u3*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F43*F45i*c5*f4*u3*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F43*F44i*c4*f4*u3*u4*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F43*F44i*c4*f4*u3*u4*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F43*F43i*c3*f4*u3**2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F43*F43i*c3*f4*u3**2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F43*F42i*c2*f4*u2*u3*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F43*F42i*c2*f4*u2*u3*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F43*F41i*c1*f4*u1*u3*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F43*F41i*c1*f4*u1*u3*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F43*F40i*c0*f4*u0*u3*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F43*F40i*c0*f4*u0*u3*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F42*F45i*c5*f4*u2*u5*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*qq*svp*svm*F42*F45i*c5*f4*u2*u5*
     & w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F42*F44i*c4*f4*u2*u4*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F42*F44i*c4*f4*u2*u4*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F42*F43i*c3*f4*u2*u3*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F42*F43i*c3*f4*u2*u3*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F42*F42i*c2*f4*u2**2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F42*F42i*c2*f4*u2**2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F42*F41i*c1*f4*u1*u2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F42*F41i*c1*f4*u1*u2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F42*F40i*c0*f4*u0*u2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F42*F40i*c0*f4*u0*u2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F41*F45i*c5*f4*u1*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F41*F45i*c5*f4*u1*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F41*F44i*c4*f4*u1*u4*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*qq*svp*svm*F41*F44i*c4*f4*u1*u4*
     & w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F41*F43i*c3*f4*u1*u3*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F41*F43i*c3*f4*u1*u3*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F41*F42i*c2*f4*u1*u2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F41*F42i*c2*f4*u1*u2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F41*F41i*c1*f4*u1**2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F41*F41i*c1*f4*u1**2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F41*F40i*c0*f4*u0*u1*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F41*F40i*c0*f4*u0*u1*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F40*F45i*c5*f4*u0*u5*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F40*F45i*c5*f4*u0*u5*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F40*F44i*c4*f4*u0*u4*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F40*F44i*c4*f4*u0*u4*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F40*F43i*c3*f4*u0*u3*w2
     &
      traza1 = traza1 - 4.D0*q_uh*PC2*qq*svp*svm*F40*F43i*c3*f4*u0*u3*
     & w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F40*F42i*c2*f4*u0*u2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F40*F42i*c2*f4*u0*u2*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F40*F41i*c1*f4*u0*u1*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F40*F41i*c1*f4*u0*u1*w1
     &  + 4.D0*q_uh*PC2*qq*svp*svm*F40*F40i*c0*f4*u0**2*w2
     &  - 4.D0*q_uh*PC2*qq*svp*svm*F40*F40i*c0*f4*u0**2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F15*F15r*c5*f1*u5**2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F15*F15r*c5*f1*u5**2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F15*F14r*c4*f1*u4*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F15*F14r*c4*f1*u4*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F15*F13r*c3*f1*u3*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F15*F13r*c3*f1*u3*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F15*F12r*c2*f1*u2*u5*w2
     &
      traza1 = traza1 - 4.D0*q_uh*i_*svp*svm*F15*F12r*c2*f1*u2*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F15*F11r*c1*f1*u1*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F15*F11r*c1*f1*u1*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F15*F10r*c0*f1*u0*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F15*F10r*c0*f1*u0*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F14*F15r*c5*f1*u4*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F14*F15r*c5*f1*u4*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F14*F14r*c4*f1*u4**2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F14*F14r*c4*f1*u4**2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F14*F13r*c3*f1*u3*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F14*F13r*c3*f1*u3*u4*w1
     &  + 4.D0*q_uh*i_*svp*svm*F14*F12r*c2*f1*u2*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F14*F12r*c2*f1*u2*u4*w1
     &  + 4.D0*q_uh*i_*svp*svm*F14*F11r*c1*f1*u1*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F14*F11r*c1*f1*u1*u4*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*svp*svm*F14*F10r*c0*f1*u0*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F14*F10r*c0*f1*u0*u4*w1
     &  + 4.D0*q_uh*i_*svp*svm*F13*F15r*c5*f1*u3*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F13*F15r*c5*f1*u3*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F13*F14r*c4*f1*u3*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F13*F14r*c4*f1*u3*u4*w1
     &  + 4.D0*q_uh*i_*svp*svm*F13*F13r*c3*f1*u3**2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F13*F13r*c3*f1*u3**2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F13*F12r*c2*f1*u2*u3*w2
     &  - 4.D0*q_uh*i_*svp*svm*F13*F12r*c2*f1*u2*u3*w1
     &  + 4.D0*q_uh*i_*svp*svm*F13*F11r*c1*f1*u1*u3*w2
     &  - 4.D0*q_uh*i_*svp*svm*F13*F11r*c1*f1*u1*u3*w1
     &  + 4.D0*q_uh*i_*svp*svm*F13*F10r*c0*f1*u0*u3*w2
     &  - 4.D0*q_uh*i_*svp*svm*F13*F10r*c0*f1*u0*u3*w1
     &  + 4.D0*q_uh*i_*svp*svm*F12*F15r*c5*f1*u2*u5*w2
     &
      traza1 = traza1 - 4.D0*q_uh*i_*svp*svm*F12*F15r*c5*f1*u2*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F12*F14r*c4*f1*u2*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F12*F14r*c4*f1*u2*u4*w1
     &  + 4.D0*q_uh*i_*svp*svm*F12*F13r*c3*f1*u2*u3*w2
     &  - 4.D0*q_uh*i_*svp*svm*F12*F13r*c3*f1*u2*u3*w1
     &  + 4.D0*q_uh*i_*svp*svm*F12*F12r*c2*f1*u2**2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F12*F12r*c2*f1*u2**2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F12*F11r*c1*f1*u1*u2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F12*F11r*c1*f1*u1*u2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F12*F10r*c0*f1*u0*u2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F12*F10r*c0*f1*u0*u2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F11*F15r*c5*f1*u1*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F11*F15r*c5*f1*u1*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F11*F14r*c4*f1*u1*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F11*F14r*c4*f1*u1*u4*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*svp*svm*F11*F13r*c3*f1*u1*u3*w2
     &  - 4.D0*q_uh*i_*svp*svm*F11*F13r*c3*f1*u1*u3*w1
     &  + 4.D0*q_uh*i_*svp*svm*F11*F12r*c2*f1*u1*u2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F11*F12r*c2*f1*u1*u2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F11*F11r*c1*f1*u1**2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F11*F11r*c1*f1*u1**2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F11*F10r*c0*f1*u0*u1*w2
     &  - 4.D0*q_uh*i_*svp*svm*F11*F10r*c0*f1*u0*u1*w1
     &  + 4.D0*q_uh*i_*svp*svm*F10*F15r*c5*f1*u0*u5*w2
     &  - 4.D0*q_uh*i_*svp*svm*F10*F15r*c5*f1*u0*u5*w1
     &  + 4.D0*q_uh*i_*svp*svm*F10*F14r*c4*f1*u0*u4*w2
     &  - 4.D0*q_uh*i_*svp*svm*F10*F14r*c4*f1*u0*u4*w1
     &  + 4.D0*q_uh*i_*svp*svm*F10*F13r*c3*f1*u0*u3*w2
     &  - 4.D0*q_uh*i_*svp*svm*F10*F13r*c3*f1*u0*u3*w1
     &  + 4.D0*q_uh*i_*svp*svm*F10*F12r*c2*f1*u0*u2*w2
     &
      traza1 = traza1 - 4.D0*q_uh*i_*svp*svm*F10*F12r*c2*f1*u0*u2*w1
     &  + 4.D0*q_uh*i_*svp*svm*F10*F11r*c1*f1*u0*u1*w2
     &  - 4.D0*q_uh*i_*svp*svm*F10*F11r*c1*f1*u0*u1*w1
     &  + 4.D0*q_uh*i_*svp*svm*F10*F10r*c0*f1*u0**2*w2
     &  - 4.D0*q_uh*i_*svp*svm*F10*F10r*c0*f1*u0**2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F25*F25r*c5*f2*u5**2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F25*F25r*c5*f2*u5**2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F25*F24r*c4*f2*u4*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F25*F24r*c4*f2*u4*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F25*F23r*c3*f2*u3*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F25*F23r*c3*f2*u3*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F25*F22r*c2*f2*u2*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F25*F22r*c2*f2*u2*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F25*F21r*c1*f2*u1*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F25*F21r*c1*f2*u1*u5*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*svp*svm*F25*F20r*c0*f2*u0*u5*
     & w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F25*F20r*c0*f2*u0*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F24*F25r*c5*f2*u4*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F24*F25r*c5*f2*u4*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F24*F24r*c4*f2*u4**2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F24*F24r*c4*f2*u4**2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F24*F23r*c3*f2*u3*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F24*F23r*c3*f2*u3*u4*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F24*F22r*c2*f2*u2*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F24*F22r*c2*f2*u2*u4*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F24*F21r*c1*f2*u1*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F24*F21r*c1*f2*u1*u4*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F24*F20r*c0*f2*u0*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F24*F20r*c0*f2*u0*u4*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*svp*svm*F23*F25r*c5*f2*u3*u5*
     & w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F23*F25r*c5*f2*u3*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F23*F24r*c4*f2*u3*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F23*F24r*c4*f2*u3*u4*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F23*F23r*c3*f2*u3**2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F23*F23r*c3*f2*u3**2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F23*F22r*c2*f2*u2*u3*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F23*F22r*c2*f2*u2*u3*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F23*F21r*c1*f2*u1*u3*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F23*F21r*c1*f2*u1*u3*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F23*F20r*c0*f2*u0*u3*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F23*F20r*c0*f2*u0*u3*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F22*F25r*c5*f2*u2*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F22*F25r*c5*f2*u2*u5*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*svp*svm*F22*F24r*c4*f2*u2*u4*
     & w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F22*F24r*c4*f2*u2*u4*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F22*F23r*c3*f2*u2*u3*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F22*F23r*c3*f2*u2*u3*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F22*F22r*c2*f2*u2**2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F22*F22r*c2*f2*u2**2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F22*F21r*c1*f2*u1*u2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F22*F21r*c1*f2*u1*u2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F22*F20r*c0*f2*u0*u2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F22*F20r*c0*f2*u0*u2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F21*F25r*c5*f2*u1*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F21*F25r*c5*f2*u1*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F21*F24r*c4*f2*u1*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F21*F24r*c4*f2*u1*u4*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*svp*svm*F21*F23r*c3*f2*u1*u3*
     & w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F21*F23r*c3*f2*u1*u3*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F21*F22r*c2*f2*u1*u2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F21*F22r*c2*f2*u1*u2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F21*F21r*c1*f2*u1**2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F21*F21r*c1*f2*u1**2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F21*F20r*c0*f2*u0*u1*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F21*F20r*c0*f2*u0*u1*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F20*F25r*c5*f2*u0*u5*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F20*F25r*c5*f2*u0*u5*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F20*F24r*c4*f2*u0*u4*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F20*F24r*c4*f2*u0*u4*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F20*F23r*c3*f2*u0*u3*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F20*F23r*c3*f2*u0*u3*w1
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*svp*svm*F20*F22r*c2*f2*u0*u2*
     & w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F20*F22r*c2*f2*u0*u2*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F20*F21r*c1*f2*u0*u1*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F20*F21r*c1*f2*u0*u1*w1
     &  + 4.D0*q_uh*i_*PC2*svp*svm*F20*F20r*c0*f2*u0**2*w2
     &  - 4.D0*q_uh*i_*PC2*svp*svm*F20*F20r*c0*f2*u0**2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F45*F25r*c5*f2*u5**2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F45*F24r*c4*f2*u4*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F45*F23r*c3*f2*u3*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F45*F22r*c2*f2*u2*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F45*F21r*c1*f2*u1*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F45*F20r*c0*f2*u0*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F44*F25r*c5*f2*u4*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F44*F24r*c4*f2*u4**2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*ssm*svp*F44*F23r*c3*f2*u3*u4*
     & w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F44*F22r*c2*f2*u2*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F44*F21r*c1*f2*u1*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F44*F20r*c0*f2*u0*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F43*F25r*c5*f2*u3*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F43*F24r*c4*f2*u3*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F43*F23r*c3*f2*u3**2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F43*F22r*c2*f2*u2*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F43*F21r*c1*f2*u1*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F43*F20r*c0*f2*u0*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F42*F25r*c5*f2*u2*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F42*F24r*c4*f2*u2*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F42*F23r*c3*f2*u2*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F42*F22r*c2*f2*u2**2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*ssm*svp*F42*F21r*c1*f2*u1*u2*
     & w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F42*F20r*c0*f2*u0*u2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F41*F25r*c5*f2*u1*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F41*F24r*c4*f2*u1*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F41*F23r*c3*f2*u1*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F41*F22r*c2*f2*u1*u2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F41*F21r*c1*f2*u1**2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F41*F20r*c0*f2*u0*u1*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F40*F25r*c5*f2*u0*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F40*F24r*c4*f2*u0*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F40*F23r*c3*f2*u0*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F40*F22r*c2*f2*u0*u2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F40*F21r*c1*f2*u0*u1*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F40*F20r*c0*f2*u0**2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*ssm*svp*F25*F45r*c5*f4*u5**2*
     & w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F25*F44r*c4*f4*u4*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F25*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F25*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F25*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F25*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F24*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F24*F44r*c4*f4*u4**2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F24*F43r*c3*f4*u3*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F24*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F24*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F24*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F23*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F23*F44r*c4*f4*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*ssm*svp*F23*F43r*c3*f4*u3**2*
     & w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F23*F42r*c2*f4*u2*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F23*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F23*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F22*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F22*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F22*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F22*F42r*c2*f4*u2**2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F22*F41r*c1*f4*u1*u2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F22*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F21*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F21*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F21*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F21*F42r*c2*f4*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*ssm*svp*F21*F41r*c1*f4*u1**2*
     & w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F21*F40r*c0*f4*u0*u1*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F20*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F20*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F20*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F20*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F20*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*q_uh*i_*PC2*ssm*svp*F20*F40r*c0*f4*u0**2*w1
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F45*F25r*c5*f2*u5**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F45*F24r*c4*f2*u4*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F45*F23r*c3*f2*u3*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F45*F22r*c2*f2*u2*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F45*F21r*c1*f2*u1*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F45*F20r*c0*f2*u0*u5*w2
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*ssp*svm*F44*F25r*c5*f2*u4*u5*
     & w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F44*F24r*c4*f2*u4**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F44*F23r*c3*f2*u3*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F44*F22r*c2*f2*u2*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F44*F21r*c1*f2*u1*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F44*F20r*c0*f2*u0*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F43*F25r*c5*f2*u3*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F43*F24r*c4*f2*u3*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F43*F23r*c3*f2*u3**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F43*F22r*c2*f2*u2*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F43*F21r*c1*f2*u1*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F43*F20r*c0*f2*u0*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F42*F25r*c5*f2*u2*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F42*F24r*c4*f2*u2*u4*w2
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*ssp*svm*F42*F23r*c3*f2*u2*u3*
     & w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F42*F22r*c2*f2*u2**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F42*F21r*c1*f2*u1*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F42*F20r*c0*f2*u0*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F41*F25r*c5*f2*u1*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F41*F24r*c4*f2*u1*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F41*F23r*c3*f2*u1*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F41*F22r*c2*f2*u1*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F41*F21r*c1*f2*u1**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F41*F20r*c0*f2*u0*u1*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F40*F25r*c5*f2*u0*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F40*F24r*c4*f2*u0*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F40*F23r*c3*f2*u0*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F40*F22r*c2*f2*u0*u2*w2
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*ssp*svm*F40*F21r*c1*f2*u0*u1*
     & w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F40*F20r*c0*f2*u0**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F25*F45r*c5*f4*u5**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F25*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F25*F43r*c3*f4*u3*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F25*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F25*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F25*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F24*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F24*F44r*c4*f4*u4**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F24*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F24*F42r*c2*f4*u2*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F24*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F24*F40r*c0*f4*u0*u4*w2
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*ssp*svm*F23*F45r*c5*f4*u3*u5*
     & w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F23*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F23*F43r*c3*f4*u3**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F23*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F23*F41r*c1*f4*u1*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F23*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F22*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F22*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F22*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F22*F42r*c2*f4*u2**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F22*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F22*F40r*c0*f4*u0*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F21*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F21*F44r*c4*f4*u1*u4*w2
     &
      traza1 = traza1 + 4.D0*q_uh*i_*PC2*ssp*svm*F21*F43r*c3*f4*u1*u3*
     & w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F21*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F21*F41r*c1*f4*u1**2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F21*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F20*F45r*c5*f4*u0*u5*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F20*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F20*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F20*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F20*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*q_uh*i_*PC2*ssp*svm*F20*F40r*c0*f4*u0**2*w2
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F45r*c5*f4*u5**2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F45r*c5*f4*u5**2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F44r*c4*f4*u4*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F44r*c4*f4*u4*u5*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F43r*c3*f4*u3*
     & u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F43r*c3*f4*u3*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F42r*c2*f4*u2*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F42r*c2*f4*u2*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F41r*c1*f4*u1*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F41r*c1*f4*u1*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F40r*c0*f4*u0*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F45*F40r*c0*f4*u0*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F45r*c5*f4*u4*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F45r*c5*f4*u4*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F44r*c4*f4*u4**2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F44r*c4*f4*u4**2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F43r*c3*f4*u3*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F43r*c3*f4*u3*u4*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F42r*c2*f4*u2*
     & u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F42r*c2*f4*u2*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F41r*c1*f4*u1*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F41r*c1*f4*u1*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F40r*c0*f4*u0*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F44*F40r*c0*f4*u0*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F45r*c5*f4*u3*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F45r*c5*f4*u3*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F44r*c4*f4*u3*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F44r*c4*f4*u3*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F43r*c3*f4*u3**2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F43r*c3*f4*u3**2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F42r*c2*f4*u2*u3*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F42r*c2*f4*u2*u3*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F41r*c1*f4*u1*
     & u3*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F41r*c1*f4*u1*u3*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F40r*c0*f4*u0*u3*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F43*F40r*c0*f4*u0*u3*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F45r*c5*f4*u2*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F45r*c5*f4*u2*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F44r*c4*f4*u2*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F44r*c4*f4*u2*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F43r*c3*f4*u2*u3*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F43r*c3*f4*u2*u3*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F42r*c2*f4*u2**2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F42r*c2*f4*u2**2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F41r*c1*f4*u1*u2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F41r*c1*f4*u1*u2*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F40r*c0*f4*u0*
     & u2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F42*F40r*c0*f4*u0*u2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F45r*c5*f4*u1*u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F45r*c5*f4*u1*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F44r*c4*f4*u1*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F44r*c4*f4*u1*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F43r*c3*f4*u1*u3*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F43r*c3*f4*u1*u3*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F42r*c2*f4*u1*u2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F42r*c2*f4*u1*u2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F41r*c1*f4*u1**2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F41r*c1*f4*u1**2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F40r*c0*f4*u0*u1*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F41*F40r*c0*f4*u0*u1*w1
     &
      traza1 = traza1 - 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F45r*c5*f4*u0*
     & u5*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F45r*c5*f4*u0*u5*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F44r*c4*f4*u0*u4*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F44r*c4*f4*u0*u4*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F43r*c3*f4*u0*u3*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F43r*c3*f4*u0*u3*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F42r*c2*f4*u0*u2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F42r*c2*f4*u0*u2*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F41r*c1*f4*u0*u1*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F41r*c1*f4*u0*u1*w1
     &  - 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F40r*c0*f4*u0**2*w2
     &  + 4.D0*q_uh*i_*PC2*qq*svp*svm*F40*F40r*c0*f4*u0**2*w1
     &



      
   
       if(print.eqv..true.)then
        write(2,*)"traza-variables_____________________________________"       
       write(2,*)"qq",qq
       write(2,*)"F"
c       do  103 i =1,4
c       write(2,*)(F(i,j),j=0,4)
c103    continue
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
