       subroutine kabT (print,w1,pc,p,q,k,svp,ssp,svm,
     .  ssm,g2,alpha,beta,tr)
       implicit none
c       include 'common.f'

       double precision  pp,qq,DG,g2,w1,w2,root5,root3,
     . root2,root6
       Integer alpha,beta,gamma,rho,mu,nu
       double precision  k(4),p(4),q(4),lpv,QM
       double complex    PC(4),p_pc,PC_q,p_q,k_p,k_pc,k_q,PM,
     . kab(8,8),i_,k_pt,k_qt,p_qt,
     . pc_pt,pc_qt,pp2,pq2, q_pt,q_qt,qt_pt,ssm,ssp,svm,svp,
     . kk,pt(4),qt(4),f2,f3,mn
       double complex  sigmasp,sigmasm,sigmavp,sigmavm,delta
       double complex pcp,a1,a2,a3,a4,pc2,tr,k11,z,
     . k12,k13,k14,k15,k16,k17,k18,k21,k22,k23,k24,k25,k26,k27,k28,
     . k31,k32,k33,k34,k35,k36,k37,k38,k41,k42,k43,k44,k45,k46,k47,k48,
     . k51,k52,k53,k54,k55,k56,k57,k58,k61,k62,k63,k64,k65,k66,k67,k68,
     . k71,k72,k73,k74,k75,k76,k77,k78,k81,k82,k83,k84,k85,k86,k87,k88
       integer i,j
       logical print

        i_ = (0,1d0)
        w2=1d0-w1

        root2= dsqrt(2d0)
        root3= dsqrt(3d0)
        root5= dsqrt(5d0)
        root6= dsqrt(6d0)

       pc2 =0d0
       pp  =0d0
       qq  =0d0
       kk  =0d0
       p_pc=0d0
       pc_q=0d0
       p_q =0d0
       k_p =0d0
       k_pc=0d0
       k_q =0d0
       tr = 0d0    

       k_pt= 0d0 !(----------------------- hay que definir estas variables)
       k_qt= 0d0
       mn =  0d0
       p_qt= 0d0
       pc_pt = 0d0
       pc_qt = 0d0
       pp2   = 0d0
       pq2   = 0d0
       q_pt  = 0d0
       q_qt  = 0d0
       qt_pt = 0d0
       z     = 0d0

       pc2=pc(4)*pc(4)
       p_pc=p(4)*pc(4)

       do 102 i = 1,4
c       pc2 =pc(i)*pc(i)+pc2
       pp=p(i)*p(i)    +pp
       qq=q(i)*q(i)   +qq
       kk=k(i)*k(i)    +kk
c       p_pc=p(i)*pc(i)+p_pc
       pc_q=pc(i)*q(i)+pc_q
       p_q =p(i)*q(i) +p_q
       k_p =k(i)*p(i) +k_p
       k_pc=k(i)*pc(i)+k_pc
       k_q =k(i)*q(i) +k_q

c       pt=p(i)-PC(i)*(p.PC)/(PC.PC)
c       qt=q(i)-PC(i)*(q.PC)/(PC.PC)
       pt(i)=p(i)-PC(i)*(p_pc)/(PC2)
       qt(i)=q(i)-PC(i)*(p_pc)/(PC2)

       k_pt=k(i)*pt(i)+k_pt
       k_qt=k(i)*qt(i)+k_qt
       p_qt= p(i)*qt(i)+p_qt
       pc_pt = pc(i)*pt(i)+pc_pt
       pc_qt = pc(i)*qt(i)+pc_qt
       pp2   = (p(i)*p(i))*(p(i)*p(i))
       pq2   = (p(i)*q(i))*(p(i)*q(i))
       q_pt  = q(i)*pt(i)+q_pt
       q_qt  = q(i)*qt(i)+q_qt
       qt_pt = qt(i)*pt(i)+qt_pt
       
       QM = dsqrt(qq)
       mn = pc(4)/i_
       PM = i_*MN

c       do 103 j = 1,4
c       A(i,j) = 0d0
c       kab(i,j)=0d0
c       103    continue
102    continue

c       mn = dsqrt(abs(pc2))
       z=pc_q/(dsqrt(qq)*i_*mn)

       p_pc= p_pc+(0d0,1d0)*1d-14
       kk = kk+ 1d-14
       LPV = 2000d0
       dg = -4d0/3d0*g2*1d0/kk*1d0/(1d0+kk/LPV**2)


        pcp =p_pc !dsqrt(pp)*pc(4)*zp
c       p_pc=pcp


c        f3= dsqrt(abs(4d0/3d0*(1-z*z)))
c        f2= dsqrt(abs(8d0/5d0*(1-z*z)*(1-z*z)))  

         f3 = 4d0/3d0*(1-z*z)
         f2 = 8d0/5d0*(1-z*z)*(1-z*z)
  
c       A(1,1) = -1d0
c       A(2,2) = -PC(4)**2
c       A(2,3) =  -pcp**2
c       A(3,2) =  -pcp**2
c       A(3,3) = -pp*pcp**2
c       A(4,4) =  pp*pc(4)**2-(pcp)**2

c       delta = A(3,3)*A(2,2)-A(3,2)*A(2,3)+1d-14

c       do 100 i = 1,4
c       do 101 j = 1,4
c       Pabi(i,j) = 0d0
c       101    continue
c       100    continue
c       pabi(1,1) = 1d0/A(1,1)
c       pabi(2,2) =  A(3,3)/delta
c       pabi(2,3) = -A(2,3)/delta
c       pabi(3,2) = -A(3,2)/delta
c       pabi(3,3) =  A(2,2)/delta
c       pabi(4,4) =  1d0/A(4,4)



c       P11 = -1d0
c       P22 = -pp*p_PC/delta
c       P23 =  p_PC**2/delta
c       P32 =  p_PC**2/delta
c       P33 = -PC2/delta
c       P44 = 1d0/(pp*PC2-p_PC**2)

c       P11 = -1d0
c       P22 = -pp*p_PC/delta
c       P23 =  p_PC**2/delta
c       P32 =  p_PC**2/delta
c       P33 = -PC2/delta
c       P44 = 1d0/(pp*PC2-p_PC**2)

c       P1 = P11*A1
c       P2 = P22*A2 + P23*A3
c       P3 = P32*A2 + P33*A3
c       P4 = P44*A4

        if((alpha.eq.1).and.(beta.eq.1))then

      K11 = - 6.D0*PC2**(-1)*DG*svp*svm*Pq2
     &  - 7.D0*DG*svp*svm*w1*w2*mn**2
     &  + 5.D0*DG*ssp*ssm
     &  - qq*DG*svp*svm
     &  + 6.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm
     &  - 5.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2
     &  + k_PC*k_q*kk**(-1)*DG*svp*svm*w1
     &  + 4.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2
     &  + 7.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*ssp*ssm
     &  - 5.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*svp*svm
     &  - k_PC**2*kk**(-1)*DG*svp*svm*w1*w2
     &  + 7.D0*PC_q*DG*svp*svm*w2
     &  - 7.D0*PC_q*DG*svp*svm*w1
     &

        elseif((alpha.eq.1).and.(beta.eq.2))then

      K12 = + 24.D0*PC2**(-2)*qq**(-1)*DG*svp*svm*Pq2**2*f2**(-1)*
     . root5**(-1)
     &  - 8.D0*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-1)*root5**(-1)
     &  - 44.D0*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*root5**(-1)
     &  - 4.D0*qq**(-1)*DG*svp*svm*w1*w2*Pq2*f2**(-1)*root5**(-1)
     &  - 4.D0*DG*svp*svm*w1*w2*mn**2*f2**(-1)*root5**(-1)
     &  + 8.D0*DG*ssp*ssm*f2**(-1)*root5**(-1)
     &  + 20.D0*qq*DG*svp*svm*f2**(-1)*root5**(-1)
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*
     & Pq2*f2**(-1)*root5**(-1)
     &  + 48.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*
     & f2**(-1)*root5**(-1)
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-1)*
     & root5**(-1)
     &  + 48.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*root5**(-1)
     &
      K12 = K12 + 8.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w2*Pq2*f2**(-1)*root5**(-1)
     &  - 16.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*
     & f2**(-1)*root5**(-1)
     &  + 40.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*root5**(-1)
     &  - 32.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-1)*root5**(-1)
     &  - 32.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*f2**(-1)*
     & root5**(-1)
     &  + 24.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*f2**(-1)*
     & root5**(-1)
     &  - 16.D0*k_PC**2*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*f2**(-1)*
     & root5**(-1)
     &
      K12 = K12 - 32.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*w2*Pq2*f2**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & root5**(-1)
     &  + 16.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f2**(-1)*
     & root5**(-1)
     &  + 8.D0*k_PC**2*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*root5**(-1)
     &  - 24.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*f2**(-1)*
     & root5**(-1)
     &  + 24.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*f2**(-1)*
     & root5**(-1)
     &  + 48.D0*k_q**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*Pq2*
     & f2**(-1)*root5**(-1)
     &  + 24.D0*k_q**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*root5**(-1)
     &
      K12 = K12 - 24.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*ssm*f2**(-1)*
     & root5**(-1)
     &  - 24.D0*k_q**2*kk**(-1)*DG*svp*svm*f2**(-1)*root5**(-1)
     &  - 4.D0*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & root5**(-1)
     &  + 4.D0*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & root5**(-1)
     &  + 4.D0*PC_q*DG*svp*svm*w2*f2**(-1)*root5**(-1)
     &  - 4.D0*PC_q*DG*svp*svm*w1*f2**(-1)*root5**(-1)
     &

        elseif((alpha.eq.1).and.(beta.eq.3))then

      K13 = + 4.D0*k_PC*k_q*PC_q*PC_qt*kk**(-1)*PC2**(-1)*DG*svp*svm*
     . f3**(-1)*PM**(-1)*QM**(-1)
     &  - 2.D0*k_PC*k_q*PC_qt*kk**(-1)*DG*svp*svm*w2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 2.D0*k_PC*k_q*PC_qt*kk**(-1)*DG*svp*svm*w1*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 2.D0*k_PC*k_qt*PC_q*kk**(-1)*DG*svp*svm*w2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  - 2.D0*k_PC*k_qt*PC_q*kk**(-1)*DG*svp*svm*w1*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  - 4.D0*k_PC*k_qt*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  - 10.D0*k_PC*k_qt*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  - 10.D0*k_PC*k_qt*kk**(-1)*DG*ssp*ssm*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &
      K13 = K13 + 10.D0*k_PC*k_qt*kk**(-1)*qq*DG*svp*svm*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  + 4.D0*k_PC**2*PC_q*q_qt*kk**(-1)*PC2**(-1)*DG*svp*svm*f3**(-1)
     & *PM**(-1)*QM**(-1)
     &  + 4.D0*k_PC**2*PC_qt*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  - 4.D0*k_PC**2*PC_qt*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  - 4.D0*k_PC**2*PC_qt*kk**(-1)*DG*svp*svm*w1*w2*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  - 2.D0*k_PC**2*q_qt*kk**(-1)*DG*svp*svm*w2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 2.D0*k_PC**2*q_qt*kk**(-1)*DG*svp*svm*w1*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  - 16.D0*k_q*k_qt*PC_q*kk**(-1)*DG*svp*svm*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &
      K13 = K13 - 8.D0*k_q*k_qt*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  + 8.D0*k_q*k_qt*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 12.D0*PC_q*q_qt*DG*svp*svm*f3**(-1)*PM**(-1)*QM**(-1)
     &  + 6.D0*PC_qt*DG*svp*svm*w1*w2*mn**2*f3**(-1)*PM**(-1)*QM**(-1)
     &  + 6.D0*PC_qt*DG*ssp*ssm*f3**(-1)*PM**(-1)*QM**(-1)
     &  - 6.D0*PC_qt*qq*DG*svp*svm*f3**(-1)*PM**(-1)*QM**(-1)
     &  + 6.D0*q_qt*DG*svp*svm*w2*mn**2*f3**(-1)*PM**(-1)*QM**(-1)
     &  - 6.D0*q_qt*DG*svp*svm*w1*mn**2*f3**(-1)*PM**(-1)*QM**(-1)
     &

        elseif((alpha.eq.1).and.(beta.eq.4))then

      K14 = + 4.D0*i_*DG*svp*svm*w2*Pq2*f3**(-1)*PM**(-1)*QM**(-1)*
     . root2
     &  + 4.D0*i_*DG*svp*svm*w1*Pq2*f3**(-1)*PM**(-1)*QM**(-1)*root2
     &  + 4.D0*i_*qq*DG*svp*svm*w2*mn**2*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2
     &  + 4.D0*i_*qq*DG*svp*svm*w1*mn**2*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2
     &  - 2.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  + 2.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w1*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  - 6.D0*k_PC*k_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  + 2.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*ssm*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &
      K14 = K14 + 10.D0*k_PC*k_q*i_*kk**(-1)*qq*DG*svp*svm*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  - 2.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  - 2.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2
     &  - 6.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  + 4.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2
     &  - 6.D0*k_PC**2*i_*kk**(-1)*qq*DG*svp*svm*w2*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  + 2.D0*k_PC**2*i_*kk**(-1)*qq*DG*svp*svm*w1*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  - 8.D0*k_q**2*PC_q*i_*kk**(-1)*DG*svp*svm*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &
      K14 = K14 - 4.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  + 4.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2
     &

        elseif((alpha.eq.1).and.(beta.eq.5))then

      K15 = - 6.D0*i_*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-1)*QM**(-1)
     &  + 6.D0*i_*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-1)*QM**(-1)
     &  + 6.D0*i_*qq*DG*ssm*svp*f3**(-1)*QM**(-1)
     &  - 6.D0*i_*qq*DG*ssp*svm*f3**(-1)*QM**(-1)
     &  + 6.D0*k_PC*k_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-1)*
     & QM**(-1)
     &  - 6.D0*k_PC*k_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-1)*
     & QM**(-1)
     &  - 10.D0*k_PC*k_q*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-1)*QM**(-1)
     &  - 10.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-1)*QM**(-1)
     &  + 10.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*
     & f3**(-1)*QM**(-1)
     &  + 10.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*
     & f3**(-1)*QM**(-1)
     &  + 2.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*qq*DG*ssm*svp*f3**(-1)*
     & QM**(-1)
     &
      K15 = K15 - 2.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*qq*DG*ssp*svm*
     & f3**(-1)*QM**(-1)
     &  - 8.D0*k_q**2*i_*kk**(-1)*DG*ssm*svp*f3**(-1)*QM**(-1)
     &  + 8.D0*k_q**2*i_*kk**(-1)*DG*ssp*svm*f3**(-1)*QM**(-1)
     &

        elseif((alpha.eq.1).and.(beta.eq.6))then

      K16 = - 6.D0*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-1)*QM**(-1)*
     . root2**(-1)
     &  - 6.D0*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-1)*QM**(-1)*root2**(-1)
     &  + 6.D0*qq*DG*ssm*svp*f3**(-1)*QM**(-1)*root2**(-1)
     &  + 6.D0*qq*DG*ssp*svm*f3**(-1)*QM**(-1)*root2**(-1)
     &  - 4.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  - 12.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  - 4.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*w1*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*w2*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  + 4.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &
      K16 = K16 - 12.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*
     & f3**(-1)*QM**(-1)*root2**(-1)
     &  - 8.D0*k_PC**2*kk**(-1)*PC2**(-2)*DG*ssm*svp*Pq2*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  + 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*ssm*svp*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  + 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*ssp*svm*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &

        elseif((alpha.eq.1).and.(beta.eq.7))then

      K17 = + 6.D0*qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-1)*PM**(-1)*
     . root2**(-1)*root5**(-1)*root6
     &  - 6.D0*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-1)*PM**(-1)*root2**(-1)
     & *root5**(-1)*root6
     &  + 6.D0*DG*ssm*svp*w1*mn**2*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 12.D0*DG*ssm*svp*w1*mn**2*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  + 12.D0*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  - 6.D0*DG*ssp*svm*w2*mn**2*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 12.D0*DG*ssp*svm*w2*mn**2*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  - 12.D0*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  - 16.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K17 = K17 + 16.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 6.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 6.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 10.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-1)*PM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  - 2.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  + 2.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &  - 10.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-1)*PM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  + 6.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &
      K17 = K17 - 6.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*z**2*f2**(-1)*
     & PM**(-1)*root3*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root3*root5**(-1)
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*z**2*f2**(-1)
     & *PM**(-1)*root3*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 10.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 10.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 2.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*PM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &
      K17 = K17 + 6.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*PM**(-1)
     & *root3*root5**(-1)
     &  - 6.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &  + 2.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*PM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  - 6.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  + 6.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &  + 8.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K17 = K17 + 8.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 6.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*PM**(-1)
     & *root2**(-1)*root5**(-1)*root6
     &  + 6.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*PM**(-1)
     & *root2**(-1)*root5**(-1)*root6
     &  - 6.D0*PC_q*DG*ssm*svp*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 12.D0*PC_q*DG*ssm*svp*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  - 12.D0*PC_q*DG*ssm*svp*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  - 6.D0*PC_q*DG*ssp*svm*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 12.D0*PC_q*DG*ssp*svm*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  - 12.D0*PC_q*DG*ssp*svm*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &

        elseif((alpha.eq.1).and.(beta.eq.8))then

      K18 = - 6.D0*qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-1)*PM**(-1)*
     . root5**(-1)*root6
     &  + 6.D0*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-1)*PM**(-1)*root5**(-1)
     & *root6
     &  - 6.D0*DG*ssm*svp*w1*mn**2*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 6.D0*DG*ssp*svm*w2*mn**2*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 16.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 16.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 10.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-1)*PM**(-1)*
     & root5**(-1)*root6
     &
      K18 = K18 + 10.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-1)*PM**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 10.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 10.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 2.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*PM**(-1)*
     & root5**(-1)*root6
     &  - 2.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*PM**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &
      K18 = K18 - 8.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 8.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 6.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*PM**(-1)
     & *root5**(-1)*root6
     &  - 6.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*PM**(-1)
     & *root5**(-1)*root6
     &  + 6.D0*PC_q*DG*ssm*svp*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 6.D0*PC_q*DG*ssp*svm*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.2).and.(beta.eq.1))then

      K21 = - 8.D0*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*f2**(-1)*
     . root5**(-1)
     &  + 4.D0*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*f2**(-1)*
     & root5**(-1)
     &  - 4.D0*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-1)*root5**(-1)
     &  - 4.D0*DG*svp*svm*w1*w2*mn**2*f2**(-1)*root5**(-1)
     &  + 8.D0*DG*ssp*ssm*f2**(-1)*root5**(-1)
     &  - 4.D0*qq*DG*svp*svm*f2**(-1)*root5**(-1)
     &  + 42.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*w2*f2**(-1)*root5**(-1)
     &  - 18.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*w1*f2**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*
     & Pq2*f2**(-1)*root5**(-1)
     &  + 36.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*
     & f2**(-1)*root5**(-1)
     &
      K21 = K21 - 12.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG
     & *svp*svm*f2**(-1)*root5**(-1)
     &  + 36.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*f2**(-1)*root5**(-1)
     &  - 6.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*f2**(-1)*
     & root5**(-1)
     &  - 18.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*f2**(-1)*
     & root5**(-1)
     &  - 18.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*f2**(-1)*
     & root5**(-1)
     &  + 18.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*f2**(-1)*
     & root5**(-1)
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*
     & f2**(-1)*root5**(-1)
     &
      K21 = K21 + 18.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2
     & *f2**(-1)*root5**(-1)
     &  - 18.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & root5**(-1)
     &  + 6.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)*
     & root5**(-1)
     &  + 12.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*
     & Pp2*f2**(-1)*root5**(-1)
     &  - 12.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-1)*
     & root5**(-1)
     &  - 10.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pp2*
     & f2**(-1)*root5**(-1)
     &  + 2.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pp2*
     & f2**(-1)*root5**(-1)
     &  + 10.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*root5**(-1)
     &
      K21 = K21 - 2.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*
     & svp*svm*f2**(-1)*root5**(-1)
     &  + 6.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w2*f2**(-1)*root5**(-1)
     &  + 18.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w1*f2**(-1)*root5**(-1)
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*w2*
     & Pp2*f2**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*f2**(-1)*
     & root5**(-1)
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*ssm*Pp2*
     & f2**(-1)*root5**(-1)
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*svm*Pp2*
     & f2**(-1)*root5**(-1)
     &
      K21 = K21 - 20.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*w2*Pp2*f2**(-1)*root5**(-1)
     &  - 14.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & root5**(-1)
     &  + 10.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f2**(-1)*
     & root5**(-1)
     &  + 2.D0*k_PC**2*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*root5**(-1)
     &  - 4.D0*PC_q*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pp2*f2**(-1)*
     & root5**(-1)
     &  + 4.D0*PC_q*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pp2*f2**(-1)*
     & root5**(-1)
     &  + 4.D0*PC_q*DG*svp*svm*w2*f2**(-1)*root5**(-1)
     &  - 4.D0*PC_q*DG*svp*svm*w1*f2**(-1)*root5**(-1)
     &

        elseif((alpha.eq.2).and.(beta.eq.2))then

      K22 = - 16.D0*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*Pq2*Pp2*
     . f2**(-2)*root5**(-2)
     &  + 8.D0*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*f2**(-2)*
     & root5**(-2)
     &  - 8.D0*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*Pp2*
     & f2**(-2)*root5**(-2)
     &  + 16.D0*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*root5**(-2)
     &  - 8.D0*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*
     & root5**(-2)
     &  + 8.D0*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)*root5**(-2)
     &  + 16.D0*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*root5**(-2)
     &  - 8.D0*PC2**(-1)*DG*svp*svm*Pq2*f2**(-2)*root5**(-2)
     &  + 8.D0*qq**(-1)*DG*svp*svm*w1*w2*Pq2*f2**(-2)*root5**(-2)
     &  + 8.D0*DG*svp*svm*w1*w2*mn**2*f2**(-2)*root5**(-2)
     &  - 16.D0*DG*ssp*ssm*f2**(-2)*root5**(-2)
     &
      K22 = K22 + 8.D0*qq*DG*svp*svm*f2**(-2)*root5**(-2)
     &  + 48.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)
     & *DG*svp*svm*w2*Pq2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*w2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-3)*qq**(-1)*DG*
     & svp*svm*Pq2**2*f2**(-2)*root5**(-2)
     &  - 144.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*
     & Pq2*f2**(-2)*root5**(-2)
     &  + 96.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm
     & *f2**(-2)*root5**(-2)
     &  + 48.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*svp*svm*Pq2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*f2**(-2)*root5**(-2)
     &
      K22 = K22 - 48.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*f2**(-2)*
     & root5**(-2)
     &  - 48.D0*k_p**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*
     & Pq2**2*f2**(-2)*root5**(-2)
     &  + 96.D0*k_p**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*
     & f2**(-2)*root5**(-2)
     &  - 48.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-2)*
     & root5**(-2)
     &  - 48.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-3)*qq**(-1)*DG*
     & svp*svm*Pq2*Pp2*f2**(-2)*root5**(-2)
     &  + 96.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*
     & ssp*ssm*Pp2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*
     & Pp2*f2**(-2)*root5**(-2)
     &
      K22 = K22 + 96.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*
     & Pq2*f2**(-2)*root5**(-2)
     &  - 96.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*
     & f2**(-2)*root5**(-2)
     &  + 48.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-2)*
     & root5**(-2)
     &  - 96.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & f2**(-2)*root5**(-2)
     &  + 16.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp*
     & svm*w2*Pq2*Pp2*f2**(-2)*root5**(-2)
     &  - 32.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp*
     & svm*w1*Pq2*Pp2*f2**(-2)*root5**(-2)
     &  + 80.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pp2*
     & f2**(-2)*root5**(-2)
     &
      K22 = K22 - 64.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w1*Pp2*f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-2)*root5**(-2)
     &  + 32.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*
     & f2**(-2)*root5**(-2)
     &  - 80.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*f2**(-2)*root5**(-2)
     &  + 64.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*f2**(-2)*root5**(-2)
     &  - 48.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-3)*
     & qq**(-1)*DG*svp*svm*Pq2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*
     & svp*svm*f2**(-2)*root5**(-2)
     &  + 48.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*svp*svm*w2*Pq2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w2*f2**(-2)*root5**(-2)
     &
      K22 = K22 - 32.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-3)*
     & qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*w2*
     & Pp2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*w1*
     & Pp2*f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*root5**(-2)
     &  + 64.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*f2**(-2)*
     & root5**(-2)
     &  - 48.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*f2**(-2)*
     & root5**(-2)
     &  - 32.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-3)*qq**(-1)*DG*ssp*ssm
     & *Pq2*Pp2*f2**(-2)*root5**(-2)
     &  + 64.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-3)*DG*svp*svm*Pq2*Pp2*
     & f2**(-2)*root5**(-2)
     &
      K22 = K22 - 64.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG
     & *svp*svm*w1*w2*Pq2*Pp2*f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*ssm*Pp2*
     & f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*svm*Pp2*
     & f2**(-2)*root5**(-2)
     &  + 16.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*w2*
     & Pp2*f2**(-2)*root5**(-2)
     &  + 32.D0*k_PC**2*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC**2*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*f2**(-2)*
     & root5**(-2)
     &  + 64.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & Pq2*f2**(-2)*root5**(-2)
     &  + 16.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-2)*
     & root5**(-2)
     &
      K22 = K22 - 32.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f2**(-2)*root5**(-2)
     &  - 16.D0*k_PC**2*kk**(-1)*DG*svp*svm*w1*w2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_q**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp
     & *svm*w2*Pp2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_q**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp
     & *svm*w1*Pp2*f2**(-2)*root5**(-2)
     &  + 48.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*f2**(-2)*
     & root5**(-2)
     &  - 48.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*f2**(-2)*
     & root5**(-2)
     &  + 96.D0*k_q**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*
     & Pq2*Pp2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_q**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*
     & Pp2*f2**(-2)*root5**(-2)
     &
      K22 = K22 - 48.D0*k_q**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & Pp2*f2**(-2)*root5**(-2)
     &  - 48.D0*k_q**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pp2*
     & f2**(-2)*root5**(-2)
     &  - 96.D0*k_q**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*Pq2*
     & f2**(-2)*root5**(-2)
     &  - 48.D0*k_q**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*root5**(-2)
     &  + 48.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*ssm*f2**(-2)*
     & root5**(-2)
     &  + 48.D0*k_q**2*kk**(-1)*DG*svp*svm*f2**(-2)*root5**(-2)
     &  - 8.D0*PC_q*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*
     & f2**(-2)*root5**(-2)
     &  + 8.D0*PC_q*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*w1*Pq2*Pp2*
     & f2**(-2)*root5**(-2)
     &
      K22 = K22 + 8.D0*PC_q*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pp2*
     & f2**(-2)*root5**(-2)
     &  - 8.D0*PC_q*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*
     & root5**(-2)
     &  + 8.D0*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*
     & root5**(-2)
     &  - 8.D0*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*
     & root5**(-2)
     &  - 8.D0*PC_q*DG*svp*svm*w2*f2**(-2)*root5**(-2)
     &  + 8.D0*PC_q*DG*svp*svm*w1*f2**(-2)*root5**(-2)
     &

        elseif((alpha.eq.2).and.(beta.eq.3))then

      K23 = + 48.D0*k_p*k_PC*p_PC*PC_q*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)
     . *DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 36.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 36.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 36.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_PC*q_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_PC*p_PC*q_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 + 12.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_q*p_PC*PC_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_q*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_q*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_q*p_qt*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_q*p_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_q*p_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_qt*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 - 12.D0*k_p*k_qt*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm
     & *w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_qt*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_qt*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_qt*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_qt*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 24.D0*k_p**2*PC_q*q_qt*kk**(-1)*pp**(-1)*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 + 12.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p**2*q_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p**2*q_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 24.D0*k_PC*k_q*p_PC*p_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_PC*k_q*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC*k_q*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_PC*k_q*PC_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*
     & svm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC*k_q*PC_q*PC_qt*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 + 8.D0*k_PC*k_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC*k_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1
     & *Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 4.D0*k_PC*k_q*PC_qt*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  - 4.D0*k_PC*k_q*PC_qt*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  - 24.D0*k_PC*k_qt*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC*k_qt*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_PC*k_qt*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2
     & *Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 + 8.D0*k_PC*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 4.D0*k_PC*k_qt*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  + 4.D0*k_PC*k_qt*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC*k_qt*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 20.D0*k_PC*k_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 20.D0*k_PC*k_qt*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 20.D0*k_PC*k_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC*k_qt*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 + 20.D0*k_PC*k_qt*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 20.D0*k_PC*k_qt*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  - 20.D0*k_PC*k_qt*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp
     & *svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_PC**2*PC_q*q_qt*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*
     & svm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*PC_q*q_qt*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 - 16.D0*k_PC**2*PC_qt*kk**(-1)*pp**(-1)*PC2**(-2)*DG*
     & ssp*ssm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC**2*PC_qt*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*svm
     & *Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC**2*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1
     & *w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*PC_qt*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*PC_qt*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*PC_qt*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &
      K23 = K23 + 4.D0*k_PC**2*q_qt*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 4.D0*k_PC**2*q_qt*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  - 32.D0*k_q*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 32.D0*k_q*k_qt*PC_q*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q*k_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_q*k_qt*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q*k_qt*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_q*k_qt*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root5**(-1)
     &

        elseif((alpha.eq.2).and.(beta.eq.4))then

      K24 = - 4.D0*i_*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-1)
     . *f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*i_*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*i_*pp**(-1)*qq*DG*svp*svm*w2*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*i_*pp**(-1)*qq*DG*svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*i_*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root2*root5**(-1)
     &  + 4.D0*i_*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root2*root5**(-1)
     &  + 4.D0*i_*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2*root5**(-1)
     &  + 4.D0*i_*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 - 12.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*
     & root5**(-1)
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG
     & *svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 36.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 - 12.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*
     & svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & mn**2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 - 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & mn**2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*k_PC*k_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*k_PC*k_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*k_PC*k_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 20.D0*k_PC*k_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_PC*k_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 + 12.D0*k_PC*k_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 20.D0*k_PC*k_q*i_*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG
     & *svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*
     & root5**(-1)
     &  - 12.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG
     & *svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*
     & root5**(-1)
     &  - 12.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 12.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 + 12.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*svp*
     & svm*w1*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*ssm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*
     & svm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 24.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 16.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*w2*
     & Pq2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 + 16.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*
     & root5**(-1)
     &  - 8.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*k_PC**2*i_*kk**(-1)*qq*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*k_PC**2*i_*kk**(-1)*qq*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 16.D0*k_q**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 16.D0*k_q**2*PC_q*i_*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_q**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K24 = K24 - 8.D0*k_q**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 8.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2*root5**(-1)
     &

        elseif((alpha.eq.2).and.(beta.eq.5))then

      K25 = + 12.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*
     . DG*ssm*svp*w1*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssp*svm*w2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp
     & *Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*svm
     & *Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssm*
     & svp*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssp*
     & svm*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm
     & *svp*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp
     & *svm*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &
      K25 = K25 - 12.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*w1
     & *f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*ssm*svp*f2**(-1)*
     & f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC*k_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp
     & *Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_PC*k_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*svm
     & *Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &
      K25 = K25 - 12.D0*k_PC*k_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC*k_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*QM**(-1)*root5**(-1)
     &  - 20.D0*k_PC*k_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*w1*
     & Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 20.D0*k_PC*k_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*
     & Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 20.D0*k_PC*k_q*i_*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*f3**(-1)*
     & QM**(-1)*root5**(-1)
     &  + 20.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*f3**(-1)*
     & QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG
     & *ssm*svp*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 12.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG
     & *ssp*svm*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &
      K25 = K25 + 12.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)
     & *DG*ssm*svp*w1*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*w2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*
     & w1*Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*svm*
     & w2*Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 20.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 20.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*ssm*svp*Pp2
     & *f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*ssp*svm*Pp2
     & *f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &
      K25 = K25 - 4.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*qq*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 4.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*qq*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_q**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q**2*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & QM**(-1)*root5**(-1)
     &  - 16.D0*k_q**2*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & QM**(-1)*root5**(-1)
     &

        elseif((alpha.eq.2).and.(beta.eq.6))then

      K26 = + 12.D0*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pq2*Pp2*f2**(-1)*
     . f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 12.D0*pp**(-1)*PC2**(-2)*DG*ssp*svm*Pq2*Pp2*f2**(-1)*f3**(-1)
     & *QM**(-1)*root2**(-1)*root5**(-1)
     &  - 12.D0*pp**(-1)*PC2**(-1)*qq*DG*ssm*svp*Pp2*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  - 12.D0*pp**(-1)*PC2**(-1)*qq*DG*ssp*svm*Pp2*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  - 12.D0*PC2**(-1)*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)
     &  - 12.D0*PC2**(-1)*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)
     &  + 12.D0*qq*DG*ssm*svp*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)
     &  + 12.D0*qq*DG*ssp*svm*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)
     &
      K26 = K26 - 24.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & DG*ssm*svp*w1*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*w2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 72.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 48.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssm*svp
     & *f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 48.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssp*svm
     & *f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K26 = K26 + 24.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_p**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_p**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*ssm*svp*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 8.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*svm*
     & Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K26 = K26 + 8.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_PC*k_q*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 8.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*w1*Pp2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 8.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*
     & ssm*svp*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*
     & ssp*svm*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K26 = K26 - 24.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG
     & *ssm*svp*w1*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm
     & *w2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 32.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*w1*
     & Pp2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 32.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-3)*DG*ssm*svp*Pq2*Pp2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 16.D0*k_PC**2*kk**(-1)*PC2**(-2)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*ssm*svp*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K26 = K26 - 24.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq*DG*ssp*svm*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &

        elseif((alpha.eq.2).and.(beta.eq.7))then

      K27 = - 12.D0*pp**(-1)*DG*ssm*svp*w1*Pp2*f2**(-2)*PM**(-1)*root3*
     . root5**(-2)
     &  + 12.D0*pp**(-1)*DG*ssm*svp*w1*Pp2*z**2*f2**(-2)*PM**(-1)*root3
     & *root5**(-2)
     &  + 12.D0*pp**(-1)*DG*ssp*svm*w2*Pp2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  - 12.D0*pp**(-1)*DG*ssp*svm*w2*Pp2*z**2*f2**(-2)*PM**(-1)*root3
     & *root5**(-2)
     &  - 12.D0*DG*ssm*svp*w1*mn**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 12.D0*DG*ssm*svp*w1*mn**2*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 12.D0*DG*ssp*svm*w2*mn**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*DG*ssp*svm*w2*mn**2*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)
     & *DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     & qq**(-1)*DG*ssp*svm*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*
     & root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 84.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 84.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 36.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 36.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*
     & root5**(-2)*root6
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*svm*w2*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 72.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 72.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 72.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 + 72.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*svm*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*z**2*f2**(-2)
     & *PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 36.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 - 36.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssm
     & *svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp
     & *svm*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 36.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 36.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 + 36.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 36.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 36.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 36.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 + 36.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 32.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*w1*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 32.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*svm*w2*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 32.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 32.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*
     & svp*Pq2*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*
     & svm*Pq2*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 20.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 - 4.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pp2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 4.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*z**2
     & *f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 20.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 20.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-2)*PM**(-1)*
     & root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 + 4.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  - 4.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  + 20.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-2)*PM**(-1)*
     & root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*ssp*svm*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp
     & *f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG
     & *ssm*svp*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp
     & *z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm
     & *f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 36.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm
     & *f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm
     & *z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 32.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 32.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 + 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp
     & *svm*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp*
     & w1*Pq2*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*svm*
     & w2*Pq2*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*w1*Pp2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 - 24.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & w1*Pp2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 24.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*w1*Pp2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 24.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 24.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 20.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 20.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 4.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-2)*PM**(-1)*
     & root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 - 12.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  - 4.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-2)*PM**(-1)*
     & root2**(-1)*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  - 12.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  + 16.D0*k_q**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssm
     & *svp*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp
     & *svm*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K27 = K27 - 16.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssm*svp*w1*Pp2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_q**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*mn**2*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*mn**2*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 12.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 12.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 12.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K27 = K27 + 12.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 12.D0*PC_q*DG*ssm*svp*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*PC_q*DG*ssm*svp*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 12.D0*PC_q*DG*ssp*svm*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*PC_q*DG*ssp*svm*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &

        elseif((alpha.eq.2).and.(beta.eq.8))then

      K28 = - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     . qq**(-1)*DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)
     & *DG*ssp*svm*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*w1*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*svm*w2*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &
      K28 = K28 - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*svm*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssm
     & *svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp
     & *svm*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &
      K28 = K28 + 12.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssm*svp*w1
     & *Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 32.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*w1*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 32.D0*k_PC*k_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*svm*w2*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 32.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 32.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &
      K28 = K28 + 12.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*ssm*svp*Pq2*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*
     & svm*Pq2*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 20.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 20.D0*k_PC*k_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 20.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-2)*PM**(-1)*
     & root5**(-2)*root6
     &  - 20.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-2)*PM**(-1)*
     & root5**(-2)*root6
     &
      K28 = K28 + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     & qq**(-1)*DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*ssp*svm*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp
     & *f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm
     & *f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*svm*Pp2
     & *f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &
      K28 = K28 - 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*
     & ssm*svp*w1*Pq2*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*svm*
     & w2*Pq2*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*w1*Pp2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 20.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 20.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 4.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-2)*PM**(-1)*
     & root5**(-2)*root6
     &  + 4.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-2)*PM**(-1)*
     & root5**(-2)*root6
     &
      K28 = K28 - 16.D0*k_q**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*ssm*svp*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_q**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp
     & *svm*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_q**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssm*svp*w1*Pp2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*mn**2*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*mn**2*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &

        elseif((alpha.eq.3).and.(beta.eq.1))then

      K31 = - 8.D0*k_PC*k_pt*PC_q*kk**(-1)*DG*svp*svm*w2*f3**(-1)*
     . PM**(-1)*QM**(-1)
     &  + 6.D0*k_PC*k_pt*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-1)*
     & PM**(-1)*QM**(-1)
     &  - 6.D0*k_PC*k_pt*kk**(-1)*DG*ssp*ssm*f3**(-1)*PM**(-1)*QM**(-1)
     &  + 2.D0*k_PC*k_pt*kk**(-1)*qq*DG*svp*svm*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 4.D0*k_q*k_pt*PC_q*kk**(-1)*DG*svp*svm*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  - 2.D0*k_q*k_pt*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  - 6.D0*k_q*k_pt*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 8.D0*PC_q*PC_pt*DG*svp*svm*w2*f3**(-1)*PM**(-1)*QM**(-1)
     &  - 4.D0*PC_q*q_pt*DG*svp*svm*f3**(-1)*PM**(-1)*QM**(-1)
     &
      K31 = K31 - 6.D0*PC_pt*DG*svp*svm*w1*w2*mn**2*f3**(-1)*PM**(-1)*
     & QM**(-1)
     &  + 6.D0*PC_pt*DG*ssp*ssm*f3**(-1)*PM**(-1)*QM**(-1)
     &  - 2.D0*PC_pt*qq*DG*svp*svm*f3**(-1)*PM**(-1)*QM**(-1)
     &  + 2.D0*q_pt*DG*svp*svm*w2*mn**2*f3**(-1)*PM**(-1)*QM**(-1)
     &  + 6.D0*q_pt*DG*svp*svm*w1*mn**2*f3**(-1)*PM**(-1)*QM**(-1)
     &

        elseif((alpha.eq.3).and.(beta.eq.2))then

      K32 = - 16.D0*k_PC*k_pt*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*
     . svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC*k_pt*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC*k_pt*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_PC*k_pt*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_q*k_pt*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q*k_pt*PC_q*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q*k_pt*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_q*k_pt*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root5**(-1)
     &
      K32 = K32 + 16.D0*PC_q*PC_pt*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*PC_q*PC_pt*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)
     &  + 16.D0*PC_q*q_pt*PC2**(-1)*qq**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*PC_q*q_pt*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)
     &  - 16.D0*PC_pt*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*PC_pt*qq*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root5**(-1)
     &  - 16.D0*q_pt*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*q_pt*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)
     &

        elseif((alpha.eq.3).and.(beta.eq.3))then

      K33 = - 8.D0*k_PC*k_pt*PC_q*q_qt*kk**(-1)*DG*svp*svm*f3**(-2)*
     . PM**(-2)*QM**(-2)
     &  - 8.D0*k_PC*k_pt*PC_qt*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)
     & *PM**(-2)*QM**(-2)
     &  - 8.D0*k_PC*k_pt*PC_qt*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  + 8.D0*k_PC*k_pt*PC_qt*kk**(-1)*qq*DG*svp*svm*f3**(-2)*PM**(-2)
     & *QM**(-2)
     &  - 4.D0*k_PC*k_pt*q_qt*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &  + 4.D0*k_PC*k_pt*q_qt*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &  - 8.D0*k_q*k_pt*PC_q*PC_qt*kk**(-1)*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &  - 4.D0*k_q*k_pt*PC_qt*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &
      K33 = K33 + 4.D0*k_q*k_pt*PC_qt*kk**(-1)*DG*svp*svm*w1*mn**2*
     & f3**(-2)*PM**(-2)*QM**(-2)
     &  + 4.D0*k_qt*k_pt*PC_q*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &  - 4.D0*k_qt*k_pt*PC_q*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &  + 8.D0*k_qt*k_pt*kk**(-1)*DG*svp*svm*Pq2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  - 4.D0*k_qt*k_pt*kk**(-1)*DG*svp*svm*w1*w2*mn**4*f3**(-2)*
     & PM**(-2)*QM**(-2)
     &  - 4.D0*k_qt*k_pt*kk**(-1)*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  + 4.D0*k_qt*k_pt*kk**(-1)*qq*DG*svp*svm*mn**2*f3**(-2)*PM**(-2)
     & *QM**(-2)
     &  + 8.D0*PC_q*PC_qt*q_pt*DG*svp*svm*f3**(-2)*PM**(-2)*QM**(-2)
     &
      K33 = K33 + 8.D0*PC_q*PC_pt*q_qt*DG*svp*svm*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  - 4.D0*PC_q*qt_pt*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  + 4.D0*PC_q*qt_pt*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  + 8.D0*PC_qt*PC_pt*DG*svp*svm*w1*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  + 8.D0*PC_qt*PC_pt*DG*ssp*ssm*f3**(-2)*PM**(-2)*QM**(-2)
     &  - 8.D0*PC_qt*PC_pt*qq*DG*svp*svm*f3**(-2)*PM**(-2)*QM**(-2)
     &  + 4.D0*PC_qt*q_pt*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  - 4.D0*PC_qt*q_pt*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &  + 4.D0*PC_pt*q_qt*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)
     &
      K33 = K33 - 4.D0*PC_pt*q_qt*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)
     & *QM**(-2)
     &  - 8.D0*qt_pt*DG*svp*svm*Pq2*f3**(-2)*PM**(-2)*QM**(-2)
     &  + 4.D0*qt_pt*DG*svp*svm*w1*w2*mn**4*f3**(-2)*PM**(-2)*QM**(-2)
     &  + 4.D0*qt_pt*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)
     &  - 4.D0*qt_pt*qq*DG*svp*svm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)
     &

        elseif((alpha.eq.3).and.(beta.eq.4))then

      K34 = + 4.D0*k_PC*k_pt*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     . f3**(-2)*PM**(-2)*QM**(-2)*root2
     &  + 4.D0*k_PC*k_pt*PC_q*i_*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*k_PC*k_pt*PC_q*i_*kk**(-1)*qq*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 8.D0*k_PC*k_pt*i_*kk**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 4.D0*k_PC*k_pt*i_*kk**(-1)*qq*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 4.D0*k_PC*k_pt*i_*kk**(-1)*qq*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 4.D0*k_q*k_pt*PC_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 4.D0*k_q*k_pt*PC_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &
      K34 = K34 + 4.D0*k_q*k_pt*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**4*
     & f3**(-2)*PM**(-2)*QM**(-2)*root2
     &  + 4.D0*k_q*k_pt*i_*kk**(-1)*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*k_q*k_pt*i_*kk**(-1)*qq*DG*svp*svm*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 4.D0*PC_q*PC_pt*i_*DG*svp*svm*w1*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 4.D0*PC_q*PC_pt*i_*DG*ssp*ssm*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  - 4.D0*PC_q*PC_pt*i_*qq*DG*svp*svm*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  + 4.D0*PC_q*q_pt*i_*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 4.D0*PC_q*q_pt*i_*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &
      K34 = K34 + 8.D0*PC_pt*i_*DG*svp*svm*w2*Pq2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*PC_pt*i_*qq*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*PC_pt*i_*qq*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 4.D0*q_pt*i_*DG*svp*svm*w1*w2*mn**4*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 4.D0*q_pt*i_*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  - 4.D0*q_pt*i_*qq*DG*svp*svm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &

        elseif((alpha.eq.3).and.(beta.eq.5))then

      K35 = - 4.D0*k_PC*k_pt*PC_q*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     . PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC*k_pt*PC_q*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC*k_pt*i_*kk**(-1)*qq*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 4.D0*k_PC*k_pt*i_*kk**(-1)*qq*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 4.D0*k_q*k_pt*PC_q*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_q*k_pt*PC_q*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_q*k_pt*i_*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_q*k_pt*i_*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &
      K35 = K35 + 4.D0*PC_q*PC_pt*i_*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 4.D0*PC_q*PC_pt*i_*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*PC_q*q_pt*i_*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*PC_q*q_pt*i_*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*PC_pt*i_*qq*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*PC_pt*i_*qq*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*q_pt*i_*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*q_pt*i_*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &

        elseif((alpha.eq.3).and.(beta.eq.6))then

      K36 = + 8.D0*k_PC*k_pt*PC_q*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     . PM**(-1)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC*k_pt*PC_q*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 16.D0*k_PC*k_pt*kk**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC*k_pt*kk**(-1)*qq*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC*k_pt*kk**(-1)*qq*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*k_q*k_pt*PC_q*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 8.D0*k_q*k_pt*PC_q*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 8.D0*k_q*k_pt*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &
      K36 = K36 + 8.D0*k_q*k_pt*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2**(-1)
     &  - 8.D0*PC_q*PC_pt*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*PC_q*PC_pt*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  + 8.D0*PC_q*q_pt*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*PC_q*q_pt*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  - 16.D0*PC_pt*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 8.D0*PC_pt*qq*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  + 8.D0*PC_pt*qq*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &
      K36 = K36 - 8.D0*q_pt*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*q_pt*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &

        elseif((alpha.eq.3).and.(beta.eq.7))then

      K37 = - 16.D0*k_PC*k_pt*PC_q*kk**(-1)*DG*ssm*svp*f2**(-1)*
     . f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 16.D0*k_PC*k_pt*PC_q*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_PC*k_pt*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_pt*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_pt*kk**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_PC*k_pt*kk**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_pt*kk**(-1)*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC*k_pt*kk**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K37 = K37 - 12.D0*k_PC*k_pt*kk**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_pt*kk**(-1)*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_q*k_pt*kk**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_q*k_pt*kk**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_q*k_pt*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_q*k_pt*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_q*k_pt*kk**(-1)*DG*ssm*svp*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_q*k_pt*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K37 = K37 + 12.D0*k_q*k_pt*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_q*k_pt*kk**(-1)*DG*ssp*svm*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 16.D0*PC_q*PC_pt*DG*ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  - 16.D0*PC_q*PC_pt*DG*ssm*svp*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 4.D0*PC_pt*qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*PC_pt*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*PC_pt*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*PC_pt*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &
      K37 = K37 + 12.D0*PC_pt*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*PC_pt*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*PC_pt*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  - 12.D0*PC_pt*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*q_pt*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*q_pt*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*q_pt*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*q_pt*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &
      K37 = K37 - 4.D0*q_pt*DG*ssm*svp*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*q_pt*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*q_pt*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*q_pt*DG*ssp*svm*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &

        elseif((alpha.eq.3).and.(beta.eq.8))then

      K38 = + 4.D0*k_PC*k_pt*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     . f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_pt*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_pt*kk**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_pt*kk**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_q*k_pt*kk**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_q*k_pt*kk**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_q*k_pt*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_q*k_pt*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K38 = K38 - 4.D0*PC_pt*qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*PC_pt*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*PC_pt*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  + 4.D0*PC_pt*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_pt*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_pt*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_pt*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_pt*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.4).and.(beta.eq.1))then

      K41 = - k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f3**(-1)*PM**(-1)
     . *QM**(-1)*root2
     &  - 3.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w1*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  - 2.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  + 2.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  - k_p*k_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  - 3.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2
     &  - 2.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2
     &  - k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w2*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &
      K41 = K41 - 3.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2
     &  + 2.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2
     &  + k_PC**2*p_q*i_*kk**(-1)*DG*svp*svm*w2*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  + 3.D0*k_PC**2*p_q*i_*kk**(-1)*DG*svp*svm*w1*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  - 4.D0*p_PC*PC_q*i_*DG*svp*svm*w2*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2
     &  + 4.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*Pq2*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2
     &  - 4.D0*p_q*PC_q*i_*DG*svp*svm*f3**(-1)*PM**(-1)*QM**(-1)*root2
     &  - 4.D0*p_q*i_*DG*svp*svm*w2*mn**2*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2
     &

        elseif((alpha.eq.4).and.(beta.eq.2))then

      K42 = - 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*
     . svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*svm*
     & Pq2**2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 8.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 8.D0*k_p*k_q*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K42 = K42 + 8.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*
     & qq**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2*root5**(-1)
     &  - 8.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 8.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 8.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*
     & svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 8.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2
     & *Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K42 = K42 - 8.D0*k_PC**2*p_q*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 12.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*p_PC*PC_q*i_*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2*root5**(-1)
     &  + 12.D0*p_PC*PC_q*i_*DG*svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2*root5**(-1)
     &  - 16.D0*p_PC*i_*PC2**(-2)*qq**(-1)*DG*svp*svm*Pq2**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 16.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 16.D0*p_q*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2*root5**(-1)
     &
      K42 = K42 - 16.D0*p_q*PC_q*i_*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*p_q*i_*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  + 12.D0*p_q*i_*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2*root5**(-1)
     &  - 4.D0*p_q*i_*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2*root5**(-1)
     &  + 12.D0*p_q*i_*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2*root5**(-1)
     &

        elseif((alpha.eq.4).and.(beta.eq.3))then

      K43 = - 2.D0*k_p*k_PC*PC_qt*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     . f3**(-2)*PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_p*k_PC*PC_qt*i_*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 2.D0*k_p*k_PC*PC_qt*i_*kk**(-1)*qq*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 4.D0*k_p*k_q*PC_q*PC_qt*i_*kk**(-1)*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_p*k_q*PC_qt*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_p*k_q*PC_qt*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_p*k_qt*PC_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_p*k_qt*PC_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &
      K43 = K43 + 4.D0*k_p*k_qt*i_*kk**(-1)*DG*svp*svm*Pq2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_p*k_qt*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**4*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_p*k_qt*i_*kk**(-1)*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 2.D0*k_p*k_qt*i_*kk**(-1)*qq*DG*svp*svm*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 4.D0*k_PC*k_q*p_qt*PC_q*i_*kk**(-1)*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_PC*k_q*p_qt*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_PC*k_q*p_qt*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_PC*k_qt*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f3**(-2)*PM**(-2)*QM**(-2)*root2
     &
      K43 = K43 - 2.D0*k_PC*k_qt*p_PC*i_*kk**(-1)*DG*ssp*ssm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_PC*k_qt*p_PC*i_*kk**(-1)*qq*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 4.D0*k_PC*k_qt*p_q*PC_q*i_*kk**(-1)*DG*svp*svm*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  - 2.D0*k_PC*k_qt*p_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_PC*k_qt*p_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_PC**2*p_qt*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)
     & *PM**(-2)*QM**(-2)*root2
     &  + 2.D0*k_PC**2*p_qt*i_*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 2.D0*k_PC**2*p_qt*i_*kk**(-1)*qq*DG*svp*svm*f3**(-2)*PM**(-2)
     & *QM**(-2)*root2
     &
      K43 = K43 + 4.D0*p_PC*PC_qt*i_*DG*svp*svm*w1*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2
     &  + 4.D0*p_PC*PC_qt*i_*DG*ssp*ssm*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  - 4.D0*p_PC*PC_qt*i_*qq*DG*svp*svm*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  + 8.D0*p_q*PC_q*PC_qt*i_*DG*svp*svm*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  + 4.D0*p_q*PC_qt*i_*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 4.D0*p_q*PC_qt*i_*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  - 4.D0*p_qt*PC_q*i_*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*p_qt*PC_q*i_*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &
      K43 = K43 - 8.D0*p_qt*i_*DG*svp*svm*Pq2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*p_qt*i_*DG*svp*svm*w1*w2*mn**4*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2
     &  + 4.D0*p_qt*i_*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &  - 4.D0*p_qt*i_*qq*DG*svp*svm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2
     &

        elseif((alpha.eq.4).and.(beta.eq.4))then

      K44 = - 2.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     . f3**(-2)*PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_p*k_PC*PC_q*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  - 2.D0*k_p*k_PC*PC_q*kk**(-1)*qq*DG*svp*svm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  + 2.D0*k_p*k_PC*kk**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  - 2.D0*k_p*k_PC*kk**(-1)*DG*svp*svm*w1*Pq2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  + 2.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*w1*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*w2*mn**4*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &
      K44 = K44 - 2.D0*k_p*k_q*kk**(-1)*DG*ssp*ssm*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_p*k_q*kk**(-1)*qq*DG*svp*svm*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  + 2.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*DG*svp*svm*w2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*DG*svp*svm*w1*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*k_PC*k_q*p_PC*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  - 2.D0*k_PC*k_q*p_PC*kk**(-1)*qq*DG*svp*svm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  - 2.D0*k_PC**2*p_q*PC_q*kk**(-1)*DG*svp*svm*w2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &
      K44 = K44 + 2.D0*k_PC**2*p_q*PC_q*kk**(-1)*DG*svp*svm*w1*f3**(-2)
     & *PM**(-2)*QM**(-2)*root2**2
     &  + 2.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)*
     & PM**(-2)*QM**(-2)*root2**2
     &  + 2.D0*k_PC**2*p_q*kk**(-1)*DG*ssp*ssm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  + 2.D0*k_PC**2*p_q*kk**(-1)*qq*DG*svp*svm*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  - 2.D0*p_PC*PC_q*DG*svp*svm*w1*w2*mn**2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  - 8.D0*p_PC*PC_q*DG*ssp*ssm*f3**(-2)*PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*p_PC*PC_q*qq*DG*svp*svm*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2**2
     &  + 2.D0*p_PC*DG*svp*svm*w2*Pq2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2**2
     &
      K44 = K44 - 2.D0*p_PC*DG*svp*svm*w1*Pq2*f3**(-2)*PM**(-2)*
     & QM**(-2)*root2**2
     &  + 2.D0*p_q*PC_q*DG*svp*svm*w2*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2**2
     &  - 2.D0*p_q*PC_q*DG*svp*svm*w1*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2**2
     &  - 2.D0*p_q*DG*svp*svm*w1*w2*mn**4*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2**2
     &  - 8.D0*p_q*DG*ssp*ssm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*root2**2
     &  - 2.D0*p_q*qq*DG*svp*svm*mn**2*f3**(-2)*PM**(-2)*QM**(-2)*
     & root2**2
     &

        elseif((alpha.eq.4).and.(beta.eq.5))then

      K45 = + 2.D0*k_p*k_PC*PC_q*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     . PM**(-1)*QM**(-2)*root2
     &  + 2.D0*k_p*k_PC*PC_q*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 2.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2
     &  - 2.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2
     &  - 2.D0*k_p*k_q*PC_q*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 2.D0*k_p*k_q*PC_q*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 2.D0*k_p*k_q*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 2.D0*k_p*k_q*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &
      K45 = K45 + 2.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp
     & *f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  - 2.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 2.D0*k_PC*k_q*p_PC*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 2.D0*k_PC*k_q*p_PC*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 2.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2
     &  + 2.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2
     &  - 2.D0*k_PC**2*p_q*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 2.D0*k_PC**2*p_q*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &
      K45 = K45 - 4.D0*p_PC*PC_q*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 4.D0*p_PC*PC_q*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  - 4.D0*p_PC*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*PM**(-1)*QM**(-2)
     & *root2
     &  + 4.D0*p_PC*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*PM**(-1)*QM**(-2)
     & *root2
     &  + 4.D0*p_q*PC_q*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  - 4.D0*p_q*PC_q*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  - 4.D0*p_q*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  - 4.D0*p_q*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &

        elseif((alpha.eq.4).and.(beta.eq.6))then

      K46 = + 4.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     . PM**(-1)*QM**(-2)
     &  + 4.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 4.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 4.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  + 4.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &
      K46 = K46 + 4.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*
     & svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &
      K46 = K46 - 2.D0*p_PC*PC_q*i_*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 10.D0*p_PC*PC_q*i_*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)
     &  - 2.D0*p_PC*i_*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 10.D0*p_PC*i_*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 2.D0*p_q*PC_q*i_*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 10.D0*p_q*PC_q*i_*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)
     &  - 2.D0*p_q*i_*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 10.D0*p_q*i_*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &

        elseif((alpha.eq.4).and.(beta.eq.7))then

      K47 = + 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*
     . svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 6.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssp*svm*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K47 = K47 - 2.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 6.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K47 = K47 + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*
     & root6
     &  + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K47 = K47 - 6.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssp*svm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 2.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 6.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K47 = K47 + 6.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssp*svm*z**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*p_PC*PC_q*i_*DG*ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 2.D0*p_PC*PC_q*i_*DG*ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*p_PC*PC_q*i_*DG*ssm*svp*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*p_PC*PC_q*i_*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  + 6.D0*p_PC*PC_q*i_*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &
      K47 = K47 - 6.D0*p_PC*PC_q*i_*DG*ssp*svm*z**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*p_q*i_*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*p_q*i_*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*p_q*i_*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 2.D0*p_q*i_*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*p_q*i_*DG*ssm*svp*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*p_q*i_*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  + 6.D0*p_q*i_*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &
      K47 = K47 - 6.D0*p_q*i_*DG*ssp*svm*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &

        elseif((alpha.eq.4).and.(beta.eq.8))then

      K48 = - 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*
     . svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*
     . root6
     &  - 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &
      K48 = K48 + 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 2.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &
      K48 = K48 - 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*p_PC*PC_q*i_*DG*ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*p_PC*PC_q*i_*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*p_q*i_*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*p_q*i_*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*p_q*i_*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root5**(-1)*root6
     &
      K48 = K48 - 4.D0*p_q*i_*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &

        elseif((alpha.eq.5).and.(beta.eq.1))then

      K51 = - 4.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     . f3**(-1)*QM**(-1)
     &  - 6.D0*k_p*k_PC*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-1)*QM**(-1)
     &  + 6.D0*k_p*k_PC*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-1)*QM**(-1)
     &  - 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*f3**(-1)*QM**(-1)
     &  - 6.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*f3**(-1)*QM**(-1)
     &  + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-1)*
     & QM**(-1)
     &  + 6.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-1)*
     & QM**(-1)
     &  + 4.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*DG*ssm*svp*
     & f3**(-1)*QM**(-1)
     &  + 6.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*
     & f3**(-1)*QM**(-1)
     &  - 6.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*
     & f3**(-1)*QM**(-1)
     &
      K51 = K51 - 2.D0*p_PC*PC_q*i_*PC2**(-1)*DG*ssm*svp*f3**(-1)*
     & QM**(-1)
     &  - 6.D0*p_PC*PC_q*i_*PC2**(-1)*DG*ssp*svm*f3**(-1)*QM**(-1)
     &  + 2.D0*p_q*i_*DG*ssm*svp*f3**(-1)*QM**(-1)
     &  + 6.D0*p_q*i_*DG*ssp*svm*f3**(-1)*QM**(-1)
     &

        elseif((alpha.eq.5).and.(beta.eq.2))then

      K52 = + 16.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm
     . *svp*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_p*k_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp
     & *Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-3)*qq**(-1)*DG*ssm
     & *svp*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &
      K52 = K52 - 16.D0*p_PC*PC_q*i_*PC2**(-2)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root5**(-1)
     &  + 16.D0*p_PC*PC_q*i_*PC2**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & QM**(-1)*root5**(-1)
     &  + 16.D0*p_q*i_*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root5**(-1)
     &  - 16.D0*p_q*i_*DG*ssm*svp*f2**(-1)*f3**(-1)*QM**(-1)*
     & root5**(-1)
     &

        elseif((alpha.eq.5).and.(beta.eq.3))then

      K53 = + 4.D0*k_p*k_PC*q_qt*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*
     . PM**(-1)*QM**(-2)
     &  + 4.D0*k_p*k_PC*q_qt*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_p*k_q*PC_qt*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_p*k_q*PC_qt*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 4.D0*k_p*k_qt*PC_q*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_p*k_qt*PC_q*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 4.D0*k_p*k_qt*i_*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_p*k_qt*i_*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &
      K53 = K53 + 4.D0*k_PC*k_q*p_PC*PC_qt*i_*kk**(-1)*PC2**(-1)*DG*ssm
     & *svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*PC_qt*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC*k_qt*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*k_PC*k_qt*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC*k_qt*p_PC*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC*k_qt*p_PC*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC**2*p_PC*q_qt*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC**2*p_PC*q_qt*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &
      K53 = K53 - 8.D0*p_PC*PC_q*PC_qt*i_*PC2**(-1)*DG*ssp*svm*f3**(-2)
     & *PM**(-1)*QM**(-2)
     &  + 4.D0*p_PC*PC_qt*i_*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*p_PC*PC_qt*i_*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*p_q*PC_qt*i_*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*p_q*PC_qt*i_*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*p_qt*PC_q*i_*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*p_qt*PC_q*i_*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*p_qt*i_*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*p_qt*i_*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &

        elseif((alpha.eq.5).and.(beta.eq.4))then

      K54 = - 4.D0*k_p*k_PC*PC_q*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     . PM**(-1)*QM**(-2)*root2
     &  - 4.D0*k_p*k_PC*PC_q*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 4.D0*k_p*k_PC*kk**(-1)*qq*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 4.D0*k_p*k_PC*kk**(-1)*qq*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 4.D0*k_p*k_q*PC_q*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 4.D0*k_p*k_q*PC_q*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 4.D0*k_p*k_q*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 4.D0*k_p*k_q*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &
      K54 = K54 - 4.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp
     & *f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  - 4.D0*k_PC*k_q*p_PC*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  - 4.D0*k_PC*k_q*p_PC*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*ssm*svp*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2
     &  - 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*ssp*svm*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2
     &
      K54 = K54 + 4.D0*p_PC*PC_q*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2
     &  + 4.D0*p_PC*PC_q*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*p_PC*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*PM**(-1)*QM**(-2)
     & *root2
     &  - 4.D0*p_PC*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*PM**(-1)*QM**(-2)
     & *root2
     &  - 4.D0*p_q*PC_q*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*p_q*PC_q*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*p_q*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &  + 4.D0*p_q*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*root2
     &

        elseif((alpha.eq.5).and.(beta.eq.5))then

      K55 = - 4.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)
     . *QM**(-2)
     &  - 4.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*
     & QM**(-2)
     &  + 4.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)
     &  + 8.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*
     & QM**(-2)
     &  - 4.D0*k_p*k_PC*kk**(-1)*qq*DG*svp*svm*w2*f3**(-2)*QM**(-2)
     &  - 4.D0*k_p*k_PC*kk**(-1)*qq*DG*svp*svm*w1*f3**(-2)*QM**(-2)
     &  - 4.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*w2*f3**(-2)*QM**(-2)
     &  + 4.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*w1*f3**(-2)*QM**(-2)
     &  + 4.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)*
     & QM**(-2)
     &  + 4.D0*k_p*k_q*kk**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)
     &
      K55 = K55 + 4.D0*k_p*k_q*kk**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)
     &  - 4.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)
     &  - 4.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)*
     & QM**(-2)
     &  - 4.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*
     & QM**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)
     &  + 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*DG*ssp*ssm*f3**(-2)
     & *QM**(-2)
     &  + 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)
     &
      K55 = K55 - 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*w2*f3**(-2)*QM**(-2)
     &  - 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*DG*svp*svm*w2*Pq2*
     & f3**(-2)*QM**(-2)
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)
     &  + 4.D0*p_PC*PC_q*PC2**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)
     &  + 4.D0*p_PC*PC_q*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)
     &  - 4.D0*p_PC*PC_q*DG*svp*svm*w1*w2*f3**(-2)*QM**(-2)
     &  - 4.D0*p_PC*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*QM**(-2)
     &  + 4.D0*p_PC*PC2**(-1)*DG*svp*svm*w1*Pq2*f3**(-2)*QM**(-2)
     &  + 4.D0*p_q*PC_q*DG*svp*svm*w2*f3**(-2)*QM**(-2)
     &  - 4.D0*p_q*PC_q*DG*svp*svm*w1*f3**(-2)*QM**(-2)
     &
      K55 = K55 - 4.D0*p_q*DG*svp*svm*w1*w2*mn**2*f3**(-2)*QM**(-2)
     &  - 4.D0*p_q*DG*ssp*ssm*f3**(-2)*QM**(-2)
     &  - 4.D0*p_q*qq*DG*svp*svm*f3**(-2)*QM**(-2)
     &

        elseif((alpha.eq.5).and.(beta.eq.6))then

      K56 = - 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*
     . f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  + 8.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*w2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*w1*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f3**(-2)*QM**(-2)
     & *root2**(-1)
     &  + 8.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w1*f3**(-2)*QM**(-2)
     & *root2**(-1)
     &  + 16.D0*k_p*k_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &
      K56 = K56 + 8.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*k_p*k_q*i_*kk**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  + 8.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 16.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &
      K56 = K56 + 8.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*DG*ssp*ssm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*w2
     & *f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 16.D0*p_PC*PC_q*i_*PC2**(-2)*DG*svp*svm*Pq2*f3**(-2)*QM**(-2)
     & *root2**(-1)
     &  + 8.D0*p_PC*PC_q*i_*PC2**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &
      K56 = K56 - 8.D0*p_PC*PC_q*i_*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*p_PC*PC_q*i_*DG*svp*svm*w1*w2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  + 8.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*w1*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  + 8.D0*p_q*PC_q*i_*DG*svp*svm*w2*f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*p_q*PC_q*i_*DG*svp*svm*w1*f3**(-2)*QM**(-2)*root2**(-1)
     &  - 16.D0*p_q*i_*PC2**(-1)*DG*svp*svm*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*p_q*i_*DG*svp*svm*w1*w2*mn**2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*p_q*i_*DG*ssp*ssm*f3**(-2)*QM**(-2)*root2**(-1)
     &
      K56 = K56 + 8.D0*p_q*i_*qq*DG*svp*svm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &

        elseif((alpha.eq.5).and.(beta.eq.7))then

      K57 = + 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*
     . svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     . root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K57 = K57 + 4.D0*k_p*k_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 12.D0*k_p*k_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_p*k_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_PC*i_*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_p*k_PC*i_*kk**(-1)*DG*ssp*ssm*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K57 = K57 + 4.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 16.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 16.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &
      K57 = K57 + 4.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 16.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 16.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & z**2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &
      K57 = K57 + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*
     & svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &
      K57 = K57 - 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*qq**(-1)
     & *DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & z**2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K57 = K57 - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*
     & svp*svm*w1*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &
      K57 = K57 + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC**2*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC**2*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*z**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*p_PC*PC_q*i_*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K57 = K57 + 4.D0*p_PC*PC_q*i_*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 4.D0*p_PC*PC_q*i_*DG*svp*svm*w2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*p_PC*PC_q*i_*DG*svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*p_PC*PC_q*i_*DG*svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*p_PC*PC_q*i_*DG*svp*svm*w1*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 16.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 16.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*Pq2*z**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 16.D0*p_q*PC_q*i_*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &
      K57 = K57 - 16.D0*p_q*PC_q*i_*DG*svp*svm*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*p_q*i_*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*p_q*i_*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*p_q*i_*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*p_q*i_*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  - 4.D0*p_q*i_*DG*svp*svm*w2*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*p_q*i_*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*p_q*i_*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &
      K57 = K57 + 12.D0*p_q*i_*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &

        elseif((alpha.eq.5).and.(beta.eq.8))then

      K58 = - 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*
     . svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &
      K58 = K58 + 4.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*i_*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &
      K58 = K58 + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*qq**(-1)
     & *DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &
      K58 = K58 - 4.D0*k_PC**2*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*p_PC*PC_q*i_*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*p_PC*PC_q*i_*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*p_PC*PC_q*i_*DG*svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*p_q*i_*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*p_q*i_*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*p_q*i_*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)*root6
     &
      K58 = K58 - 4.D0*p_q*i_*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.6).and.(beta.eq.1))then

      K61 = + 12.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     . f3**(-1)*QM**(-1)*root2**(-1)
     &  - 4.D0*k_p*k_PC*kk**(-1)*DG*ssm*svp*w1*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  - 4.D0*k_p*k_PC*kk**(-1)*DG*ssp*svm*w2*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  - 4.D0*k_p*k_q*kk**(-1)*DG*ssm*svp*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  - 8.D0*k_p*k_q*kk**(-1)*DG*ssp*svm*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  + 4.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  + 8.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  + 8.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*DG*ssm*svp*f3**(-1)
     & *QM**(-1)*root2**(-1)
     &
      K61 = K61 + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*
     & f3**(-1)*QM**(-1)*root2**(-1)
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  - 8.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  - 12.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-1)*
     & QM**(-1)*root2**(-1)
     &  + 4.D0*p_PC*PC_q*PC2**(-1)*DG*ssp*svm*f3**(-1)*QM**(-1)*
     & root2**(-1)
     &  - 4.D0*p_q*DG*ssp*svm*f3**(-1)*QM**(-1)*root2**(-1)
     &

        elseif((alpha.eq.6).and.(beta.eq.2))then

      K62 = + 24.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*
     . svp*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 24.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 24.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 20.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 4.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 20.D0*k_p*k_PC*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  - 4.D0*k_p*k_PC*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &
      K62 = K62 - 44.D0*k_p*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 20.D0*k_p*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 44.D0*k_p*k_q*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)
     &  - 20.D0*k_p*k_q*kk**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)
     &  + 36.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*
     & svp*w1*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 36.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*
     & svm*w2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 80.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 16.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K62 = K62 - 44.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 20.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 36.D0*k_PC*k_q*p_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*
     & svp*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 36.D0*k_PC*k_q*p_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*
     & svm*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 36.D0*k_PC*k_q*p_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 36.D0*k_PC*k_q*p_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 32.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-3)*qq**(-1)*DG*ssm*
     & svp*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K62 = K62 - 36.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*DG*ssp*svm
     & *f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 16.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp*w1*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 32.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*svm*w2*
     & Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 20.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 36.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp
     & *w1*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 36.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm
     & *w2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 8.D0*k_PC**2*p_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &
      K62 = K62 + 24.D0*k_PC**2*p_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*
     & svm*Pq2*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 28.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 12.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 48.D0*p_PC*PC_q*PC2**(-2)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 56.D0*p_PC*PC_q*PC2**(-2)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  + 48.D0*p_PC*PC_q*PC2**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  - 56.D0*p_PC*PC_q*PC2**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)
     &  + 48.D0*p_q*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)
     & *QM**(-1)*root2**(-1)*root5**(-1)
     &
      K62 = K62 - 56.D0*p_q*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     &  - 48.D0*p_q*DG*ssm*svp*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)
     &  + 56.D0*p_q*DG*ssp*svm*f2**(-1)*f3**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)
     &

        elseif((alpha.eq.6).and.(beta.eq.3))then

      K63 = + 8.D0*k_p*k_PC*PC_q*PC_qt*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     . f3**(-2)*PM**(-1)*QM**(-2)*root2**(-1)
     &  - 4.D0*k_p*k_PC*PC_qt*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 4.D0*k_p*k_PC*PC_qt*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 4.D0*k_p*k_PC*q_qt*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 4.D0*k_p*k_PC*q_qt*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 12.D0*k_p*k_qt*PC_q*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 12.D0*k_p*k_qt*PC_q*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 12.D0*k_p*k_qt*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)
     & *QM**(-2)*root2**(-1)
     &
      K63 = K63 - 12.D0*k_p*k_qt*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC*k_q*p_PC*PC_qt*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC*k_q*p_PC*PC_qt*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2**(-1)
     &  - 4.D0*k_PC*k_q*p_qt*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 4.D0*k_PC*k_q*p_qt*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 28.D0*k_PC*k_qt*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2**(-1)
     &  - 4.D0*k_PC*k_qt*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)*root2**(-1)
     &  - 12.D0*k_PC*k_qt*p_PC*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)
     & *QM**(-2)*root2**(-1)
     &
      K63 = K63 - 12.D0*k_PC*k_qt*p_PC*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2**(-1)
     &  + 16.D0*k_PC*k_qt*p_q*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 16.D0*k_PC*k_qt*p_q*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC**2*p_PC*q_qt*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-2)
     & *PM**(-1)*QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC**2*p_PC*q_qt*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-2)
     & *PM**(-1)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC**2*p_q*PC_qt*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-2)
     & *PM**(-1)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC**2*p_q*PC_qt*kk**(-1)*PC2**(-1)*DG*ssp*svm*f3**(-2)
     & *PM**(-1)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC**2*p_qt*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f3**(-2)
     & *PM**(-1)*QM**(-2)*root2**(-1)
     &
      K63 = K63 + 4.D0*k_PC**2*p_qt*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     & PM**(-1)*QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC**2*p_qt*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 24.D0*p_PC*PC_q*PC_qt*PC2**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*p_PC*PC_q*PC_qt*PC2**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)*root2**(-1)
     &  + 16.D0*p_PC*PC_qt*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  + 16.D0*p_PC*PC_qt*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*p_q*PC_qt*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  - 8.D0*p_q*PC_qt*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &
      K63 = K63 - 16.D0*p_qt*PC_q*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     & *root2**(-1)
     &  + 16.D0*p_qt*PC_q*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  + 16.D0*p_qt*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &  + 16.D0*p_qt*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)*
     & root2**(-1)
     &

        elseif((alpha.eq.6).and.(beta.eq.4))then

      K64 = - 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     . PM**(-1)*QM**(-2)
     &  + 2.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  + 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 6.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 2.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 2.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 2.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &
      K64 = K64 - 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssm*svp*w1*mn**2*f3**(-2)
     & *PM**(-1)*QM**(-2)
     &  + 2.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*svm*w2*mn**2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  - 8.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 8.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  + 2.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*
     & PM**(-1)*QM**(-2)
     &  + 6.D0*k_PC*k_q*p_q*i_*kk**(-1)*DG*ssm*svp*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 6.D0*k_PC*k_q*p_q*i_*kk**(-1)*DG*ssp*svm*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &
      K64 = K64 + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*
     & svm*w2*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 6.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 2.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  - 4.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  + 4.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*svm*
     & f3**(-2)*PM**(-1)*QM**(-2)
     &  + 2.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 10.D0*k_PC**2*p_q*i_*kk**(-1)*DG*ssp*svm*w2*f3**(-2)*PM**(-1)
     & *QM**(-2)
     &  + 16.D0*p_PC*PC_q*i_*DG*ssm*svp*w1*f3**(-2)*PM**(-1)*QM**(-2)
     &
      K64 = K64 + 12.D0*p_PC*PC_q*i_*DG*ssp*svm*w2*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  + 16.D0*p_PC*i_*PC2**(-1)*DG*ssm*svp*Pq2*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 12.D0*p_PC*i_*PC2**(-1)*DG*ssp*svm*Pq2*f3**(-2)*PM**(-1)*
     & QM**(-2)
     &  - 16.D0*p_q*PC_q*i_*DG*ssm*svp*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 12.D0*p_q*PC_q*i_*DG*ssp*svm*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 16.D0*p_q*i_*DG*ssm*svp*w1*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &  + 12.D0*p_q*i_*DG*ssp*svm*w2*mn**2*f3**(-2)*PM**(-1)*QM**(-2)
     &

        elseif((alpha.eq.6).and.(beta.eq.5))then

      K65 = + 12.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*
     . f3**(-2)*QM**(-2)*root2**(-1)
     &  + 12.D0*k_p*k_PC*PC_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 12.D0*k_p*k_PC*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 16.D0*k_p*k_PC*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*w2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 4.D0*k_p*k_PC*i_*kk**(-1)*qq*DG*svp*svm*w1*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  + 12.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &
      K65 = K65 - 12.D0*k_p*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w1*f3**(-2)
     & *QM**(-2)*root2**(-1)
     &  - 12.D0*k_p*k_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 12.D0*k_p*k_q*i_*kk**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 12.D0*k_p*k_q*i_*kk**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 24.D0*k_PC*k_q*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 12.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)
     & *QM**(-2)*root2**(-1)
     &  + 12.D0*k_PC*k_q*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 12.D0*k_PC*k_q*p_PC*i_*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &
      K65 = K65 + 12.D0*k_PC*k_q*p_q*i_*kk**(-1)*DG*svp*svm*w2*f3**(-2)
     & *QM**(-2)*root2**(-1)
     &  + 12.D0*k_PC*k_q*p_q*i_*kk**(-1)*DG*svp*svm*w1*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*DG*ssp*ssm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-2)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 8.D0*k_PC**2*p_PC*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*w2
     & *f3**(-2)*QM**(-2)*root2**(-1)
     &  + 16.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-2)*DG*svp*svm*w2*Pq2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC**2*p_PC*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &
      K65 = K65 - 8.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*
     & svm*w2*f3**(-2)*QM**(-2)*root2**(-1)
     &  - 16.D0*k_PC**2*p_q*PC_q*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-1)
     &  - 4.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 4.D0*k_PC**2*p_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-2)
     & *QM**(-2)*root2**(-1)
     &  + 4.D0*k_PC**2*p_q*i_*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 16.D0*p_PC*PC_q*i_*PC2**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 16.D0*p_PC*PC_q*i_*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  + 16.D0*p_PC*PC_q*i_*DG*svp*svm*w1*w2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &
      K65 = K65 + 16.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-1)
     &  - 16.D0*p_PC*i_*PC2**(-1)*DG*svp*svm*w1*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  - 16.D0*p_q*PC_q*i_*DG*svp*svm*w2*f3**(-2)*QM**(-2)*root2**(-1)
     &  + 16.D0*p_q*PC_q*i_*DG*svp*svm*w1*f3**(-2)*QM**(-2)*root2**(-1)
     &  + 16.D0*p_q*i_*DG*svp*svm*w1*w2*mn**2*f3**(-2)*QM**(-2)*
     & root2**(-1)
     &  + 16.D0*p_q*i_*DG*ssp*ssm*f3**(-2)*QM**(-2)*root2**(-1)
     &  + 16.D0*p_q*i_*qq*DG*svp*svm*f3**(-2)*QM**(-2)*root2**(-1)
     &

        elseif((alpha.eq.6).and.(beta.eq.6))then

      K66 = - 16.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*
     . f3**(-2)*QM**(-2)*root2**(-2)
     &  - 12.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 12.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-2)
     & *QM**(-2)*root2**(-2)
     &  + 4.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 12.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 4.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  - 8.D0*k_p*k_PC*kk**(-1)*qq*DG*svp*svm*w2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  - 8.D0*k_p*k_PC*kk**(-1)*qq*DG*svp*svm*w1*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &
      K66 = K66 - 4.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*w2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 4.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*w1*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  + 16.D0*k_p*k_q*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 4.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 12.D0*k_p*k_q*kk**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  - 12.D0*k_p*k_q*kk**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &  - 4.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &
      K66 = K66 - 16.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2
     & *f3**(-2)*QM**(-2)*root2**(-2)
     &  - 12.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 12.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-2)
     & *QM**(-2)*root2**(-2)
     &  + 4.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  - 8.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*DG*ssp*ssm*f3**(-2)
     & *QM**(-2)*root2**(-2)
     &  + 8.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*qq*DG*svp*svm*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &  - 8.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*w2*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &  - 16.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*DG*svp*svm*w2*Pq2*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &
      K66 = K66 + 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w2
     & *f3**(-2)*QM**(-2)*root2**(-2)
     &  + 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &  + 4.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &  - 4.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & f3**(-2)*QM**(-2)*root2**(-2)
     &  + 16.D0*k_PC**2*p_q*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 20.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  - 20.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  + 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w1*w2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &
      K66 = K66 - 16.D0*p_PC*PC_q*PC2**(-2)*DG*svp*svm*Pq2*f3**(-2)*
     & QM**(-2)*root2**(-2)
     &  - 4.D0*p_PC*PC_q*PC2**(-1)*DG*ssp*ssm*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  + 4.D0*p_PC*PC_q*PC2**(-1)*qq*DG*svp*svm*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  + 12.D0*p_PC*PC_q*DG*svp*svm*w1*w2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  + 12.D0*p_PC*PC2**(-1)*DG*svp*svm*w2*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  - 12.D0*p_PC*PC2**(-1)*DG*svp*svm*w1*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  - 12.D0*p_q*PC_q*DG*svp*svm*w2*f3**(-2)*QM**(-2)*root2**(-2)
     &  + 12.D0*p_q*PC_q*DG*svp*svm*w1*f3**(-2)*QM**(-2)*root2**(-2)
     &  + 16.D0*p_q*PC2**(-1)*DG*svp*svm*Pq2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &
      K66 = K66 + 12.D0*p_q*DG*svp*svm*w1*w2*mn**2*f3**(-2)*QM**(-2)*
     & root2**(-2)
     &  + 4.D0*p_q*DG*ssp*ssm*f3**(-2)*QM**(-2)*root2**(-2)
     &  - 4.D0*p_q*qq*DG*svp*svm*f3**(-2)*QM**(-2)*root2**(-2)
     &

        elseif((alpha.eq.6).and.(beta.eq.7))then

      K67 = + 8.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm
     . *w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     . root5**(-1)*root6
     &  + 16.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1
     & *Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)
     & *root6
     &  - 8.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 12.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &
      K67 = K67 - 12.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  - 4.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 16.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  - 4.D0*k_p*k_PC*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &
      K67 = K67 - 4.D0*k_p*k_PC*kk**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*k_p*k_q*PC_q*kk**(-1)*DG*svp*svm*z**2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*k_p*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 12.D0*k_p*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &
      K67 = K67 - 12.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 12.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K67 = K67 - 12.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*
     & svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 16.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  - 16.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  + 12.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 24.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2
     & *Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)
     & *root6
     &
      K67 = K67 - 12.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 12.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*p_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K67 = K67 + 12.D0*k_PC*k_q*p_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*
     & w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 12.D0*k_PC*k_q*p_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC*k_q*p_q*kk**(-1)*qq**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 12.D0*k_PC*k_q*p_q*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*svp*
     & svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K67 = K67 + 12.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm
     & *w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*
     & w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &
      K67 = K67 - 4.D0*k_PC**2*p_PC*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 12.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  + 12.D0*k_PC**2*p_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K67 = K67 + 8.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*
     & svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 16.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 12.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &
      K67 = K67 + 12.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w1*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  - 16.D0*p_PC*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*p_PC*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 16.D0*p_PC*PC_q*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*p_PC*PC_q*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*p_PC*PC_q*DG*svp*svm*w2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*p_PC*PC_q*DG*svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-2)*root5**(-1)*root6
     &
      K67 = K67 + 12.D0*p_PC*PC_q*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*p_PC*PC_q*DG*svp*svm*w1*z**2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*p_PC*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*p_PC*PC2**(-1)*DG*svp*svm*Pq2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*p_q*PC_q*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root2**(-1)*root3*root5**(-1)
     &  + 16.D0*p_q*PC_q*DG*svp*svm*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*p_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 16.D0*p_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &
      K67 = K67 + 16.D0*p_q*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*p_q*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*p_q*DG*svp*svm*w2*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*p_q*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 12.D0*p_q*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 12.D0*p_q*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2**(-1)*root3*root5**(-1)
     &

        elseif((alpha.eq.6).and.(beta.eq.8))then

      K68 = - 8.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm
     . *w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     . root5**(-1)*root6
     &  - 16.D0*k_p*k_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1
     & *Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     & *root6
     &  + 8.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 16.D0*k_p*k_PC*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p*k_PC*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K68 = K68 + 4.D0*k_p*k_PC*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*kk**(-1)*qq*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K68 = K68 + 12.D0*k_p*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 12.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC*k_q*p_PC*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 24.D0*k_PC*k_q*p_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2
     & *Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     & *root6
     &  + 12.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K68 = K68 + 12.D0*k_PC*k_q*p_PC*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_PC*k_q*p_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC*k_q*p_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC*k_q*p_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 12.D0*k_PC*k_q*p_q*kk**(-1)*qq**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*k_PC*k_q*p_q*kk**(-1)*DG*svp*svm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K68 = K68 + 16.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*
     & DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 12.D0*k_PC**2*p_PC*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-2)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*
     & w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &
      K68 = K68 + 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*DG*ssp*ssm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_PC**2*p_PC*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 12.D0*k_PC**2*p_q*PC_q*kk**(-1)*PC2**(-1)*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*k_PC**2*p_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K68 = K68 - 8.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*
     & svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 16.D0*k_PC**2*p_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 16.D0*p_PC*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 16.D0*p_PC*PC_q*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*p_PC*PC_q*DG*svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K68 = K68 - 16.D0*p_PC*PC_q*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*p_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*p_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*p_q*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*p_q*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.7).and.(beta.eq.1))then

      K71 = - 6.D0*pp**(-1)*DG*ssm*svp*w1*Pp2*f2**(-1)*PM**(-1)*
     . root2**(-1)*root5**(-1)*root6
     &  + 6.D0*pp**(-1)*DG*ssp*svm*w2*Pp2*f2**(-1)*PM**(-1)*root2**(-1)
     & *root5**(-1)*root6
     &  - 6.D0*DG*ssm*svp*w1*mn**2*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 2.D0*DG*ssm*svp*w1*mn**2*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  - 2.D0*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  + 6.D0*DG*ssp*svm*w2*mn**2*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 22.D0*DG*ssp*svm*w2*mn**2*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  + 22.D0*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  + 10.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K71 = K71 + 6.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & DG*ssp*svm*f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 2.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 6.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 6.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K71 = K71 - 6.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 2.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  - 2.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &  - 2.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  + 2.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &  - 2.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 6.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K71 = K71 - 12.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*
     & f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*z**2*
     & f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  - 6.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*w1*Pp2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 6.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-1)*PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 10.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  + 10.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &  + 2.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &  - 2.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*z**2*f2**(-1)*PM**(-1)*
     & root3*root5**(-1)
     &
      K71 = K71 - 6.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*f2**(-1)*
     & PM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 6.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*f2**(-1)*PM**(-1)
     & *root2**(-1)*root5**(-1)*root6
     &  + 6.D0*PC_q*DG*ssm*svp*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 2.D0*PC_q*DG*ssm*svp*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  + 2.D0*PC_q*DG*ssm*svp*z**2*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  + 6.D0*PC_q*DG*ssp*svm*f2**(-1)*PM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 22.D0*PC_q*DG*ssp*svm*f2**(-1)*PM**(-1)*root3*root5**(-1)
     &  + 22.D0*PC_q*DG*ssp*svm*z**2*f2**(-1)*PM**(-1)*root3*
     & root5**(-1)
     &

        elseif((alpha.eq.7).and.(beta.eq.2))then

      K72 = + 32.D0*qq**(-1)*DG*ssm*svp*w1*Pq2*f2**(-2)*PM**(-1)*root3*
     . root5**(-2)
     &  - 32.D0*qq**(-1)*DG*ssm*svp*w1*Pq2*z**2*f2**(-2)*PM**(-1)*root3
     & *root5**(-2)
     &  + 8.D0*qq**(-1)*DG*ssp*svm*w2*Pq2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  - 8.D0*qq**(-1)*DG*ssp*svm*w2*Pq2*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 32.D0*DG*ssm*svp*w1*mn**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 32.D0*DG*ssm*svp*w1*mn**2*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 8.D0*DG*ssp*svm*w2*mn**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 8.D0*DG*ssp*svm*w2*mn**2*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  - 16.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)
     & *DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K72 = K72 + 16.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & DG*ssm*svp*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 16.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 72.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 72.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 72.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 72.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 104.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &
      K72 = K72 + 104.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*
     & svp*Pq2*z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 40.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 40.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 32.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  - 32.D0*k_PC*k_q*kk**(-1)*DG*ssm*svp*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  - 32.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 32.D0*k_PC*k_q*kk**(-1)*DG*ssp*svm*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  - 16.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &
      K72 = K72 + 16.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG
     & *ssm*svp*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-3)*qq**(-1)*DG*
     & ssm*svp*Pq2*Pp2*f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *f2**(-2)*PM**(-1)*root2**(-1)*root5**(-2)*root6
     &  + 48.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp*Pq2
     & *f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 48.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-2)*qq**(-1)*DG*ssm*svp*Pq2
     & *z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 12.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 12.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssm*svp*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K72 = K72 - 36.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*DG*ssp*svm*z**2
     & *f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 16.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 16.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssm*svp*w1*Pq2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 32.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 32.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*svm*w2*Pq2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 20.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  - 20.D0*k_PC**2*kk**(-1)*DG*ssm*svp*w1*z**2*f2**(-2)*PM**(-1)*
     & root3*root5**(-2)
     &  - 4.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &
      K72 = K72 + 4.D0*k_PC**2*kk**(-1)*DG*ssp*svm*w2*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 36.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssm*svp*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 36.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*ssp*svm*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 36.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*mn**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssm*svp*w1*mn**2*z**2*
     & f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 36.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*mn**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &
      K72 = K72 - 36.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*svm*w2*mn**2*
     & z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 32.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 32.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssm*svp*Pq2*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 8.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*f2**(-2)*PM**(-1)
     & *root3*root5**(-2)
     &  + 8.D0*PC_q*PC2**(-1)*qq**(-1)*DG*ssp*svm*Pq2*z**2*f2**(-2)*
     & PM**(-1)*root3*root5**(-2)
     &  - 32.D0*PC_q*DG*ssm*svp*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  + 32.D0*PC_q*DG*ssm*svp*z**2*f2**(-2)*PM**(-1)*root3*
     & root5**(-2)
     &  + 8.D0*PC_q*DG*ssp*svm*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &  - 8.D0*PC_q*DG*ssp*svm*z**2*f2**(-2)*PM**(-1)*root3*root5**(-2)
     &

        elseif((alpha.eq.7).and.(beta.eq.3))then

      K73 = + 8.D0*k_p*k_PC*p_PC*PC_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)
     . *DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*
     . root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_PC*q_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_PC*q_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K73 = K73 + 8.D0*k_p*k_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_qt*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p*k_q*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_q*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K73 = K73 - 4.D0*k_p*k_q*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p*k_q*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_qt*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_qt*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_qt*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_qt*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p**2*PC_q*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K73 = K73 - 4.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & mn**2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 4.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p**2*q_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p**2*q_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K73 = K73 + 4.D0*k_PC*k_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssp*svm*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 4.D0*k_PC*k_qt*p_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_qt*p_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_PC*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 24.D0*k_PC*k_qt*PC_q*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &
      K73 = K73 - 24.D0*k_PC*k_qt*PC_q*kk**(-1)*DG*ssm*svp*z**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_qt*kk**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_qt*kk**(-1)*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_qt*kk**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_qt*kk**(-1)*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 8.D0*k_PC**2*p_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*p_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &
      K73 = K73 + 8.D0*k_PC**2*p_PC*p_qt*PC_q*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K73 = K73 + 12.D0*k_q*k_qt*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_q*k_qt*kk**(-1)*DG*ssm*svp*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_q*k_qt*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_q*k_qt*kk**(-1)*DG*ssp*svm*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  - 8.D0*PC_q*PC_qt*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 24.D0*PC_q*PC_qt*DG*ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 24.D0*PC_q*PC_qt*DG*ssm*svp*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 8.D0*PC_q*PC_qt*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K73 = K73 + 4.D0*PC_qt*pp**(-1)*DG*ssm*svp*w1*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*PC_qt*pp**(-1)*DG*ssp*svm*w2*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*PC_qt*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*PC_qt*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  - 12.D0*PC_qt*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*PC_qt*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 12.D0*PC_qt*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  - 12.D0*PC_qt*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root3*root5**(-1)
     &
      K73 = K73 + 4.D0*q_qt*pp**(-1)*DG*ssm*svp*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*q_qt*pp**(-1)*DG*ssp*svm*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*q_qt*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*q_qt*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*q_qt*DG*ssm*svp*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 4.D0*q_qt*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*q_qt*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*q_qt*DG*ssp*svm*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root3*root5**(-1)
     &

        elseif((alpha.eq.7).and.(beta.eq.4))then

      K74 = - 4.D0*i_*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*Pp2*f2**(-1)*
     . f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*pp**(-1)*qq*DG*ssm*svp*Pp2*f2**(-1)*f3**(-1)*PM**(-2)
     & *QM**(-1)*root5**(-1)*root6
     &  - 4.D0*i_*pp**(-1)*qq*DG*ssp*svm*Pp2*f2**(-1)*f3**(-1)*PM**(-2)
     & *QM**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root5**(-1)*root6
     &  - 20.D0*i_*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root2*root3*root5**(-1)
     &  + 20.D0*i_*DG*ssm*svp*Pq2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*i_*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root5**(-1)*root6
     &
      K74 = K74 + 16.D0*i_*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  - 16.D0*i_*DG*ssp*svm*Pq2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*i_*qq*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 20.D0*i_*qq*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  + 20.D0*i_*qq*DG*ssm*svp*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*i_*qq*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  + 16.D0*i_*qq*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &  - 16.D0*i_*qq*DG*ssp*svm*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root3*root5**(-1)
     &
      K74 = K74 - 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssm*
     & svp*w1*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*ssp*svm*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K74 = K74 - 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & mn**2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*ssm*svp*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*ssp*svm*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 8.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*ssm*svp*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K74 = K74 + 8.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 8.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*ssp*svm*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*k_PC*k_q*i_*kk**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*k_PC*k_q*i_*kk**(-1)*DG*ssm*svp*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*svm*w2*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K74 = K74 - 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp
     & *w1*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & w1*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & w2*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*ssm*svp*w1*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*ssm*svp*w1*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*ssp*svm*w2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*ssp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K74 = K74 + 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssm
     & *svp*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssp*svm*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_PC**2*i_*kk**(-1)*qq*DG*ssm*svp*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*k_PC**2*i_*kk**(-1)*qq*DG*ssm*svp*z**2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 2.D0*k_PC**2*i_*kk**(-1)*qq*DG*ssp*svm*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 2.D0*k_PC**2*i_*kk**(-1)*qq*DG*ssp*svm*z**2*f2**(-1)*f3**(-1)
     & *PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  + 6.D0*k_q**2*i_*kk**(-1)*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 6.D0*k_q**2*i_*kk**(-1)*DG*ssm*svp*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &
      K74 = K74 + 6.D0*k_q**2*i_*kk**(-1)*DG*ssp*svm*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &  - 6.D0*k_q**2*i_*kk**(-1)*DG*ssp*svm*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root3*root5**(-1)
     &

        elseif((alpha.eq.7).and.(beta.eq.5))then

      K75 = - 4.D0*i_*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-1)
     . *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*i_*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*pp**(-1)*qq*DG*svp*svm*w2*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*pp**(-1)*qq*DG*svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root2**(-1)*root5**(-1)*root6
     &  - 12.D0*i_*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*i_*DG*svp*svm*w2*Pq2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  + 4.D0*i_*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root2**(-1)*root5**(-1)*root6
     &
      K75 = K75 - 12.D0*i_*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root3*root5**(-1)
     &  + 12.D0*i_*DG*svp*svm*w1*Pq2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  + 4.D0*i_*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*i_*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*i_*qq*DG*svp*svm*w2*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*i_*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 12.D0*i_*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root3*root5**(-1)
     &  + 12.D0*i_*qq*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &
      K75 = K75 - 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 12.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &
      K75 = K75 - 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm
     & *w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K75 = K75 + 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*ssm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K75 = K75 + 36.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 36.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w1*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*PC_q*i_*kk**(-1)*DG*svp*svm*w1*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_q*i_*kk**(-1)*DG*svp*svm*w1*w2*mn**2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_q*i_*kk**(-1)*DG*ssp*ssm*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root3*root5**(-1)
     &
      K75 = K75 - 12.D0*k_PC*k_q*i_*kk**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC*k_q*i_*kk**(-1)*qq*DG*svp*svm*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &
      K75 = K75 + 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm
     & *w1*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*ssm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*
     & svm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 12.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*DG*ssp*ssm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &
      K75 = K75 + 12.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*
     & svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC**2*PC_q*i_*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_PC**2*PC_q*i_*kk**(-1)*DG*svp*svm*w1*w2*z**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 8.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*w2*Pq2
     & *Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)
     & *root6
     &  + 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K75 = K75 - 24.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 24.D0*k_PC**2*i_*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w2*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  + 12.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &  - 12.D0*k_q**2*i_*kk**(-1)*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root3*root5**(-1)
     &

        elseif((alpha.eq.7).and.(beta.eq.6))then

      K76 = - 8.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-1)*
     . f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 8.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 8.D0*pp**(-1)*qq*DG*svp*svm*w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 8.D0*pp**(-1)*qq*DG*svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 8.D0*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-2)*root5**(-1)*root6
     &  - 36.D0*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root3*root5**(-1)
     &  + 36.D0*DG*svp*svm*w2*Pq2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 8.D0*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-2)*root5**(-1)*root6
     &
      K76 = K76 - 4.D0*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*DG*svp*svm*w1*Pq2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 8.D0*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 36.D0*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 36.D0*qq*DG*svp*svm*w2*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 8.D0*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 4.D0*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 4.D0*qq*DG*svp*svm*w1*mn**2*z**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root3*root5**(-1)
     &
      K76 = K76 - 16.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     & DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 24.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K76 = K76 - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 8.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K76 = K76 - 16.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 16.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 16.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &
      K76 = K76 - 8.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 8.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 8.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 8.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &
      K76 = K76 - 16.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*w1*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  + 16.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*k_PC*k_q*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  + 16.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*w2*mn**2*z**2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 8.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &
      K76 = K76 + 8.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-2)*root5**(-1)*root6
     &  + 16.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm
     & *Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)
     & *root6
     &  + 8.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 8.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*
     & svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 8.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &
      K76 = K76 - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*
     & ssm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*
     & root5**(-1)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*qq*DG*svp*svm*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w1*w2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &
      K76 = K76 - 16.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w1*w2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  + 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  - 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*
     & root6
     &  + 16.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 16.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*z**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*
     & root5**(-1)
     &  - 16.D0*k_PC**2*kk**(-1)*qq*DG*svp*svm*w1*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &
      K76 = K76 + 16.D0*k_PC**2*kk**(-1)*qq*DG*svp*svm*w1*z**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*PC_q*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*PC_q*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  - 16.D0*PC_q*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-2)*root5**(-1)*root6
     &  + 32.D0*PC_q*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  - 32.D0*PC_q*PC2**(-1)*DG*svp*svm*Pq2*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &  + 16.D0*PC_q*qq*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-2)*root5**(-1)*root6
     &  - 32.D0*PC_q*qq*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root3*root5**(-1)
     &
      K76 = K76 + 32.D0*PC_q*qq*DG*svp*svm*z**2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root2**(-1)*root3*root5**(-1)
     &

        elseif((alpha.eq.7).and.(beta.eq.7))then

      K77 = + 4.D0*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*Pp2*
     . f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  - 16.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*f2**(-2)*PM**(-2)
     & *root2**(-1)*root3*root5**(-2)*root6
     &  + 16.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*Pp2*f2**(-2)*
     & PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 - 4.D0*pp**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*pp**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root3*root5**(-2)*root6
     &  - 12.D0*pp**(-1)*DG*ssp*ssm*Pp2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*PM**(-2)*root2**(-2)
     & *root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*PM**(-2)*root2**(-1)
     & *root3*root5**(-2)*root6
     &  - 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 12.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*
     & root5**(-2)*root6**2
     &  + 12.D0*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root3*root5**(-2)*root6
     &  - 12.D0*qq**(-1)*DG*ssp*ssm*Pq2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)
     & *root6**2
     &  + 4.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root3*
     & root5**(-2)*root6
     &  - 32.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root3**2*root5**(-2)
     &  - 4.D0*DG*svp*svm*Pq2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*
     & root5**(-2)*root6
     &
      K77 = K77 + 64.D0*DG*svp*svm*Pq2*z**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  - 32.D0*DG*svp*svm*Pq2*z**4*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  + 4.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root2**(-2)*
     & root5**(-2)*root6**2
     &  - 24.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root2**(-1)*
     & root3*root5**(-2)*root6
     &  + 24.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  + 24.D0*DG*svp*svm*w1*w2*mn**4*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  - 48.D0*DG*svp*svm*w1*w2*mn**4*z**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  + 24.D0*DG*svp*svm*w1*w2*mn**4*z**4*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &
      K77 = K77 - 4.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root2**(-2)*
     & root5**(-2)*root6**2
     &  + 24.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*
     & root5**(-2)*root6
     &  - 24.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root3**2*root5**(-2)
     &  - 24.D0*DG*ssp*ssm*mn**2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root3*root5**(-2)*root6
     &  + 48.D0*DG*ssp*ssm*mn**2*z**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  - 24.D0*DG*ssp*ssm*mn**2*z**4*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  + 4.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root2**(-2)*
     & root5**(-2)*root6**2
     &  - 8.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*
     & root5**(-2)*root6
     &
      K77 = K77 - 8.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  + 8.D0*qq*DG*svp*svm*mn**2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root3*root5**(-2)*root6
     &  + 16.D0*qq*DG*svp*svm*mn**2*z**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  - 8.D0*qq*DG*svp*svm*mn**2*z**4*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)
     & *DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*
     & root6**2
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*
     & DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*
     & root6**2
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &
      K77 = K77 - 20.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm
     & *w2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 20.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*z**2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*z**2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 - 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*Pq2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)
     & *root6
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & w2*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*z**2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-2)*
     & PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 16.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*f2**(-2)
     & *PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 16.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &
      K77 = K77 + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & mn**2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*z**2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &
      K77 = K77 + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & mn**2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2
     & *Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &
      K77 = K77 - 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*f2**(-2)
     & *PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 16.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 16.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*f2**(-2)
     & *PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*ssm*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 24.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 16.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 32.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*z**2*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  - 16.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*z**4*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  - 36.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 36.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 - 12.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*
     & z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)
     & *root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 8.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  - 4.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*z**4*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  - 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)
     & *root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)
     & *root3**2*root5**(-2)
     &
      K77 = K77 + 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 24.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*z**4*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  - 16.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 16.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*
     & root6
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*
     & root6**2
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*
     & root6**2
     &
      K77 = K77 + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2
     & *f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp
     & *svm*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*
     & root6**2
     &
      K77 = K77 - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*w2*Pp2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 24.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 24.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 - 16.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*z**2*
     & f2**(-2)*PM**(-2)*root3**2*root5**(-2)
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*z**4*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  - 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w1*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*
     & Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &
      K77 = K77 + 12.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*
     & Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*z**2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)*
     & PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)
     & *PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*z**2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*z**2*f2**(-2)
     & *PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*f2**(-2)
     & *PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*DG*svp*svm*w1*w2*mn**2*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &  + 24.D0*k_PC**2*kk**(-1)*DG*svp*svm*w1*w2*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root3**2*root5**(-2)
     &
      K77 = K77 - 12.D0*k_PC**2*kk**(-1)*DG*svp*svm*w1*w2*mn**2*z**4*
     & f2**(-2)*PM**(-2)*root3**2*root5**(-2)
     &  + 12.D0*k_PC**2*kk**(-1)*DG*ssp*ssm*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  - 24.D0*k_PC**2*kk**(-1)*DG*ssp*ssm*z**2*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 12.D0*k_PC**2*kk**(-1)*DG*ssp*ssm*z**4*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 4.D0*k_PC**2*kk**(-1)*qq*DG*svp*svm*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  - 8.D0*k_PC**2*kk**(-1)*qq*DG*svp*svm*z**2*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 4.D0*k_PC**2*kk**(-1)*qq*DG*svp*svm*z**4*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  - 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*
     & mn**2*z**2*f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**4*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**4*z**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*ssm*mn**2*z**2*f2**(-2)
     & *PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*kk**(-1)*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 - 12.D0*k_q**2*kk**(-1)*DG*svp*svm*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-2)*root5**(-2)*root6**2
     &  - 4.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &
      K77 = K77 + 12.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*z**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  + 12.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  - 12.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 4.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-2)*root5**(-2)*root6**2
     &  - 12.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 12.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  - 4.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root2**(-2)*
     & root5**(-2)*root6**2
     &
      K77 = K77 + 24.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  - 24.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &  - 24.D0*PC_q*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root3*root5**(-2)*root6
     &  + 48.D0*PC_q*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  - 24.D0*PC_q*DG*svp*svm*w2*mn**2*z**4*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 4.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*root2**(-2)*
     & root5**(-2)*root6**2
     &  - 24.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root3*root5**(-2)*root6
     &  + 24.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*root3**2*
     & root5**(-2)
     &
      K77 = K77 + 24.D0*PC_q*DG*svp*svm*w1*mn**2*z**2*f2**(-2)*PM**(-2)
     & *root2**(-1)*root3*root5**(-2)*root6
     &  - 48.D0*PC_q*DG*svp*svm*w1*mn**2*z**2*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &  + 24.D0*PC_q*DG*svp*svm*w1*mn**2*z**4*f2**(-2)*PM**(-2)*
     & root3**2*root5**(-2)
     &

        elseif((alpha.eq.7).and.(beta.eq.8))then

      K78 = - 4.D0*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*Pp2*
     . f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*Pp2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*PM**(-2)*root2**(-1)
     & *root5**(-2)*root6**2
     &  - 4.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &
      K78 = K78 - 12.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 12.D0*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 12.D0*qq**(-1)*DG*ssp*ssm*Pq2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)
     & *root6**2
     &  + 12.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root3*root5**(-2)*
     & root6
     &  - 12.D0*DG*svp*svm*Pq2*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)
     & *root6
     &  - 4.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &
      K78 = K78 + 12.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 12.D0*DG*svp*svm*w1*w2*mn**4*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 4.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 12.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*
     & root6
     &  + 12.D0*DG*ssp*ssm*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  + 12.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)
     & *root6
     &  - 12.D0*qq*DG*svp*svm*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &
      K78 = K78 + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*
     & DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & w2*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K78 = K78 + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2
     & *mn**2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K78 = K78 + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*
     & svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2
     & *Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*f2**(-2)
     & *PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K78 = K78 + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2*
     & z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*qq**(-1)*DG*ssp*ssm*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 24.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  + 24.D0*k_PC*k_q*PC_q*kk**(-1)*DG*svp*svm*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 36.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &
      K78 = K78 - 36.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*
     & z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)
     & *root3*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)
     & *root3*root5**(-2)*root6
     &  - 12.D0*k_PC*k_q*kk**(-1)*DG*svp*svm*w1*mn**2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &
      K78 = K78 + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp
     & *svm*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 24.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 24.D0*k_PC**2*PC_q*kk**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &
      K78 = K78 + 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w1*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*PC_q*kk**(-1)*DG*svp*svm*w1*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*
     & Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K78 = K78 - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*
     & svm*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*z**2
     & *f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*PC2**(-1)*DG*svp*svm*Pq2*z**2*f2**(-2)
     & *PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*f2**(-2)
     & *PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &
      K78 = K78 + 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*
     & mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w2*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*PC_q*kk**(-1)*qq**(-1)*DG*svp*svm*w1*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**4*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**4*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_q**2*kk**(-1)*qq**(-1)*DG*ssp*ssm*mn**2*z**2*f2**(-2)
     & *PM**(-2)*root3*root5**(-2)*root6
     &
      K78 = K78 - 12.D0*k_q**2*kk**(-1)*DG*svp*svm*mn**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_q**2*kk**(-1)*DG*svp*svm*mn**2*z**2*f2**(-2)*PM**(-2)
     & *root3*root5**(-2)*root6
     &  - 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root3
     & *root5**(-2)*root6
     &
      K78 = K78 + 12.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root3
     & *root5**(-2)*root6
     &  - 12.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*z**2*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  + 4.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 12.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 12.D0*PC_q*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &
      K78 = K78 + 12.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  - 12.D0*PC_q*DG*svp*svm*w1*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.1))then

      K81 = + 6.D0*pp**(-1)*DG*ssm*svp*w1*Pp2*f2**(-1)*PM**(-1)*
     . root5**(-1)*root6
     &  - 6.D0*pp**(-1)*DG*ssp*svm*w2*Pp2*f2**(-1)*PM**(-1)*root5**(-1)
     & *root6
     &  + 6.D0*DG*ssm*svp*w1*mn**2*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  - 6.D0*DG*ssp*svm*w2*mn**2*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  - 10.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  - 6.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & svm*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  - 12.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  + 12.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 2.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &
      K81 = K81 - 6.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  - 6.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*f2**(-1)*
     & PM**(-1)*root5**(-1)*root6
     &  + 2.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 6.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &
      K81 = K81 + 6.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & w1*Pp2*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  - 6.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*w2*Pp2*
     & f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  + 6.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*f2**(-1)*PM**(-1)
     & *root5**(-1)*root6
     &  + 6.D0*PC_q*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*f2**(-1)*PM**(-1)
     & *root5**(-1)*root6
     &  - 6.D0*PC_q*DG*ssm*svp*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &  - 6.D0*PC_q*DG*ssp*svm*f2**(-1)*PM**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.2))then

      K82 = + 16.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     . qq**(-1)*DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*
     & svp*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssm*svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-2)*
     & PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*
     & DG*ssm*svp*Pq2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp
     & *f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  - 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-3)*qq**(-1)*DG*
     & ssm*svp*Pq2*Pp2*f2**(-2)*PM**(-1)*root5**(-2)*root6
     &  + 16.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssm*svp*Pp2
     & *f2**(-2)*PM**(-1)*root5**(-2)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.3))then

      K83 = - 8.D0*k_p*k_PC*p_PC*PC_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)
     . *DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*q_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*q_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_qt*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K83 = K83 - 4.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & mn**2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*p_PC*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p*k_q*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_qt*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_qt*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K83 = K83 - 4.D0*k_p*k_qt*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_qt*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p**2*PC_q*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*PC_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*q_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*q_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*mn**2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_q*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K83 = K83 + 4.D0*k_PC*k_q*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_qt*p_PC*p_q*kk**(-1)*pp**(-1)*DG*ssm*svp*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC*k_qt*p_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*svm*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC*k_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_PC**2*p_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K83 = K83 + 8.D0*k_PC**2*p_PC*p_q*PC_qt*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*p_PC*p_qt*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_qt*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*q_qt*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*PC_q*PC_qt*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root5**(-1)*root6
     &
      K83 = K83 - 8.D0*PC_q*PC_qt*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)
     & *QM**(-1)*root5**(-1)*root6
     &  - 4.D0*PC_qt*pp**(-1)*DG*ssm*svp*w1*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*PC_qt*pp**(-1)*DG*ssp*svm*w2*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-2)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*PC_qt*DG*ssm*svp*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*PC_qt*DG*ssp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_qt*pp**(-1)*DG*ssm*svp*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_qt*pp**(-1)*DG*ssp*svm*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*q_qt*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root5**(-1)*root6
     &
      K83 = K83 - 4.D0*q_qt*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)
     & *QM**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.4))then

      K84 = + 4.D0*i_*pp**(-1)*PC2**(-1)*DG*ssm*svp*Pq2*Pp2*f2**(-1)*
     . f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*i_*pp**(-1)*PC2**(-1)*DG*ssp*svm*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*i_*pp**(-1)*qq*DG*ssm*svp*Pp2*f2**(-1)*f3**(-1)*PM**(-2)
     & *QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*i_*pp**(-1)*qq*DG*ssp*svm*Pp2*f2**(-1)*f3**(-1)*PM**(-2)
     & *QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*i_*DG*ssm*svp*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root2*root5**(-1)*root6
     &  + 4.D0*i_*DG*ssp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*
     & root2*root5**(-1)*root6
     &  - 4.D0*i_*qq*DG*ssm*svp*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*i_*qq*DG*ssp*svm*mn**2*f2**(-1)*f3**(-1)*PM**(-2)*
     & QM**(-1)*root2*root5**(-1)*root6
     &
      K84 = K84 + 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssm*
     & svp*w1*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*
     & root6
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*ssp*svm*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &
      K84 = K84 + 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*ssp*
     & svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*ssm*svp*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*ssp*svm*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & ssm*svp*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*
     & root6
     &
      K84 = K84 - 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*ssp*svm*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*
     & root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssm*svp*w1*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssm*svp*
     & w1*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*
     & root6
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*svm*
     & w2*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*
     & root6
     &  - 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssm*svp*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*root6
     &
      K84 = K84 + 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*ssp
     & *svm*Pp2*f2**(-1)*f3**(-1)*PM**(-2)*QM**(-1)*root2*root5**(-1)*
     & root6
     &

        elseif((alpha.eq.8).and.(beta.eq.5))then

      K85 = + 4.D0*i_*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-1)
     . *f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*i_*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*i_*pp**(-1)*qq*DG*svp*svm*w2*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*i_*pp**(-1)*qq*DG*svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*
     & PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*i_*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root5**(-1)*root6
     &  - 4.D0*i_*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)
     & *root5**(-1)*root6
     &  - 4.D0*i_*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)*root6
     &  - 4.D0*i_*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root5**(-1)*root6
     &
      K85 = K85 + 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*
     & PC2**(-1)*DG*ssp*ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root5**(-1)*root6
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 12.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm
     & *w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &
      K85 = K85 + 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*
     & svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_p*k_PC*p_q*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &
      K85 = K85 + 4.D0*k_p**2*i_*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & mn**2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*i_*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*ssm*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*svp*
     & svm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &
      K85 = K85 + 4.D0*k_PC**2*PC_q*i_*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*w1*w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root5**(-1)*root6
     &  + 8.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*w2*Pq2
     & *Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &  - 4.D0*k_PC**2*i_*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*
     & Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.6))then

      K86 = + 8.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-1)*
     . f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*pp**(-1)*qq*DG*svp*svm*w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*pp**(-1)*qq*DG*svp*svm*w1*Pp2*f2**(-1)*f3**(-1)*PM**(-1)
     & *QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*DG*svp*svm*w2*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  + 8.D0*DG*svp*svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &  - 8.D0*qq*DG*svp*svm*w2*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*qq*DG*svp*svm*w1*mn**2*f2**(-1)*f3**(-1)*PM**(-1)*
     & QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K86 = K86 + 16.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     & DG*svp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*
     & ssm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 24.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*
     & svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K86 = K86 + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*w1*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K86 = K86 + 16.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  + 8.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 16.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K86 = K86 + 8.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w2*Pq2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 8.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*Pq2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w2*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 8.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*w1*mn**2*f2**(-1)
     & *f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 8.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*w2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp
     & *svm*w1*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &
      K86 = K86 - 16.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG
     & *svp*svm*Pq2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*
     & svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*DG*ssp*ssm*Pp2*
     & f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &
      K86 = K86 + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq*DG*
     & svp*svm*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w1*
     & w2*Pp2*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*
     & root5**(-1)*root6
     &  - 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w2*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  + 8.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*w1*Pp2
     & *f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*
     & root6
     &  - 16.D0*PC_q*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  + 16.D0*PC_q*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*f2**(-1)*
     & f3**(-1)*PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &
      K86 = K86 + 16.D0*PC_q*PC2**(-1)*DG*svp*svm*Pq2*f2**(-1)*f3**(-1)
     & *PM**(-1)*QM**(-1)*root2**(-1)*root5**(-1)*root6
     &  - 16.D0*PC_q*qq*DG*svp*svm*f2**(-1)*f3**(-1)*PM**(-1)*QM**(-1)*
     & root2**(-1)*root5**(-1)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.7))then

      K87 = - 4.D0*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*Pp2*
     . f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 16.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*f2**(-2)*PM**(-2)
     & *root3*root5**(-2)*root6
     &  - 16.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*Pp2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  - 12.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &
      K87 = K87 + 4.D0*pp**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*pp**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 12.D0*pp**(-1)*DG*ssp*ssm*Pp2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*PM**(-2)*root2**(-1)
     & *root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &
      K87 = K87 - 4.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 16.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root3*root5**(-2)*
     & root6
     &  + 16.D0*DG*svp*svm*Pq2*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)
     & *root6
     &  - 4.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  + 12.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 12.D0*DG*svp*svm*w1*w2*mn**4*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 4.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 12.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*
     & root6
     &
      K87 = K87 + 12.D0*DG*ssp*ssm*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 4.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*
     & root6
     &  + 4.D0*qq*DG*svp*svm*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)
     & *DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*
     & DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K87 = K87 + 20.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm
     & *w2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 20.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*z**2
     & *f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*z**2
     & *f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*
     & ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*
     & Pq2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &
      K87 = K87 + 16.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*Pq2*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & w2*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 24.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &
      K87 = K87 + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 16.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*f2**(-2)
     & *PM**(-2)*root3*root5**(-2)*root6
     &  + 16.D0*k_p*k_PC*p_q*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K87 = K87 - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*z**2
     & *f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K87 = K87 + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & mn**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2
     & *Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K87 = K87 + 16.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 16.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*z**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*f2**(-2)
     & *PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*z**2*f2**(-2)
     & *PM**(-2)*root3*root5**(-2)*root6
     &
      K87 = K87 + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*f2**(-2)*
     & PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 16.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 16.D0*k_PC**2*p_PC*p_q*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*
     & svp*svm*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &
      K87 = K87 - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2
     & *f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 12.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*svp
     & *svm*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*
     & root6**2
     &
      K87 = K87 + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*
     & svm*w2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*
     & Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &
      K87 = K87 + 12.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*
     & Pp2*z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & z**2*f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*f2**(-2)
     & *PM**(-2)*root3*root5**(-2)*root6
     &  + 12.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*z**2*
     & f2**(-2)*PM**(-2)*root3*root5**(-2)*root6
     &  - 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &
      K87 = K87 + 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*Pp2*f2**(-2)*PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  - 12.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-2)*PM**(-2)*root3
     & *root5**(-2)*root6
     &  + 12.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*z**2*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  - 4.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &  + 12.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*PM**(-2)*root3
     & *root5**(-2)*root6
     &  - 12.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*z**2*f2**(-2)*PM**(-2)*
     & root3*root5**(-2)*root6
     &  + 4.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*
     & root2**(-1)*root5**(-2)*root6**2
     &
      K87 = K87 - 4.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*
     & PM**(-2)*root2**(-1)*root5**(-2)*root6**2
     &  + 4.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  - 12.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  + 12.D0*PC_q*DG*svp*svm*w2*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 4.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*root2**(-1)*
     & root5**(-2)*root6**2
     &  + 12.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &  - 12.D0*PC_q*DG*svp*svm*w1*mn**2*z**2*f2**(-2)*PM**(-2)*root3*
     & root5**(-2)*root6
     &

        elseif((alpha.eq.8).and.(beta.eq.8))then

      K88 = + 4.D0*pp**(-1)*PC2**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*Pp2*
     . f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2*Pp2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*Pq2*Pp2*f2**(-2)*
     & PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*pp**(-1)*DG*svp*svm*w1*w2*mn**2*Pp2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  - 4.D0*pp**(-1)*DG*ssp*ssm*Pp2*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &  + 4.D0*pp**(-1)*qq*DG*svp*svm*Pp2*f2**(-2)*PM**(-2)*root5**(-2)
     & *root6**2
     &  + 4.D0*qq**(-1)*DG*svp*svm*w1*w2*mn**2*Pq2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  - 4.D0*qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &
      K88 = K88 + 4.D0*DG*svp*svm*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &  + 4.D0*DG*svp*svm*w1*w2*mn**4*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &  - 4.D0*DG*ssp*ssm*mn**2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*qq*DG*svp*svm*mn**2*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &  - 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)
     & *DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*
     & DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 12.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p*k_PC*p_PC*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &
      K88 = K88 - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*
     & qq**(-1)*DG*ssp*ssm*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*Pq2
     & *f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & w2*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*DG*ssp*ssm*f2**(-2)*
     & PM**(-2)*root5**(-2)*root6**2
     &  - 8.D0*k_p*k_PC*p_PC*kk**(-1)*pp**(-1)*qq*DG*svp*svm*f2**(-2)*
     & PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*
     & Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*
     & Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &
      K88 = K88 + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & mn**2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p*k_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2
     & *f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2
     & *f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*mn**2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*PC_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*mn**2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*svp*svm*w1*w2*mn**2
     & *Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq**(-1)*DG*ssp*ssm*Pq2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &
      K88 = K88 - 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*Pq2*f2**(-2)
     & *PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*mn**4*f2**(-2)
     & *PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_p**2*kk**(-1)*pp**(-1)*DG*ssp*ssm*mn**2*f2**(-2)*
     & PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_p**2*kk**(-1)*pp**(-1)*qq*DG*svp*svm*mn**2*f2**(-2)*
     & PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG
     & *svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*p_PC*p_q*kk**(-1)*pp**(-1)*DG*svp*svm*w1*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &
      K88 = K88 + 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-2)*
     & qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &  - 8.D0*k_PC**2*PC_q*kk**(-1)*pp**(-1)*PC2**(-1)*DG*svp*svm*w2*
     & Pp2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*qq**(-1)*DG*ssp*ssm*
     & Pq2*Pp2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-2)*DG*svp*svm*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*
     & w1*w2*Pq2*Pp2*f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*DG*ssp*ssm*Pp2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*k_PC**2*kk**(-1)*pp**(-1)*PC2**(-1)*qq*DG*svp*svm*Pp2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &
      K88 = K88 - 4.D0*k_PC**2*kk**(-1)*pp**(-1)*DG*svp*svm*w1*w2*Pp2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  + 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w2*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*PC_q*pp**(-1)*PC2**(-1)*qq**(-1)*DG*svp*svm*w1*Pq2*Pp2*
     & f2**(-2)*PM**(-2)*root5**(-2)*root6**2
     &  - 4.D0*PC_q*pp**(-1)*DG*svp*svm*w2*Pp2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  + 4.D0*PC_q*pp**(-1)*DG*svp*svm*w1*Pp2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  - 4.D0*PC_q*qq**(-1)*DG*svp*svm*w2*Pq2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  + 4.D0*PC_q*qq**(-1)*DG*svp*svm*w1*Pq2*f2**(-2)*PM**(-2)*
     & root5**(-2)*root6**2
     &  - 4.D0*PC_q*DG*svp*svm*w2*mn**2*f2**(-2)*PM**(-2)*root5**(-2)*
     & root6**2
     &
      K88 = K88 + 4.D0*PC_q*DG*svp*svm*w1*mn**2*f2**(-2)*PM**(-2)
     & *root5**(-2)*root6**2
     &
       
         else
         endif

          if((alpha.eq.1).and.(beta.eq.1))then
          tr =  k11
          elseif((alpha.eq.1).and.(beta.eq.2))then
          tr =  k12
          elseif((alpha.eq.1).and.(beta.eq.3))then
          tr =  k13
          elseif((alpha.eq.1).and.(beta.eq.4))then
          tr =  k14
          elseif((alpha.eq.1).and.(beta.eq.5))then
          tr =  k15
          elseif((alpha.eq.1).and.(beta.eq.6))then
          tr =  k16
          elseif((alpha.eq.1).and.(beta.eq.7))then
          tr =  k17
          elseif((alpha.eq.1).and.(beta.eq.8))then
          tr =  k18
          elseif((alpha.eq.2).and.(beta.eq.1))then
          tr =  k21
          elseif((alpha.eq.2).and.(beta.eq.2))then
          tr =  k22
          elseif((alpha.eq.2).and.(beta.eq.3))then
          tr =  k23
          elseif((alpha.eq.2).and.(beta.eq.4))then
          tr =  k24
          elseif((alpha.eq.2).and.(beta.eq.5))then
          tr =  k25
          elseif((alpha.eq.2).and.(beta.eq.6))then
          tr =  k26
          elseif((alpha.eq.2).and.(beta.eq.7))then
          tr =  k27
          elseif((alpha.eq.2).and.(beta.eq.8))then
          tr =  k28
          elseif((alpha.eq.3).and.(beta.eq.1))then
          tr =  k31
          elseif((alpha.eq.3).and.(beta.eq.2))then
          tr =  k32
          elseif((alpha.eq.3).and.(beta.eq.3))then
          tr =  k33
          elseif((alpha.eq.3).and.(beta.eq.4))then
          tr =  k34
          elseif((alpha.eq.3).and.(beta.eq.5))then
          tr =  k35
          elseif((alpha.eq.3).and.(beta.eq.6))then
          tr =  k36
          elseif((alpha.eq.3).and.(beta.eq.7))then
          tr =  k37
          elseif((alpha.eq.3).and.(beta.eq.8))then
          tr =  k38
          elseif((alpha.eq.4).and.(beta.eq.1))then
          tr =  k41
          elseif((alpha.eq.4).and.(beta.eq.2))then
          tr =  k42
          elseif((alpha.eq.4).and.(beta.eq.3))then
          tr =  k43
          elseif((alpha.eq.4).and.(beta.eq.4))then
          tr =  k44
          elseif((alpha.eq.4).and.(beta.eq.5))then
          tr =  k45
          elseif((alpha.eq.4).and.(beta.eq.6))then
          tr =  k46
          elseif((alpha.eq.4).and.(beta.eq.7))then
          tr =  k47
          elseif((alpha.eq.4).and.(beta.eq.8))then
          tr =  k48
          elseif((alpha.eq.5).and.(beta.eq.1))then
          tr =  k51
          elseif((alpha.eq.5).and.(beta.eq.2))then
          tr =  k52
          elseif((alpha.eq.5).and.(beta.eq.3))then
          tr =  k53
          elseif((alpha.eq.5).and.(beta.eq.4))then
          tr =  k54
          elseif((alpha.eq.5).and.(beta.eq.5))then
          tr =  k55
          elseif((alpha.eq.5).and.(beta.eq.6))then
          tr =  k56
          elseif((alpha.eq.5).and.(beta.eq.7))then
          tr =  k57
          elseif((alpha.eq.5).and.(beta.eq.8))then
          tr =  k58
          elseif((alpha.eq.6).and.(beta.eq.1))then
          tr =  k61
          elseif((alpha.eq.6).and.(beta.eq.2))then
          tr =  k62
          elseif((alpha.eq.6).and.(beta.eq.3))then
          tr =  k63
          elseif((alpha.eq.6).and.(beta.eq.4))then
          tr =  k64
          elseif((alpha.eq.6).and.(beta.eq.5))then
          tr =  k65
          elseif((alpha.eq.6).and.(beta.eq.6))then
          tr =  k66
          elseif((alpha.eq.6).and.(beta.eq.7))then
          tr =  k67
          elseif((alpha.eq.6).and.(beta.eq.8))then
          tr =  k68
          elseif((alpha.eq.7).and.(beta.eq.1))then
          tr =  k71
          elseif((alpha.eq.7).and.(beta.eq.2))then
          tr =  k72
          elseif((alpha.eq.7).and.(beta.eq.3))then
          tr =  k73
          elseif((alpha.eq.7).and.(beta.eq.4))then
          tr =  k74
          elseif((alpha.eq.7).and.(beta.eq.5))then
          tr =  k75
          elseif((alpha.eq.7).and.(beta.eq.6))then
          tr =  k76
          elseif((alpha.eq.7).and.(beta.eq.7))then
          tr =  k77
          elseif((alpha.eq.7).and.(beta.eq.8))then
          tr =  k78
          elseif((alpha.eq.8).and.(beta.eq.1))then
          tr =  k81
          elseif((alpha.eq.8).and.(beta.eq.2))then
          tr =  k82
          elseif((alpha.eq.8).and.(beta.eq.3))then
          tr =  k83
          elseif((alpha.eq.8).and.(beta.eq.4))then
          tr =  k84
          elseif((alpha.eq.8).and.(beta.eq.5))then
          tr =  k85
          elseif((alpha.eq.8).and.(beta.eq.6))then
          tr =  k86
          elseif((alpha.eq.8).and.(beta.eq.7))then
          tr =  k87
          elseif((alpha.eq.8).and.(beta.eq.8))then
          tr =  k88
          else
          endif

           kab(alpha,beta)=tr


       if(print.eqv..true.)then
         write(2,*)'================= Kabt variables==========='
       write(2,*)'kk,pp,DG,g2', kk,pp,DG,g2
       write(2,*)' alpha,beta,gamma,rho,mu,nu',alpha,beta,gamma,
     . rho,mu,nu
        write(2,*)"k(i)"
       write(2,*)(k(i),i=1,4)
       write(2,*) "p(i)"
       write(2,*)(p(i),i=1,4)
       write(2,*) "pc(i)"
       write(2,*)(pc(i),i=1,4)
       write(2,*) "q(i)"
       write(2,*)(q(i),i=1,4)
        write(2,*)',p_pc,PC_q,p_q,k_p,k_pc,k_q',p_pc,PC_q,p_q,
     .  k_p,k_pc,k_q
c        write(2,*)'p11,p22,p23,p32,p33,p44,p1,p2,p3,p4',p11,p22,
c     .  p23,p32,p33,p44,p1,p2,p3,p4
        write(2,*)'kab'
        do 104 i = 1,4
        write(2,*)(kab(i,j),j=1,8)
104     continue
c        write(2,*)'A(i,j)'
c        do 105 i = 1,4
c        write(2,*)(A(i,j),j=1,4)
c       105     continue
c        write(2,*)'pab(i,j)'
c        do 106 i = 1,4
c        write(2,*)(pabi(i,j),j=1,4)
c       106     continue


        write(2,*)'sigmasp,sigmasm,sigmavp,sigmavm,qq,delta',
     .             sigmasp,sigmasm,sigmavp,sigmavm,qq,delta
        write(2,*)'pcp,a1,a2,a3,a4,pc2,tr',
     .   pcp,a1,a2,a3,a4,pc2,tr



       print = .false.
       else
       endif


      return
      end
