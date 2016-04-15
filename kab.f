       subroutine kabT (print,pc,p,q,k,sigmavp,sigmasp,sigmavm,
     .  sigmasm,g2,alpha,beta,tr)
       implicit none
c       include 'common.f'

       double precision  kk,pp,DG,g2
       Integer alpha,beta,gamma,rho,mu,nu
       double precision  k(4),p(4),q(4) 
       double complex    PC(4),p_pc,PC_q,p_q,k_p,k_pc,k_q,
     . p11,p22,p23,p32,p33,p44,p1,p2,p3,p4,kab(4,4)  
       double complex  sigmasp,sigmasm,sigmavp,sigmavm,qq,delta 
       double complex A(4,4),pcp,pabi(4,4),a1,a2,a3,a4,pc2,tr,k11,
     . k12,k13,k14,k21,k22,k23,k24,k31,k32,k33,k34,k41,k42,k43,k44
       integer i,j
       logical print

      



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


       do 102 i = 1,4
       pc2 =pc(i)*pc(i)+pc2
       pp=p(i)*p(i)    +pp
       qq=q(i)*q(i)   +qq
       kk=k(i)*k(i)    +kk
       p_pc=p(i)*pc(i)+p_pc
       pc_q=pc(i)*q(i)+pc_q
       p_q =p(i)*q(i) +p_q
       k_p =k(i)*p(i) +k_p
       k_pc=k(i)*pc(i)+k_pc
       k_q =k(i)*q(i) +k_q

       do 103 j = 1,4
       A(i,j) = 0d0
       kab(i,j)=0d0
103    continue
102    continue
   

        dg = -4d0/3d0*g2/kk   

   
        pcp =p_pc !dsqrt(pp)*pc(4)*zp
c       p_pc=pcp

       A(1,1) = -1d0
       A(2,2) = -PC(4)**2
       A(2,3) =  -pcp**2
       A(3,2) =  -pcp**2
       A(3,3) = -pp*pcp**2
       A(4,4) =  pp*pc(4)**2-(pcp)**2 

       delta = A(3,3)*A(2,2)-A(3,2)*A(2,3)+1d-10
                

       do 100 i = 1,4
       do 101 j = 1,4
       Pabi(i,j) = 0d0
101    continue
100    continue
       pabi(1,1) = 1d0/A(1,1)
       pabi(2,2) =  A(3,3)/delta
       pabi(2,3) = -A(2,3)/delta
       pabi(3,2) = -A(3,2)/delta
       pabi(3,3) =  A(2,2)/delta
       pabi(4,4) =  1d0/A(4,4)



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

                K11 = - 3.D0*DG*sigmasp*sigmasm
     &  - 3.D0*qq*DG*sigmavp*sigmavm
     &  + 3.D0/4.D0*PC2*DG*sigmavp*sigmavm
     &

          tr= k11
        elseif((alpha.eq.1).and.(beta.eq.2))then

      K12 = - 3.D0/2.D0*PC2*DG*sigmasm*sigmavp
     &  - 3.D0/2.D0*PC2*DG*sigmasp*sigmavm
     &  - 3.D0*PC_q*DG*sigmasm*sigmavp
     &  + 3.D0*PC_q*DG*sigmasp*sigmavm
     &

           tr= k12
         elseif((alpha.eq.1).and.(beta.eq.3))then

      K13 = - 3.D0*PC_q*qq*DG*sigmasm*sigmavp
     &  + 3.D0*PC_q*qq*DG*sigmasp*sigmavm
     &  - 3.D0/2.D0*PC_q**2*DG*sigmasm*sigmavp
     &  - 3.D0/2.D0*PC_q**2*DG*sigmasp*sigmavm
     &

         tr= k13
         elseif((alpha.eq.1).and.(beta.eq.4))then

      K14 = + 3.D0*PC2*qq*DG*sigmavp*sigmavm
     &  - 3.D0*PC_q**2*DG*sigmavp*sigmavm
     &


             tr= k14

          elseif((alpha.eq.2).and.(beta.eq.1))then

      K21 = + k_p*k_PC*p_PC**3*kk**(-1)*delta**(-1)*DG*sigmasm*sigmavp
     &  + k_p*k_PC*p_PC**3*kk**(-1)*delta**(-1)*DG*sigmasp*sigmavm
     &  + 2.D0*k_p*k_q*p_PC**3*kk**(-1)*delta**(-1)*DG*sigmasm*sigmavp
     &  - 2.D0*k_p*k_q*p_PC**3*kk**(-1)*delta**(-1)*DG*sigmasp*sigmavm
     &  - 2.D0*k_PC*k_q*p_PC**2*kk**(-1)*pp*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  + 2.D0*k_PC*k_q*p_PC**2*kk**(-1)*pp*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - k_PC**2*p_PC**2*kk**(-1)*pp*delta**(-1)*DG*sigmasm*sigmavp
     &  - k_PC**2*p_PC**2*kk**(-1)*pp*delta**(-1)*DG*sigmasp*sigmavm
     &  - p_PC**2*PC_q*pp*delta**(-1)*DG*sigmasm*sigmavp
     &  + p_PC**2*PC_q*pp*delta**(-1)*DG*sigmasp*sigmavm
     &  - 1.D0/2.D0*p_PC**2*pp*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  - 1.D0/2.D0*p_PC**2*pp*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + p_PC**3*p_q*delta**(-1)*DG*sigmasm*sigmavp
     &
      K21 = K21 - p_PC**3*p_q*delta**(-1)*DG*sigmasp*sigmavm
     &  + 1.D0/2.D0*p_PC**4*delta**(-1)*DG*sigmasm*sigmavp
     &  + 1.D0/2.D0*p_PC**4*delta**(-1)*DG*sigmasp*sigmavm
     &



          tr= k21

         elseif((alpha.eq.2).and.(beta.eq.2))then  




      K22 = - 2.D0*k_p*k_PC*p_PC**3*kk**(-1)*delta**(-1)*DG*sigmasp*
     . sigmasm
     &  + 2.D0*k_p*k_PC*p_PC**3*kk**(-1)*qq*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  + 1.D0/2.D0*k_p*k_PC*p_PC**3*kk**(-1)*PC2*delta**(-1)*DG*
     & sigmavp*sigmavm
     &  - 4.D0*k_p*k_q*p_PC**3*PC_q*kk**(-1)*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  + 4.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*pp*delta**(-1)*DG*sigmavp
     & *sigmavm
     &  + 2.D0*k_PC**2*p_PC**2*kk**(-1)*pp*delta**(-1)*DG*sigmasp*
     & sigmasm
     &  - 2.D0*k_PC**2*p_PC**2*kk**(-1)*pp*qq*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  - 1.D0/2.D0*k_PC**2*p_PC**2*kk**(-1)*pp*PC2*delta**(-1)*DG*
     & sigmavp*sigmavm
     &
      K22 = K22 + 2.D0*p_PC**2*PC_q**2*pp*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  + p_PC**2*pp*PC2*delta**(-1)*DG*sigmasp*sigmasm
     &  - p_PC**2*pp*PC2*qq*delta**(-1)*DG*sigmavp*sigmavm
     &  - 1.D0/4.D0*p_PC**2*pp*PC2**2*delta**(-1)*DG*sigmavp*sigmavm
     &  - 2.D0*p_PC**3*p_q*PC_q*delta**(-1)*DG*sigmavp*sigmavm
     &  - p_PC**4*delta**(-1)*DG*sigmasp*sigmasm
     &  + p_PC**4*qq*delta**(-1)*DG*sigmavp*sigmavm
     &  + 1.D0/4.D0*p_PC**4*PC2*delta**(-1)*DG*sigmavp*sigmavm
     &


        tr= k22
       elseif((alpha.eq.2).and.(beta.eq.3))then

      K23 = + k_p*k_PC*p_PC**3*PC_q**2*kk**(-1)*delta**(-1)*DG*sigmavp*
     . sigmavm
     &  - 2.D0*k_p*k_q*p_PC**3*PC_q*kk**(-1)*delta**(-1)*DG*sigmasp*
     & sigmasm
     &  - 2.D0*k_p*k_q*p_PC**3*PC_q*kk**(-1)*qq*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  - 1.D0/2.D0*k_p*k_q*p_PC**3*PC_q*kk**(-1)*PC2*delta**(-1)*DG*
     & sigmavp*sigmavm
     &  + 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*pp*delta**(-1)*DG*sigmasp
     & *sigmasm
     &  + 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*pp*qq*delta**(-1)*DG*
     & sigmavp*sigmavm
     &  + 1.D0/2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*pp*PC2*delta**(-1)*
     & DG*sigmavp*sigmavm
     &  - k_PC**2*p_PC**2*PC_q**2*kk**(-1)*pp*delta**(-1)*DG*sigmavp*
     & sigmavm
     &
      K23 = K23 + p_PC**2*PC_q**2*pp*delta**(-1)*DG*sigmasp*sigmasm
     &  + p_PC**2*PC_q**2*pp*qq*delta**(-1)*DG*sigmavp*sigmavm
     &  - 1.D0/4.D0*p_PC**2*PC_q**2*pp*PC2*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  - p_PC**3*p_q*PC_q*delta**(-1)*DG*sigmasp*sigmasm
     &  - p_PC**3*p_q*PC_q*qq*delta**(-1)*DG*sigmavp*sigmavm
     &  - 1.D0/4.D0*p_PC**3*p_q*PC_q*PC2*delta**(-1)*DG*sigmavp*sigmavm
     &  + 1.D0/2.D0*p_PC**4*PC_q**2*delta**(-1)*DG*sigmavp*sigmavm
     &

           tr= k23
         elseif((alpha.eq.2).and.(beta.eq.4))then


      K24 = + k_p*k_PC*p_PC**3*PC_q*kk**(-1)*delta**(-1)*DG*sigmasm*
     . sigmavp
     &  - k_p*k_PC*p_PC**3*PC_q*kk**(-1)*delta**(-1)*DG*sigmasp*sigmavm
     &  + 2.D0*k_p*k_PC*p_PC**3*kk**(-1)*qq*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  + 2.D0*k_p*k_PC*p_PC**3*kk**(-1)*qq*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - 2.D0*k_p*k_q*p_PC**3*PC_q*kk**(-1)*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  - 2.D0*k_p*k_q*p_PC**3*PC_q*kk**(-1)*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - k_p*k_q*p_PC**3*kk**(-1)*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  + k_p*k_q*p_PC**3*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*pp*delta**(-1)*DG*sigmasm
     & *sigmavp
     &
      K24 = K24 + 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*pp*delta**(-1)*DG
     & *sigmasp*sigmavm
     &  + k_PC*k_q*p_PC**2*kk**(-1)*pp*PC2*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  - k_PC*k_q*p_PC**2*kk**(-1)*pp*PC2*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - k_PC**2*p_PC**2*PC_q*kk**(-1)*pp*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  + k_PC**2*p_PC**2*PC_q*kk**(-1)*pp*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - 2.D0*k_PC**2*p_PC**2*kk**(-1)*pp*qq*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  - 2.D0*k_PC**2*p_PC**2*kk**(-1)*pp*qq*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  + p_PC**2*PC_q**2*pp*delta**(-1)*DG*sigmasm*sigmavp
     &
      K24 = K24 + p_PC**2*PC_q**2*pp*delta**(-1)*DG*sigmasp*sigmavm
     &  - p_PC**2*pp*PC2*qq*delta**(-1)*DG*sigmasm*sigmavp
     &  - p_PC**2*pp*PC2*qq*delta**(-1)*DG*sigmasp*sigmavm
     &  - p_PC**3*p_q*PC_q*delta**(-1)*DG*sigmasm*sigmavp
     &  - p_PC**3*p_q*PC_q*delta**(-1)*DG*sigmasp*sigmavm
     &  - 1.D0/2.D0*p_PC**3*p_q*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  + 1.D0/2.D0*p_PC**3*p_q*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + 1.D0/2.D0*p_PC**4*PC_q*delta**(-1)*DG*sigmasm*sigmavp
     &  - 1.D0/2.D0*p_PC**4*PC_q*delta**(-1)*DG*sigmasp*sigmavm
     &  + p_PC**4*qq*delta**(-1)*DG*sigmasm*sigmavp
     &  + p_PC**4*qq*delta**(-1)*DG*sigmasp*sigmavm
     &

             tr= k24
         elseif((alpha.eq.3).and.(beta.eq.1))then


      K31 = - k_p*k_PC*p_PC*kk**(-1)*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  - k_p*k_PC*p_PC*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  - 2.D0*k_p*k_q*p_PC*kk**(-1)*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  + 2.D0*k_p*k_q*p_PC*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + 2.D0*k_PC*k_q*p_PC**2*kk**(-1)*delta**(-1)*DG*sigmasm*sigmavp
     &  - 2.D0*k_PC*k_q*p_PC**2*kk**(-1)*delta**(-1)*DG*sigmasp*sigmavm
     &  + k_PC**2*p_PC**2*kk**(-1)*delta**(-1)*DG*sigmasm*sigmavp
     &  + k_PC**2*p_PC**2*kk**(-1)*delta**(-1)*DG*sigmasp*sigmavm
     &  - p_PC*p_q*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  + p_PC*p_q*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + p_PC**2*PC_q*delta**(-1)*DG*sigmasm*sigmavp
     &  - p_PC**2*PC_q*delta**(-1)*DG*sigmasp*sigmavm
     &

                  tr= k31
          elseif((alpha.eq.3).and.(beta.eq.2))then

      K32 = + 2.D0*k_p*k_PC*p_PC*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*
     . sigmasm
     &  - 2.D0*k_p*k_PC*p_PC*kk**(-1)*PC2*qq*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  - 1.D0/2.D0*k_p*k_PC*p_PC*kk**(-1)*PC2**2*delta**(-1)*DG*
     & sigmavp*sigmavm
     &  + 4.D0*k_p*k_q*p_PC*PC_q*kk**(-1)*PC2*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  - 4.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  - 2.D0*k_PC**2*p_PC**2*kk**(-1)*delta**(-1)*DG*sigmasp*sigmasm
     &  + 2.D0*k_PC**2*p_PC**2*kk**(-1)*qq*delta**(-1)*DG*sigmavp*
     & sigmavm
     &  + 1.D0/2.D0*k_PC**2*p_PC**2*kk**(-1)*PC2*delta**(-1)*DG*sigmavp
     & *sigmavm
     &
      K32 = K32 + 2.D0*p_PC*p_q*PC_q*PC2*delta**(-1)*DG*sigmavp*sigmavm
     &  - 2.D0*p_PC**2*PC_q**2*delta**(-1)*DG*sigmavp*sigmavm
     &



             tr= k32
         elseif((alpha.eq.3).and.(beta.eq.3))then

      K33 = - k_p*k_PC*p_PC*PC_q**2*kk**(-1)*PC2*delta**(-1)*DG*sigmavp
     . *sigmavm
     &  + 2.D0*k_p*k_q*p_PC*PC_q*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*
     & sigmasm
     &  + 2.D0*k_p*k_q*p_PC*PC_q*kk**(-1)*PC2*qq*delta**(-1)*DG*sigmavp
     & *sigmavm
     &  + 1.D0/2.D0*k_p*k_q*p_PC*PC_q*kk**(-1)*PC2**2*delta**(-1)*DG*
     & sigmavp*sigmavm
     &  - 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*delta**(-1)*DG*sigmasp*
     & sigmasm
     &  - 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*qq*delta**(-1)*DG*sigmavp
     & *sigmavm
     &  - 1.D0/2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*PC2*delta**(-1)*DG*
     & sigmavp*sigmavm
     &  + k_PC**2*p_PC**2*PC_q**2*kk**(-1)*delta**(-1)*DG*sigmavp*
     & sigmavm
     &
      K33 = K33 + p_PC*p_q*PC_q*PC2*delta**(-1)*DG*sigmasp*sigmasm
     &  + p_PC*p_q*PC_q*PC2*qq*delta**(-1)*DG*sigmavp*sigmavm
     &  + 1.D0/4.D0*p_PC*p_q*PC_q*PC2**2*delta**(-1)*DG*sigmavp*sigmavm
     &  - p_PC**2*PC_q**2*delta**(-1)*DG*sigmasp*sigmasm
     &  - p_PC**2*PC_q**2*qq*delta**(-1)*DG*sigmavp*sigmavm
     &  - 1.D0/4.D0*p_PC**2*PC_q**2*PC2*delta**(-1)*DG*sigmavp*sigmavm
     &

                   tr= k33
               elseif((alpha.eq.3).and.(beta.eq.4))then

      K34 = - k_p*k_PC*p_PC*PC_q*kk**(-1)*PC2*delta**(-1)*DG*sigmasm*
     . sigmavp
     &  + k_p*k_PC*p_PC*PC_q*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - 2.D0*k_p*k_PC*p_PC*kk**(-1)*PC2*qq*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  - 2.D0*k_p*k_PC*p_PC*kk**(-1)*PC2*qq*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  + 2.D0*k_p*k_q*p_PC*PC_q*kk**(-1)*PC2*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  + 2.D0*k_p*k_q*p_PC*PC_q*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  + k_p*k_q*p_PC*kk**(-1)*PC2**2*delta**(-1)*DG*sigmasm*sigmavp
     &  - k_p*k_q*p_PC*kk**(-1)*PC2**2*delta**(-1)*DG*sigmasp*sigmavm
     &  - 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*delta**(-1)*DG*sigmasm*
     & sigmavp
     &
      K34 = K34 - 2.D0*k_PC*k_q*p_PC**2*PC_q*kk**(-1)*delta**(-1)*DG*
     & sigmasp*sigmavm
     &  - k_PC*k_q*p_PC**2*kk**(-1)*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  + k_PC*k_q*p_PC**2*kk**(-1)*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + k_PC**2*p_PC**2*PC_q*kk**(-1)*delta**(-1)*DG*sigmasm*sigmavp
     &  - k_PC**2*p_PC**2*PC_q*kk**(-1)*delta**(-1)*DG*sigmasp*sigmavm
     &  + 2.D0*k_PC**2*p_PC**2*kk**(-1)*qq*delta**(-1)*DG*sigmasm*
     & sigmavp
     &  + 2.D0*k_PC**2*p_PC**2*kk**(-1)*qq*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  + p_PC*p_q*PC_q*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &  + p_PC*p_q*PC_q*PC2*delta**(-1)*DG*sigmasp*sigmavm
     &  + 1.D0/2.D0*p_PC*p_q*PC2**2*delta**(-1)*DG*sigmasm*sigmavp
     &  - 1.D0/2.D0*p_PC*p_q*PC2**2*delta**(-1)*DG*sigmasp*sigmavm
     &  - 1.D0/2.D0*p_PC**2*PC_q*PC2*delta**(-1)*DG*sigmasm*sigmavp
     &
      K34 = K34 + 1.D0/2.D0*p_PC**2*PC_q*PC2*delta**(-1)*DG*sigmasp*
     & sigmavm
     &  - p_PC**2*PC_q**2*delta**(-1)*DG*sigmasm*sigmavp
     &  - p_PC**2*PC_q**2*delta**(-1)*DG*sigmasp*sigmavm
     &

            tr= k34
           elseif((alpha.eq.4).and.(beta.eq.1))then

      K41 = - 2.D0/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q*kk**(-1)*DG*sigmavp
     . *sigmavm
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_p*k_q*kk**(-1)*PC2*DG*sigmavp*
     & sigmavm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*kk**(-1)*DG*sigmavp*
     & sigmavm
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_PC**2*p_q*kk**(-1)*DG*sigmavp*
     & sigmavm
     &  + 1/(pp*PC2 - p_PC**2)*p_PC*PC_q*DG*sigmavp*sigmavm
     &  - 1/(pp*PC2 - p_PC**2)*p_q*PC2*DG*sigmavp*sigmavm
     &


            tr= k41
          elseif((alpha.eq.4).and.(beta.eq.2))then

      K42 = + 2.D0/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q*kk**(-1)*DG*sigmasm
     . *sigmavp
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q*kk**(-1)*DG*sigmasp*
     & sigmavm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_p*k_q*kk**(-1)*PC2*DG*sigmasm*
     & sigmavp
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_p*k_q*kk**(-1)*PC2*DG*sigmasp*
     & sigmavm
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*kk**(-1)*DG*sigmasm*
     & sigmavp
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*kk**(-1)*DG*sigmasp*
     & sigmavm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_PC**2*p_q*kk**(-1)*DG*sigmasm*
     & sigmavp
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_PC**2*p_q*kk**(-1)*DG*sigmasp*
     & sigmavm
     &
      K42 = K42 - 1/(pp*PC2 - p_PC**2)*p_PC*PC_q*DG*sigmasm*sigmavp
     &  - 1/(pp*PC2 - p_PC**2)*p_PC*PC_q*DG*sigmasp*sigmavm
     &  + 1/(pp*PC2 - p_PC**2)*p_q*PC2*DG*sigmasm*sigmavp
     &  + 1/(pp*PC2 - p_PC**2)*p_q*PC2*DG*sigmasp*sigmavm
     &

            tr= k42

           elseif((alpha.eq.4).and.(beta.eq.3))then

      K43 = - 1/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q**2*kk**(-1)*DG*sigmasm
     .  *sigmavp
     &  + 1/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q**2*kk**(-1)*DG*sigmasp*
     & sigmavm
     &  + 1/(pp*PC2 - p_PC**2)*k_p*k_q*PC_q*kk**(-1)*PC2*DG*sigmasm*
     & sigmavp
     &  - 1/(pp*PC2 - p_PC**2)*k_p*k_q*PC_q*kk**(-1)*PC2*DG*sigmasp*
     & sigmavm
     &  - 1/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*PC_q*kk**(-1)*DG*sigmasm*
     & sigmavp
     &  + 1/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*PC_q*kk**(-1)*DG*sigmasp*
     & sigmavm
     &  + 1/(pp*PC2 - p_PC**2)*k_PC**2*p_q*PC_q*kk**(-1)*DG*sigmasm*
     & sigmavp
     &  - 1/(pp*PC2 - p_PC**2)*k_PC**2*p_q*PC_q*kk**(-1)*DG*sigmasp*
     & sigmavm
     &
      K43 = K43 + 1.D0/2.D0/(pp*PC2 - p_PC**2)*p_PC*PC_q**2*DG*sigmasm*
     & sigmavp
     &  - 1.D0/2.D0/(pp*PC2 - p_PC**2)*p_PC*PC_q**2*DG*sigmasp*sigmavm
     &  - 1.D0/2.D0/(pp*PC2 - p_PC**2)*p_q*PC_q*PC2*DG*sigmasm*sigmavp
     &  + 1.D0/2.D0/(pp*PC2 - p_PC**2)*p_q*PC_q*PC2*DG*sigmasp*sigmavm
     &
             

          tr= k43

          elseif((alpha.eq.4).and.(beta.eq.4))then

      K44 = + 2.D0/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q*kk**(-1)*DG*sigmasp
     .  *sigmasm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q*kk**(-1)*qq*DG*sigmavp*
     & sigmavm
     &  + 1.D0/2.D0/(pp*PC2 - p_PC**2)*k_p*k_PC*PC_q*kk**(-1)*PC2*DG*
     & sigmavp*sigmavm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_p*k_q*kk**(-1)*PC2*DG*sigmasp*
     & sigmasm
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_p*k_q*kk**(-1)*PC2*qq*DG*sigmavp*
     & sigmavm
     &  - 1.D0/2.D0/(pp*PC2 - p_PC**2)*k_p*k_q*kk**(-1)*PC2**2*DG*
     & sigmavp*sigmavm
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*kk**(-1)*DG*sigmasp*
     & sigmasm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*kk**(-1)*qq*DG*sigmavp*
     & sigmavm
     &
      K44 = K44 + 1.D0/2.D0/(pp*PC2 - p_PC**2)*k_PC*k_q*p_PC*kk**(-1)*
     & PC2*DG*sigmavp*sigmavm
     &  - 2.D0/(pp*PC2 - p_PC**2)*k_PC**2*p_q*kk**(-1)*DG*sigmasp*
     & sigmasm
     &  + 2.D0/(pp*PC2 - p_PC**2)*k_PC**2*p_q*kk**(-1)*qq*DG*sigmavp*
     & sigmavm
     &  - 1.D0/2.D0/(pp*PC2 - p_PC**2)*k_PC**2*p_q*kk**(-1)*PC2*DG*
     & sigmavp*sigmavm
     &  - 1/(pp*PC2 - p_PC**2)*p_PC*PC_q*DG*sigmasp*sigmasm
     &  + 1/(pp*PC2 - p_PC**2)*p_PC*PC_q*qq*DG*sigmavp*sigmavm
     &  - 1.D0/4.D0/(pp*PC2 - p_PC**2)*p_PC*PC_q*PC2*DG*sigmavp*sigmavm
     &  + 1/(pp*PC2 - p_PC**2)*p_q*PC2*DG*sigmasp*sigmasm
     &  - 1/(pp*PC2 - p_PC**2)*p_q*PC2*qq*DG*sigmavp*sigmavm
     &  + 1.D0/4.D0/(pp*PC2 - p_PC**2)*p_q*PC2**2*DG*sigmavp*sigmavm
     &


             tr =k44           

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
        write(2,*)'p11,p22,p23,p32,p33,p44,p1,p2,p3,p4',p11,p22,
     .  p23,p32,p33,p44,p1,p2,p3,p4 
        write(2,*)'kab'
        do 104 i = 1,4       
        write(2,*)(kab(i,j),j=1,4)
104     continue
        write(2,*)'A(i,j)'         
        do 105 i = 1,4
        write(2,*)(A(i,j),j=1,4)
105     continue
        write(2,*)'pab(i,j)'
        do 106 i = 1,4
        write(2,*)(pabi(i,j),j=1,4)
106     continue


        write(2,*)'sigmasp,sigmasm,sigmavp,sigmavm,qq,delta',
     .             sigmasp,sigmasm,sigmavp,sigmavm,qq,delta  
        write(2,*)'pcp,a1,a2,a3,a4,pc2,tr',
     .   pcp,a1,a2,a3,a4,pc2,tr
       


       print = .false.
       else
       endif


       return
       end
