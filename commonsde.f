
       double precision  pp,lambda,pmin2,CC,CCin,X0,ppmil,ccc,
     . ppc,FAK1m,FAK2m,FAK3m,FAK1p,FAK2p,FAK3p,pr,xp,Wp,xv,Wv,
     . uppv,ppv,eplus,eminus,mh,pi,wzq,xzq,z2g
   
        double complex p,AA, BB,AAin,BBin,pc4
   


       integer  np,nvertex,ikernel,np2,npc,contour,nzq
       

       logical print, linear, flag2


     
       common /block8/pp(128),lambda,pmin2,CC(128),CCin(128),X0(73,73),
     .  ppmil(73),ccc(200),ppc(200),FAK1m(200,200),FAK2m(200,200),
     .  FAK3m(200,200),FAK1p(200,200),FAK2p(200,200),FAK3p(200,200),
     .  pr,xp(128),Wp(128),xv(128), Wv(128),uppv(128),ppv(128),eplus,
     .  eminus,mh,pi,wzq(128),xzq(128),z2g



 
       common /block9/np,nvertex,ikernel,np2,npc,contour,nzq

       common /block10/ print, linear,flag2

       common /block11/p,AA(128,3), BB(128,3),AAin(128,3),
     .  BBin(128,3),pc4
 
