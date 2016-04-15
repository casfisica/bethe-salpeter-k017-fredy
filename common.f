
       double precision  p,lambda,pmin2,
     . gamma1,gamma2,gamma3,gamma4,uno,gamma5,gamma0,
     . xp,Wp,xv,Wv,eplus,eminus,mh,mh2,q2E,uhat,xq,wq,
     . xs,ws,xt,wt,kbsr,kbsi,lamb_bs,z2k,xval,xval2

       integer  np,jbv,jb,np2,alpha,beta,m,n,nzq,jp,jq,ns,nt,nm,
     . npol,npol2
       

       logical bmass,print,linear,flag2,init,printkabt,ffit


       double complex cgamma,trb,aa,bb,aa2,bb2,kabR,ssigmavp,ssigmasp,
     . ssigmavm,ssigmasm

       common /block1/ 
     .   p,lambda,pmin2,gamma1(4,4),gamma2(4,4),gamma3(4,4),
     .  gamma4(4,4),uno(4,4),gamma5(4,4),gamma0,xp(128),Wp(128),
     .  xv(128), Wv(128),eplus,eminus,mh,mh2,q2E,uhat(4),
     .  xq(128), Wq(128),xs(128),ws(128),xt(128),wt(128),
     .  kbsr(8,8,0:5,0:5,64,64),kbsi(8,8,0:5,0:5,64,64),lamb_bs,
     .  z2k,xval(26),xval2(26)
     


 
       common /block2/np,jbv(7,2000,4),jb(4),np2,alpha,beta,m,n,
     . nzq,jp,jq,ns,nt,nm,npol,npol2

       common /block3/ bmass,print,linear,flag2,init,printkabt,ffit

       common /block4/ cgamma(4,4,0:5),trb(2000,4),AA(128,3),
     . BB(128,3),AA2(128,3),BB2(128,3),kabR(64,64,30,30,25),
     . ssigmavp(64,30),ssigmasp(64,30),ssigmavm(64,30),ssigmasm(64,30)

