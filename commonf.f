       double precision xp,wp,xv,wv,xt,wt,xs,ws,lambda,pmin2,
     . eplus,eminus,mh,mh2,xval,xval2
       integer np,ntr,npol,npol2,nt,ns
       double complex aa,bb, aa2,bb2
       logical ffit

       common /f1/xp(128),Wp(128), xv(128), Wv(128), xt(128),wt(128),
     . xs(128),ws(128),lambda,pmin2,eplus,eminus,mh,mh2,
     . xval(26),xval2(26)
       common /f2/np,ntr,npol,npol2,nt,ns
       common /f3/AA(128,3),BB(128,3),AA2(128,3),BB2(128,3)
       common /f4/ffit
