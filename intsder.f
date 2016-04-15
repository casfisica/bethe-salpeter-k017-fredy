* demo-fortran.F.
*********************************************************************
* This program has adapted the Cuba1.5 Integrator  to Calculate     *
* Sigma*B in pico barns. The momenta units are in TeV. We use       *
* the Cuba subroutines of integration, which calls Cteq6Pdf.f to     *
* calculate  the parton distribution functions (PDF). This program  *
* asumme only  standard model fields, hence defines eps2_R(0)=0.    *
* Contrary to zprime.f. For install the program                     *
* ./configure                                                       *
* this give a new  makefile. For compile this program is need add   *
* into the  makefile, the line:                                     *
* DEMO_F = $(demo)/demo-fortran.F $(demo)/Cteq6Pdf.f $(demo)/zprime.f*
* after this compile the program with  make,                        *
* any cuestion eduardor@fisica.unam.mx                              *
*********************************************************************
*To change some                            
*parameters see below. 31 october 2010
********************************************************************
* test program for the Cuba library
* last modified 2 Feb 05 th


	subroutine  intsder(flag12,nvertex2,ikernel2,np2,pp2,aain2,
     .   bbin2,ccin2,p22,lambda2,pmin22,z22,sigmai)
	implicit none         
        integer ndim, ncomp, mineval, maxeval, verbose, last
        integer introut 
        double precision epsrel, epsabs, sigmai
        double precision pi
       	parameter (pi = 3.14159265358979323846D0)


        double precision p22,lambda2,pmin22,z22
        integer np2,nvertex2,ikernel2,i
        double precision pp2(128),aain2(128),bbin2(128),ccin2(128)
        logical flag12
      
        double precision p2,lambda,pmin2,z2 
        integer np,nvertex,ikernel
        double precision pp,aain,bbin,ccin
        logical flag1

        common /sde1/pp(128),aain(128),bbin(128),ccin(128)
        common /sde2/p2,lambda,pmin2,z2
        common /sde3/np,nvertex,ikernel
        common /sde4/flag1
             
	parameter (ndim = 2)
	parameter (ncomp = 1)
	parameter (epsrel = 1D-3)
	parameter (epsabs = 1D-12)
	parameter (verbose = 2)
	parameter (last = 4)
	parameter (mineval = 0)
	parameter (maxeval = 50000)

	integer nstart, nincrease
	parameter (nstart = 1000000)
	parameter (nincrease = 500000)

	integer nnew
	double precision flatness
	parameter (nnew = 100000)
	parameter (flatness = 1D0)

	integer key1, key2, key3, maxpass
	double precision border, maxchisq, mindeviation
	integer ngiven, ldxgiven, nextra
	parameter (key1 = 47)
	parameter (key2 = 1)
	parameter (key3 = 1)
	parameter (maxpass = 5)
	parameter (border = 0D0)
	parameter (maxchisq = 10D0)
	parameter (mindeviation = .25D0)
	parameter (ngiven = 0)
	parameter (ldxgiven = ndim)
	parameter (nextra = 0)

	integer key
	parameter (key = 0)

	external integrand2

	double precision integral(ncomp), error(ncomp), prob(ncomp)
	integer nregions, neval, fail

	integer c
**********************************
         p2      = p22
         lambda  = lambda2
         pmin2   = pmin22
         z2      = z22
         np      = np2
         nvertex = nvertex2
         ikernel = ikernel2
          
         flag1 = flag12
         
        do 100 i=1,np2
        pp(i)   =pp2(i)
        aain(i) =aain2(i)
        bbin(i)=bbin2(i) 
        ccin(i)=ccin2(i)
100     continue

               
**************************************************
*       Introut = integration soubroutine        *
*       1 vegas, 2 suave, 3 divonne, 4 Cuhre     *
**************************************************

        introut = 4


C	 "------Integration subroutines--------------------"

           

 



        if(introut.eq.1) then
           call vegas(ndim, ncomp, integrand2,
     &    epsrel, epsabs, verbose, mineval, maxeval,
     &    nstart, nincrease,
     &    neval, fail, integral, error, prob)

        elseif(introut.eq.2) then
          call suave(ndim, ncomp, integrand2,
     &    epsrel, epsabs, verbose + last, mineval, maxeval,
     &    nnew, flatness,
     &    nregions, neval, fail, integral, error, prob)
        
        elseif(introut.eq.3) then
          call divonne(ndim, ncomp, integrand2,
     &    epsrel, epsabs, verbose, mineval, maxeval,
     &    key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    nregions, neval, fail, integral, error, prob)
        
         elseif(introut.eq.4) then
        call cuhre(ndim, ncomp, integrand2,
     &    epsrel, epsabs, verbose + last, mineval, maxeval,
     &    key,
     &    nregions, neval, fail, integral, error, prob)
        endif
          
        sigmai =  integral(1) 


	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &  (integral(c), error(c), prob(c), c = 1, ncomp)

              


 
    
       return

	end


************************************************************************

	subroutine integrand2(ndim, xx, ncomp, ff)
	implicit none
        integer ndim, ncomp 
	double precision xx(*), ff(*)
        double precision pi,z1,jac,cf,
     .  k2,q2,pk,k,tr


        double precision p2,lambda,pmin2,z2
        integer np,nvertex,ikernel
        double precision pp(128),aain(128),bbin(128),ccin(128)
        logical flag1

        common /sde1/pp,aain,bbin,ccin
        common /sde2/p2,lambda,pmin2,z2
        common /sde3/np,nvertex,ikernel
        common /sde4/flag1



        external massf
        external vertexr


        pi = 3.14159265358979323846D0
        

             


   


* aqui introduzco la llamada de pdf
*=============================================================================
*     ESTE PROGRAMA DEFINE EL INTEGRANDO PARA CUBA
*=============================================================================
*


       
       k2 = (lambda**2-pmin2)*xx(1)+pmin2
       z1 = 2*xx(2)-1d0

       jac =2*pi*2*(lambda**2-pmin2)*k2*dsqrt(1d0-z1**2) 
       cf = 4d0/3d0
       k = abs(dsqrt(k2))
       

       pk = dsqrt(p2*k2)*z1
       q2 = p2+k2-2*pk


     
        
        call vertexr (flag1,nvertex,ikernel,np,pp,aain,
     .            bbin,ccin,p2,lambda,pmin2,z2,k2,z1,tr)

      
  


        ff(1)  = tr/(2*pi)**4*jac

       
        
	end




