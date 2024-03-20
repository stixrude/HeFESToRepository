        subroutine therml(ispec,Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'
	
	logical sfix,vofix,cvfix
	integer ispec,i,j,nbm,nobm,noth,mfit,maxij,lineart,noln
	double precision vi,volnl,cp,cv,gamma,alp,ftot,ph,ent,deltas,tcal,zeta,gsh,uth,uto,thet,q,etas
	double precision pzp,akt,aktel,aktig,aktxs,apar,cvel,cvig,cvxs,d2fdv2,d2tdt2,dgdt
	double precision dfdv,dtdt,e,eel,eig,exs,f,fac,fel
	double precision fig,fmth,fn,fxs,Pi,pig,sel,sig,sxs,theta,Ti,to,vo,dfac
        double precision K,Ks,Kxs,Kig,Kel,Kelp,P,Pxs,Pel
        double precision daktdvel,betael
	double precision aliq(nparp,nparp),aliqc,cliq0,cliq1,cliq2,tee,dteedt,acof,bcof
	double precision betaig,daktdvig,Kigp,d1mach
        common /state/ apar(nspecp,nparp),Ti,Pi
	common /liqc/ aliqc(nparp,nparp),mfit,nobm,noth,sfix,vofix,cvfix,maxij,lineart,nbm,noln
	double precision, parameter :: fsmall=1.e-12
c	fac = hplanck/sqrt(2.*pirad*boltzk)
c	ee = exp(1.)

C  For fitting, get parameters from aliqc
	do 21 i=1,nparp
	 do 21 j=1,nparp
	  aliq(i,j) = aliqc(i,j)
21	continue
C  For forward code, get parameters from aliqset
        call aliqset(ispec,aliq)
	
	Vo = apar(ispec,6)
	To = apar(ispec,4)
	fmth = apar(ispec,33)
	fn = apar(ispec,1)

	f = 0.5*((Vo/Vi)**(2./3.) - 1.)
	if (f .eq. 0.) f = fsmall
	dfdv = -(2.*f + 1.)**(2.5)/(3.*Vo)
	d2fdv2 = 5.*(2.*f + 1)**4/(9.*Vo*Vo)
	theta = ((Ti/To)**fmth - 1.)
	tee = Ti/To - 1.
	dteedt = 1./To
	if (theta .eq. 0.) theta = fsmall
	dtdt = fmth*(theta + 1.)/Ti
	d2tdt2 = fmth*(fmth - 1.)*(theta + 1.)/Ti**2
c	pig = 0.001*fn*Rgas*Ti/Vi

	call thermlel(ispec,Vi,Fel,Eel,Sel,Pel,Cvel,betael,Kel,Kelp,aktel,daktdvel)
	call thermlig(ispec,Vi,Fig,Eig,Sig,Pig,Cvig,betaig,Kig,Kigp,aktig,daktdvig)

	Fxs = 0.
	do 1 i=0,nobm
	 do 1 j=0,noth
	  if (i+j .ge. maxij) go to 1
	  Fxs = Fxs + aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i)*theta**(j)
1	continue
	do 11 i=0,noln
	 Fxs = Fxs + tee*aliq(i+1,lineart)*f**i/dfac(i)
11	continue
	Ftot = 1000.*Fxs + 1000.*Fel + Fig
c	print*, 'in therml',ispec,Fxs,Fel,Fig/1000.

	Sxs = 0.
	do 2 i=0,nobm
	 do 2 j=0,noth
	  if (i+j .ge. maxij) go to 2
	  Sxs = Sxs + float(j)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i)*theta**(j-1)
2	continue
	Sxs = -dtdt*Sxs
	do 12 i=0,noln
	 Sxs = Sxs - dteedt*aliq(i+1,lineart)*f**i/dfac(i)
12	continue
	ent = 1000.*Sxs + 1000.*Sel + Sig
	write(31,*) 'Entropy',Ti,ent,1000.*Sxs+Sig

        Pxs = 0.
        do 7 i=1,nobm
         do 7 j=0,noth
          if (i+j .ge. maxij) go to 7
          Pxs = Pxs + float(i)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i-1)*theta**(j)
7       continue
	do 17 i=1,noln
	 Pxs = Pxs + tee*aliq(i+1,lineart)*float(i)*f**(i-1)/dfac(i)
17	continue
        Pxs = -dfdv*Pxs
	P = Pxs + Pig + Pel

	aktxs = 0.
	do 5 i=0,nobm
	 do 5 j=0,noth
	  if (i+j .ge. maxij) go to 5
	  aktxs = aktxs + float(i)*float(j)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i-1)*theta**(j-1)
5	continue
	aktxs = -dfdv*dtdt*aktxs
	do 15 i=1,noln
	 aktxs = aktxs - dfdv*dteedt*aliq(i+1,lineart)*float(i)*f**(i-1)/dfac(i)
15	continue
	akt = aktxs + aktig + aktel

	Kxs = 0.
	do 3 i=0,nobm
	 do 3 j=0,noth
	  if (i+j .ge. maxij) go to 3
	  Kxs = Kxs + float(i)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*theta**(j)*(d2fdv2*f**(i-1) + dfdv**2*float(i-1)*f**(i-2))
3	continue
	do 13 i=1,noln
	 Kxs = Kxs + d2fdv2*tee*aliq(i+1,lineart)*float(i)*f**(i-1)/dfac(i)
13	continue
	Kxs = Vi*Kxs
	K = Kxs + Kig + Kel

	Cvxs = 0.
	do 4 i=0,nobm
	 do 4 j=0,noth
	  if (i+j .ge. maxij) go to 4
          Cvxs = Cvxs + float(j)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i)*(d2tdt2*theta**(j-1) + dtdt**2*float(j-1)*theta**(j-2))
4	continue
	Cvxs = -Ti*Cvxs
	Cv = 1000.*Cvxs + Cvig + 1000.*Cvel

	Exs = 0.
	acof = 0.
	bcof = 0.
	do 6 i=0,nobm
	 do 6 j=0,noth
	  if (i+j .ge. maxij) go to 6
	  Exs = Exs + aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i)*theta**(j-1)*(theta - float(j)*Ti*dtdt)
	  if (j .eq. 0.) acof = acof + aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i)
	  if (j .eq. 1.) bcof = bcof + aliq(i+1,j+1)/(dfac(i)*dfac(j))*f**(i)
c	  print*, i,j,aliq(i+1,j+1),f,theta,dtdt
6	continue
	acof = acof - bcof
	bcof = bcof*(1. - fmth)/To**fmth
	do 16 i=0,noln
	 Exs = Exs - aliq(i+1,lineart)*f**i/dfac(i)
16	continue
	E = 1000.*Exs + Eig + 1000.*Eel

	alp = akt/K
	gamma = 1000.*akt*Vi/Cv
	Cp = Cv*(1. + alp*gamma*Ti)
	Ks = K*(1. + alp*gamma*Ti)

	Gsh = d1mach(3)
	dGdT = d1mach(3)
	
        return
        end
