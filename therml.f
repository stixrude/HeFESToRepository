        subroutine therml(ispec,Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ent,akt,daktdv,
     &                   beta,Kp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'
	
	integer ispec,i,j,nbm,nobm,noth,mfit,maxij,lineart,noln
	double precision vi,volnl,cp,cv,gamma,alp,ftot,ph,ent,deltas,tcal,zeta,gsh,uth,uto,thet,q,etas
	double precision pzp,akt,aktel,aktig,aktxs,apar,cvel,cvig,cvxs,d2fdv2,d3fdv3,d2tdt2,dgdt
	double precision dfdv,dtdt,e,eel,eig,exs,f,fac,fel
	double precision fig,fmth,fn,fxs,Pi,pig,sel,sig,sxs,theta,Ti,to,vo,dfac
        double precision K,Ks,Kxs,Kig,Kel,Kelp,P,Pxs,Pel,Kxsp,Kp
        double precision daktdv,daktdvxs,daktdvel,beta,betaxs,betael
	double precision aliq(nparp,nparp),aliqc,cliq0,cliq1,cliq2,tee,dteedt,acof,bcof
	double precision betaig,daktdvig,Kigp,d1mach
        double precision bpar,eta,Fo,telo,tinf,wm,xi,zelo
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /liqc/ aliqc(nspecp,nparp,nparp),mfit,nobm,noth,maxij,lineart,nbm,noln
	double precision, parameter :: fsmall=1.e-12

        call liqset(ispec,apar,aliq,fn,wm,To,Fo,Vo,Telo,eta,Tinf,zelo,xi,bpar,fmth)

	f = 0.5*((Vo/Vi)**(2./3.) - 1.)
	if (f .eq. 0.) f = fsmall
	dfdv = -(2.*f + 1.)**(2.5)/(3.*Vo)
	d2fdv2 = 5.*(2.*f + 1)**4/(9.*Vo*Vo)
	d3fdv3 = -40.*(2.*f + 1)**(5.5)/(27.*Vo*Vo*Vo)
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
c	print*, 'in therml F',ispec,Ftot,Fxs,Fel,Fig/1000.

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
c	write(31,*) 'Entropy',Ti,ent,1000.*Sxs+Sig

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

c	print*, 'akt (MPa/K) =',Ti,Pi,Vi,1000.*akt

	daktdvxs = 0.
	do 8 i=0,nobm
	 do 8 j=0,noth
	  if (i+j .ge. maxij) go to 8
	  daktdvxs = daktdvxs + float(i)*float(j)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*theta**(j-1)
     &             *(d2fdv2*f**(i-1) + dfdv**2*float(i-1)*f**(i-2))
8	continue
	daktdvxs = -dtdt*daktdvxs
	do 18 i=1,noln
	 daktdvxs = daktdvxs - dteedt*aliq(i+1,lineart)/dfac(i)*float(i)*(d2fdv2*f**(i-1) + dfdv**2*float(i-1)*f**(i-2))
18	continue
	daktdv = daktdvxs + daktdvig + daktdvel

c	print*, 'daktdv (MPa/(K cm^3/mol)) =',Ti,Pi,Vi,1000.*daktdv

	Kxs = 0.
	do 3 i=0,nobm
	 do 3 j=0,noth
	  if (i+j .ge. maxij) go to 3
	  Kxs = Kxs + float(i)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*theta**(j)*(d2fdv2*f**(i-1) + dfdv**2*float(i-1)*f**(i-2))
3	continue
	do 13 i=1,noln
c	 Kxs = Kxs + d2fdv2*tee*aliq(i+1,lineart)*float(i)*f**(i-1)/dfac(i)
	 Kxs = Kxs + tee*aliq(i+1,lineart)*float(i)/dfac(i)*(d2fdv2*f**(i-1) + dfdv**2*float(i-1)*f**(i-2))
13	continue
	Kxs = Vi*Kxs
	K = Kxs + Kig + Kel

	Kxsp = 0.
	do 9 i=0,nobm
	 do 9 j=0,noth
	  if (i+j .ge. maxij) go to 9
	  Kxsp = Kxsp + float(i)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*theta**(j)
     &         *(d3fdv3*f**(i-1) + 3.*d2fdv2*dfdv*float(i-1)*f**(i-2) + dfdv**3*float(i-1)*float(i-2)*f**(i-3))
9	continue
	do 19 i=1,noln
	 Kxsp = Kxsp + tee*aliq(i+1,lineart)*float(i)/dfac(i)
     &        *(d3fdv3*f**(i-1) + 3.*d2fdv2*dfdv*float(i-1)*f**(i-2) + dfdv**3*float(i-1)*float(i-2)*f**(i-3))
19	continue
	Kxsp = - 1. - Vi**2/Kxs*Kxsp
	Kp = (Kxs*Kxsp + Kig*Kigp + Kel*Kelp)/K

c        print*, 'K prime (-) =',Ti,Pi,Vi,Kp,Kxsp,Kigp,Kelp

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

        betaxs = 0.
        do 20 i=0,nobm
         do 20 j=0,noth
          if (i+j .ge. maxij) go to 20
          betaxs = betaxs + float(j)*aliq(i+1,j+1)/(dfac(i)*dfac(j))*float(i)*f**(i-1)*dfdv
     &           *(d2tdt2*theta**(j-1) + dtdt**2*float(j-1)*theta**(j-2))
20       continue
        betaxs = -Ti*betaxs
        beta = 1000.*betaxs + betaig + 1000.*betael

c	print*, 'beta =',beta

	alp = akt/K
	gamma = 1000.*akt*Vi/Cv
	Cp = Cv*(1. + alp*gamma*Ti)
	Ks = K*(1. + alp*gamma*Ti)

	Gsh = d1mach(3)
	dGdT = d1mach(3)

        return
        end
