	subroutine Ftotsubh(ispec,Vi,Ftot)

	include 'P1'
	include 'hydrogen.inc'
        double precision Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E
        double precision apar,Ti,Pi,volve,entve,cpve,bkve,Fpv
	integer ispec
        integer jlop,klop,jlot,klot,i,j,it,jp
	integer, save :: ncall
        integer, parameter :: mint = 2
        double precision logTi,logPi,y,dy,x1,rho,rhop,rhom,volumeh
        double precision ya(mint,mint),y2a(nt,np),y2u(nt,np),y2S(nt,np),y2rhoT(nt,np),y2rhoP(nt,np),y2ST(nt,np)
	double precision, parameter :: eps = 1.e-4
	data ncall/0/
        common /volent/ volve,entve,cpve,bkve
        common /state/ apar(nspecp,nparp),Ti,Pi
	ncall = ncall + 1
	ispec = 1
	x1 = 1.0

        logTi = log10(Ti)
        logPi = log10(Pi)
	Vi = volumeh(ispec,x1)
	rho = log10(hmass/Vi)

C --> Testing
c	Ti = 10**(logTi)
c	Pi = 10**(logPi + eps)
c	Vi = volumeh(ispec,x1)
c	rhop = log10(hmass/Vi)
c
c	Ti = 10**(logTi)
c	Pi = 10**(logPi - eps)
c	Vi = volumeh(ispec,x1)
c	rhom = log10(hmass/Vi)
c
c	Pi = 10**logPi
C <-- Testing

	y = (rhop - rhom)/(2.*eps)
	K = Pi/y

c	if (ncall .eq. 1) call splie2(temph,pressh,rhoh,nt,np,y2a)
c	call splin2(temph,pressh,rhoh,y2a,nt,np,logTi,logPi,y)
c        rho = 10**y
c        Vi = hmass/rho

	if (ncall .eq. 1) call splie2(temph,pressh,uh,nt,np,y2u)
	call splin2(temph,pressh,uh,y2u,nt,np,logTi,logPi,y)
	E = 10**y*1000.*hmass

	if (ncall .eq. 1) call splie2(temph,pressh,Sh,nt,np,y2S)
	call splin2(temph,pressh,Sh,y2S,nt,np,logTi,logPi,y)
	ent = 10**y*1000.*hmass

	if (ncall .eq. 1) call splie2(temph,pressh,dlnrhodlnT,nt,np,y2rhoT)
	call splin2(temph,pressh,dlnrhodlnT,y2rhoT,nt,np,logTi,logPi,y)
	alp = -y/Ti

	if (ncall .eq. 1) call splie2(temph,pressh,dlnrhodlnP,nt,np,y2rhoP)
	call splin2(temph,pressh,dlnrhodlnP,y2rhoP,nt,np,logTi,logPi,y)
	K = Pi/y

	if (ncall .eq. 1) call splie2(temph,pressh,dlnSdlnT,nt,np,y2ST)
	call splin2(temph,pressh,dlnSdlnT,y2ST,nt,np,logTi,logPi,y)
	Cp = y*ent

	gamma = alp*K/(0.001*Cp/Vi - alp**2*K*Ti)

	Cv = Cp/(1. + alp*gamma*Ti)
	KS = K*(1. + alp*gamma*Ti)

	P = Pi

	Ftot = E - ent*Ti
	Fpv = 1000.*Pi*Vi

	Ftot = Ftot + Fpv + hG0

        volve = Vi
        entve = ent      
        cpve = Cp
        bkve = K

	return
	end
