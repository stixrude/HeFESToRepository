        subroutine thermw(ispec,Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E)

	include 'P1'

      LOGICAL qprnt1, qprnt2, qprnt3, qwrphi, qsilent

c      INTEGER nttyo, noutpt, ncall
      INTEGER nttyo, noutpt
      INTEGER ispec

      REAL(8) tcr, rhocr, pcr, rcnstw
      REAL(8) delta, tau, rho, tempk, press
      REAL(8) ax, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx, gx,
     $ vx, pdx, ptx, adx, gdx, avx, ktx, ztx

      REAL(8) xmcapw

	double precision Vi,volnl,Cp,Cv,gamma,K,Ks,alp,Ftot,ph,ent,deltas,
     &                   tcal,zeta,Gsh,uth,uto,thet,q,etas,dGdT,pzp,Sel,Eel,Pel,Cvel,Eig,Pig,P,E
	double precision apar,Ti,Pi
	double precision, parameter :: Sconst = 84.164 - 20.82961565		!  Recovers JANAF value at 360 K in J/mol/K
	double precision, parameter :: Fconst = -236.839 - (-19.09516677)	!  Recovers JANAF value of DG_f at 300 K in kJ/mol/K

        common /state/ apar(nspecp,nparp),Ti,Pi

      DATA tcr / 647.096d+00 /, rhocr / 322.0d+00 /,
     $ pcr / 22.064d0 /, rcnstw / 0.46151805d+00 /,
     $ xmcapw / 0.018015268d0 /

c        data ncall/0/
c        ncall = ncall + 1
c        if (ncall .eq. 1) write(31,*) "INFORMATION: WATER FTR OF WAGNER AND PRUSS 2002"

	qprnt1 = .false.
	qprnt2 = .false.
	qprnt3 = .false.
	qwrphi = .false.
	qsilent = .true.
	noutpt = 999
	nttyo = 6

	rho = 1.e6*xmcapw/Vi
	tempk = Ti
	delta = rho/rhocr
	tau = tcr/tempk

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, delta, rho, px, ax, ux, sx, hx,
     $ cvx, cpx, wx, mux, dtx, bsx, gx, vx, pdx, ptx, adx,
     $ gdx, avx, ktx, ztx, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilent)

	Cp = cpx*xmcapw*1000.		! J/mol/K
	Cv = cvx*xmcapw*1000.		! J/mol/K
	gamma = avx/ktx/(rho*cvx)*1000.	! -
	K = 1./ktx/1000.		! GPa
	KS = (wx/1000.)**2*rho/1000.	! GPa
	alp = avx			! 1/K
	Ftot = 1000.*Fconst + ax*xmcapw*1000. - Sconst*Ti		! J/mol
	ent = sx*xmcapw*1000. + Sconst	! J/mol/K
	P = px/1000.			! GPa
	E = ux*xmcapw*1000.		! J/mol
	E = Ftot + ent*Ti		! J/mol

	return
	end
