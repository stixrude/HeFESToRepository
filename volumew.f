	double precision function volumew(ispec,rhog)
	include 'P1'
	double precision, parameter :: tmin = 273.15
      LOGICAL qprnt1, qprnt2, qprnt3, qwrphi, qsilent
      LOGICAL qerr,qfail,qrhog,isochor

      CHARACTER(LEN=24) udescr

      INTEGER nttyo, noutpt
      INTEGER iter
	integer ispec,jspec

      REAL(8) tcr, rhocr, pcr, rcnstw
      REAL(8) delta, tau, rho, tempk, press
      REAL(8) ax, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx, gx,
     $ vx, pdx, ptx, adx, gdx, avx, ktx, ztx
      REAL(8) betamx, bettol, btxtol, rhotol, psat, rhog

      REAL(8) xmcapw

	real(8) apar,Ti,Pi

        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ jspec,isochor

      DATA tcr / 647.096d+00 /, rhocr / 322.0d+00 /,
     $ pcr / 22.064d0 /, rcnstw / 0.46151805d+00 /,
     $ xmcapw / 0.018015268d0 /

C  Water equation of state is invalid below 273.15 K as per FDESCR (called by CALPRE)
	if (Ti .lt. tmin) then
	 volumew = -1
	 return
	end if

	qprnt1 = .false.
	qprnt2 = .false.
	qprnt3 = .false.
	qwrphi = .false.
	qsilent = .true.
	noutpt = 999
	nttyo = 6

C  Pi is in GPa.  CALPRE expects pressure in MPa.
	press = 1000.*Pi
	tempk = Ti
	rho = rhog
	delta = rho/rhocr
	tau = tcr/tempk
	qrhog = .false.

      CALL CALPRE(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qwrphi, qrhog, rhog, pcr, rcnstw, rhocr, tcr, iter,
     $ betamx, bettol, rhotol, btxtol, udescr, press, tempk, tau,
     $ delta, rho, psat, px, ux, sx, hx, cvx, cpx, wx, mux, dtx,
     $ bsx, ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx)

	volumew = 1.e6*xmcapw/rho

	return
	end
