      SUBROUTINE EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, delta, rho, px, ax, ux, sx, hx,
     $ cvx, cpx, wx, mux, dtx, bsx, gx, vx, pdx, ptx, adx,
     $ gdx, avx, ktx, ztx, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilent)

c     Evaluate the basic equation of state, which is written as
c     a function of:

c       reduced density: delta = rho/rhoc
c       inverse reduced temperature: tau = tcr/tempk

c     The pressure px is an output. To obtain properties as
c     a function of temperature and pressure, it is necessary to
c     iterate on pressure (adjusting density).

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qprnt1, qprnt2, qprnt3, qwrphi, qsilent

      INTEGER nttyo, noutpt

      REAL(8) tcr, rhocr, pcr, rcnstw
      REAL(8) delta, tau, rho, tempk, press
      REAL(8) ax, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx, gx,
     $ vx, pdx, ptx, adx, gdx, avx, ktx, ztx

      REAL(8) cin0(1:8), cigam0(4:8)

c     Original Wagner and Pruss (2002) values for cin0(1) and cin0(2).

c     DATA cin0(1:8) / -8.32044648201d+00, 6.6832105268d+00,
c    $ 3.00632d+00, 0.012436d+00, 0.97315d+00, 1.27950d+00,
c    $ 0.96956d+00, 0.24873d+00 /

c     IAPWS (2016) revised values for cin0(1) and cin0(2).

      DATA cin0(1:8) / -8.3204464837497d+00, 6.6832105275932d+00,
     $ 3.00632d+00, 0.012436d+00, 0.97315d+00, 1.27950d+00,
     $ 0.96956d+00, 0.24873d+00 /

      DATA cigam0(4:8) / 1.28728967d+00, 3.53734222d+00,
     $ 7.74073708d+00, 9.24437796d+00, 27.5075105d+00 /

c     Note: cigam0(1:3) are not used in the EOS model.

c     Coefficients for the residual part.

c     First part, indices 1-51.

      INTEGER c1c(8:54), c1d(1:54)

      REAL(8) c1t(1:54), c1n(1:54)

      DATA c1c(8:51) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $  3, 3, 3, 3, 4, 6, 6, 6, 6 /

c     Note: c1c(1:7) are not used in the EOS model.

      DATA c1d(1:51) / 1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5,
     $ 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9,
     $ 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6 /

      DATA c1t(1:51) / -0.5d0, 0.875d0, 1.0d0, 0.5d0, 0.75d0, 0.375d0,
     $ 1.0d0, 4.0d0, 6.0d0, 12.0d0, 1.0d0, 5.0d0, 4.0d0, 2.0d0,
     $ 13.0d0, 9.0d0, 3.0d0, 4.0d0, 11.0d0, 4.0d0, 13.0d0, 1.0d0,
     $ 7.0d0, 1.0d0, 9.0d0, 10.0d0, 10.0d0, 3.0d0, 7.0d0, 10.0d0,
     $ 10.0d0, 6.0d0, 10.0d0, 10.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0,
     $ 8.0d0, 6.0d0, 9.0d0, 8.0d0, 16.0d0, 22.0d0, 23.0d0, 23.0d0,
     $ 10.0d0, 50.0d0, 44.0d0, 46.0d0, 50.0d0 /

      DATA c1n(1:51) /  0.12533547935523d-01,
     $  0.78957634722828d+01, -0.87803203303561d+01,
     $  0.31802509345418d+00, -0.26145533859358d+00,
     $ -0.78199751687981d-02,  0.88089493102134d-02,
     $ -0.66856572307965d+00,  0.20433810950965d+00,
     $ -0.66212605039687d-04, -0.19232721156002d+00,
     $ -0.25709043003438d+00,  0.16074868486251d+00,
     $ -0.40092828925807d-01,  0.39343422603254d-06,
     $ -0.75941377088144d-05,  0.56250979351888d-03,
     $ -0.15608652257135d-04,  0.11537996422951d-08,
     $  0.36582165144204d-06, -0.13251180074668d-11,
     $ -0.62639586912454d-09, -0.10793600908932d+00,
     $  0.17611491008752d-01,  0.22132295167546d+00,
     $ -0.40247669763528d+00,  0.58083399985759d+00,
     $  0.49969146990806d-02, -0.31358700712549d-01,
     $ -0.74315929710341d+00,  0.47807329915480d+00,
     $  0.20527940895948d-01, -0.13636435110343d+00,
     $  0.14180634400617d-01,  0.83326504880713d-02,
     $ -0.29052336009585d-01,  0.38615085574206d-01,
     $ -0.20393486513704d-01, -0.16554050063734d-02,
     $  0.19955571979541d-02,  0.15870308324157d-03,
     $ -0.16388568342530d-04,  0.43613615723811d-01,
     $  0.34994005463765d-01, -0.76788197844621d-01,
     $  0.22446277332006d-01, -0.62689710414685d-04,
     $ -0.55711118565645d-09, -0.19905718354408d+00,
     $  0.31777497330738d+00, -0.11841182425981d+00 /

c     Second part, indices 52-54.

      INTEGER c2d(52:54), c2t(52:54), c2alph(52:54), c2beta(52:54),
     $ c2eps(52:54)

      REAL(8) c2n(52:54), c2gamm(52:54)

c     Note: c2c(52:54) are not used in the EOS model.

      DATA c2d(52:54) / 3, 3, 3 /
      DATA c2t(52:54) / 0, 1, 4 /
      DATA c2n(52:54) / -0.31306260323435d+02, 0.31546140237781d+02,
     $ -0.25213154341695d+04 /

      DATA c2alph(52:54) / 20, 20, 20 /
      DATA c2beta(52:54)  / 150, 150, 250 /
      DATA c2gamm(52:54) / 1.21d0, 1.21d0, 1.25d0 /
      DATA c2eps(52:54)   / 1, 1, 1 /

c     Third part, indices 55-56.

      INTEGER c3ccap(55:56), c3dcap(55:56)

      REAL(8) c3a(55:56), c3b(55:56), c3bcap(55:56), c3n(55:56),
     $ c3acap(55:56), c3beta(55:56)

      DATA c3a(55:56) / 3.5d0, 3.50d0 /
      DATA c3b(55:56) / 0.85d0, 0.95d0 /
      DATA c3bcap(55:56) / 0.2d0, 0.2d0 /
      DATA c3n(55:56) / -0.14874640856724d0, 0.31806110878444d0 /
      DATA c3ccap(55:56) / 28, 32 /
      DATA c3dcap(55:56) / 700, 800 /
      DATA c3acap(55:56) / 0.32d0, 0.32d0 /
      DATA c3beta(55:56) / 0.3d0, 0.3d0 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      REAL(8) phi0, phi0d, phi0dd, phi0t, phi0tt, phi0dt

      REAL(8) phir, phird, phirdd, phirt, phirtt, phirdt

      REAL(8) delsq, tausq

      REAL(8) egam0(4:8), xegam0(4:8)

      REAL(8) deldi(1:54),tauti(1:54)
      REAL(8) deldm1(1:54), deldm2(1:54)
      REAL(8) delci(8:51), dmci(8:51), edci(8:51)
      REAL(8) tautm1(1:51), tautm2(1:51)

      REAL(8) exxi(52:54),dade(52:54),tbtg(52:54)

      REAL(8) TH(55:56), HK(55:56), HKK(55:56), BEX(55:56)
      REAL(8) WW(55:56), WWd(55:56), WWdd(55:56)
      REAL(8) HH(55:56), HHd(55:56), HHdd(55:56), HHt(55:56),
     $ HHtt(55:56), HHdt(55:56)
      REAL(8) PS(55:56), PSd(55:56), PSdd(55:56), PSt(55:56),
     $ PStt(55:56), PSdt(55:56)

      REAL(8) XC3B, XC3B2

      REAL(8) dm1, dm1sq, tm1, tm1sq

      REAL(8) epxc, rtx, x1, x2, x3, x4, xsum, xxt, ztxc

      INTEGER i,ix

	real(8) nx(2),dx(2),tx(2),deltax,rhox
c	data dx / 3.0d0 /
c	data tx / 1.5d0 /
c	data nx / -2.00d+01 /
c	data tx / 1.2d0 /
c	data nx / -8.80d+00 /
c	data tx / 1.0d0 /
c	data nx / -5.05d+00 /
c	data tx / 0.5d0 /
c	data nx / -1.42d+00 /
c	data tx / 0.0d0 /
c	data nx / -4.82d-01 /
c	data nx / 0.00d-01 /

C  No additional term
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 0.0d0 , 0.0d0/
c	data nx / -0.00d+00, 0.00d-00 /
c Density Differences and RMS residual (g/cc)	-0.113109 8.87693e-05 -0.264789 	RMS residual (g/cc) 0.287936

C + dx^3
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 1.0d0 , 0.0d0/
c	data nx / -2.50d-00, 0.00d-00 /
c Density Differences and RMS residual (g/cc)	-0.031551 0.136684 -0.123214 	RMS residual (g/cc) 0.186707

C + dx^3 + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 1.0d0 , 1.0d0/
c	data nx / +3.10d-01, -1.40d-00 /
c Density Differences and RMS residual (g/cc)	-0.0743586 0.136528 -0.0941778 	RMS residual (g/cc) 0.181765

C + dx^3*T + dx^4*T
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 0.0d0 , 0.0d0/
c	data nx / +5.7d-01, -3.20d-01 /
c Density Differences and RMS residual (g/cc)	-0.155289 0.0273835 -0.202376 	RMS residual (g/cc) 0.256555

C + dx^3*T + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 0.0d0 , 1.0d0/
c	data nx / +1.70d-01, -2.80d+00 /
c Density Differences and RMS residual (g/cc)	-0.0492207 -0.0209655 0.0295288 	RMS residual (g/cc) 0.0611079

C + dx^4*T + dx^5
c	data dx / 4.0d0, 5.0d0 /
c	data tx / 0.0d0 , 1.0d0/
c	data nx / +1.50d-01, -1.70d-00 /
c Density Differences and RMS residual (g/cc)	-0.0900889 -0.136228 0.0486411 	RMS residual (g/cc) 0.170411

C + dx^3*T^1.5 + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / -0.5d0 , 1.0d0/
c	data nx / +2.30d-02, -2.30d+00 /
c Density Differences and RMS residual (g/cc)	-0.0401839 -0.00808412 0.0236153 	RMS residual (g/cc) 0.0473052

C + dx^3*T^1.5 + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / -0.5d0 , 1.0d0/
c	data nx / 2.73483396E-02 , -2.85302258/
c Density Differences / Error and RMS residual (g/cc)     -1.34218 -0.136719 2.13162      Chi Squared 6.36394

C + dx^3*T + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 0.0d0 , 1.0d0/
c	data nx / 0.169315994  ,    -3.28534436  /
c Density Differences / Error and RMS residual (g/cc)	-2.29092 0.912064 2.40161 	Chi Squared 11.8479

C + dx^4*T + dx^5
c	data dx / 4.0d0, 5.0d0 /
c	data tx / 0.0d0 , 1.0d0/
c	data nx / 8.18269998E-02 ,  -1.44293773 /
c Density Differences / Error and RMS residual (g/cc)	-8.14652 -0.343852 2.92626 	Chi Squared 75.047

C + dx^3*T^2 + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / -1.0d0 , 1.0d0/
c	data nx / 4.14401293E-03 , -2.58466506/
c Density Differences / Error and RMS residual (g/cc)	-1.31474 -0.111378 2.04633 	Chi Squared 5.92841

C + dx^3 + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 1.0d0 , 1.0d0/
c       data nx / +3.10d-01, -1.40d-00 /

C + dx^3
	data dx / 3.0d0, 4.0d0 /
	data tx / 1.0d0 , 0.0d0/
	data nx / -2.74160981   ,   0.00000000  /
c Density Differences / Error and RMS residual (g/cc)	-2.16322 5.07193 -1.80081 	Chi Squared 33.6469

C + dx^3
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 1.0d0 , 1.0d0/
c	data nx / -2.5   ,   0.00000000  /
c Density Differences / Error and RMS residual (g/cc)	-3.1551 4.55612 -2.05357 	Chi Squared 34.93

C + dx^4
c	data dx / 4.0d0, 4.0d0 /
c	data tx / 1.0d0 , 0.0d0/
c	data nx / -1.47413874       0.00000000  /
c Density Differences / Error and RMS residual (g/cc)	-6.02431 5.66495 -0.978558 	Chi Squared 69.3415

C + dx^3 + dx^4
c	data dx / 3.0d0, 4.0d0 /
c	data tx / 1.0d0 , 1.0d0/
c	data nx / -6.22489786   ,   2.46272516  /
c Density Differences / Error and RMS residual (g/cc)	-0.266995 2.49055 -3.68727 	Chi Squared 19.8701

	data rhox / 1237.39d0 /
	deltax = delta*rhocr/rhox
c	open(151,file='nxvalues.txt',status='old')
c	read(151,*) nx(1),nx(2)
c	close (151)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      epxc = 100*EPSILON(delta)

c     For 64-bit arithmetic, 100*EPS is approximately 2.2204d-14.

c     Avoid a singularity for delta too close to zero.

      IF (delta .LT. epxc) THEN
        delta = epxc
        rho = rhocr*delta
      ENDIF

c     Avoid a singularity at the critical density (at any temperature).
c     If delta = rho/rhocr is unity, then delta - 1 is zero. This
c     cannot be used with (delta - 1)**n where n is negative, or in a
c     division (x/(d -1)).

      x1 = 1.0d0 - epxc
      x2 = 1.0d0 + epxc

      IF (delta.GT.x1 .AND. delta.LT.x2) THEN

c       The density is too close to the critical density.

        IF (qprnt1) THEN
          IF (.NOT.qsilent) THEN

            WRITE(nttyo,1110)
            WRITE(noutpt,1110)

 1110       FORMAT(/3x,'TOO CLOSE TO THE CRITICAL DENSITY,',
     $      ' ADJUSTING.',/)

          ENDIF
        ENDIF

        IF (delta .LT. 1.0d0) THEN
          delta = x1
        ELSE
          delta = x2
        ENDIF

        rho = rhocr*delta

      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Zero parameters to be calculated.

      px = 0.0d0
      ax = 0.0d0
      ux = 0.0d0
      sx = 0.0d0
      hx = 0.0d0
      cvx = 0.0d0
      cpx = 0.0d0
      wx = 0.0d0
      mux = 0.0d0
      dtx = 0.0d0
      bsx = 0.0d0
      gx = 0.0d0
      vx = 0.0d0
      pdx = 0.0d0
      ptx = 0.0d0
      adx = 0.0d0
      gdx = 0.0d0
      avx = 0.0d0
      ktx = 0.0d0
      ztx = 0.0d0

c     Calculate some common pieces.

      delsq = delta**2
      tausq = tau**2
      rtx = rcnstw*tempk

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate the ideal-gas parts.

c     Calculate common pieces first.

      DO i = 4,8
        egam0(i) = DEXP( -cigam0(i)*tau )
        xegam0(i) = 1.0d0 - egam0(i)
      ENDDO

c     Now calculate phi0 and its partial derivatives.

      phi0 = DLOG(delta) + cin0(1) + cin0(2)*tau
     $ + cin0(3)*DLOG(tau)
      DO i = 4,8
        phi0 = phi0 + cin0(i)*DLOG(xegam0(i))
      ENDDO

      phi0d = 1.0d0/delta

      phi0dd = -1.0d0/delsq

      phi0t = cin0(2) + (cin0(3)/tau)
      DO i = 4,8
        phi0t = phi0t + cin0(i)*cigam0(i)*( (xegam0(i)**(-1)) -1.0d0 )
      ENDDO

      phi0tt = -cin0(3)/tausq
      DO i = 4,8
        phi0tt = phi0tt
     $  - cin0(i)*(cigam0(i)**2)*egam0(i)*( xegam0(i)**(-2) )
      ENDDO

      phi0dt = 0.0d0

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate the residual parts.

c     Calculate common pieces.

      dm1 = delta - 1.0d0
      dm1sq = dm1**2
      tm1 = tau - 1.0d0
      tm1sq = tm1**2

      DO i = 1,51
        deldi(i) = delta**c1d(i)
        tauti(i) = tau**c1t(i)
        deldm1(i) = delta**(c1d(i) - 1.0d0)
        deldm2(i) = delta**(c1d(i) - 2.0d0)
        tautm1(i) = tau**(c1t(i) - 1.0d0)
        tautm2(i) = tau**(c1t(i) - 2.0d0)
      ENDDO

      DO i = 8,51
        delci(i) = delta**c1c(i)
        edci(i) = DEXP(-delci(i))
        dmci(i) = c1d(i) - c1c(i)*delci(i)
      ENDDO

      DO i = 52,54
        deldi(i) = delta**c2d(i)
        tauti(i) = tau**c2t(i)
        deldm1(i) = delta**(c2d(i) - 1.0d0)
        deldm2(i) = delta**(c2d(i) - 2.0d0)
        exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
     $  - c2beta(i)*(tau - c2gamm(i))**2)
        dade(i) = (c2d(i)/delta) -2.0d0*c2alph(i)*(delta - c2eps(i))
        tbtg(i) = (c2t(i)/tau) -2.0d0*c2beta(i)*(tau - c2gamm(i))
      ENDDO

      DO i = 55,56
        XC3B = 1.0d0/( 2.0d0*c3beta(i) )

        BEX(i) = XC3B - 1.0d0
        TH(i) = (1.0d0 - tau) + c3acap(i)*dm1sq**XC3B
        PS(i) = DEXP( -c3ccap(i)*dm1sq - c3dcap(i)*tm1**2 )

        WW(i) = TH(i)**2 + c3bcap(i)*dm1sq**c3a(i)
        HH(i) = WW(i)**c3b(i)
        HK(i) = WW(i)**(c3b(i) - 1.0d0)
        HKK(i) = WW(i)**(c3b(i) - 2.0d0)

        WWd(i) = dm1*( c3acap(i)*TH(i)*(2.0d0/c3beta(i))*dm1sq**BEX(i)
     $  + 2.0d0*c3bcap(i)*c3a(i)*dm1sq**(c3a(i) - 1.0d0) )

        x1 = 4.0d0*c3bcap(i)*c3a(i)*(c3a(i) - 1.0d0)
     $  *dm1sq**(c3a(i) - 2.0d0)

        XC3B2 = (1.0d0/c3beta(i))**2
        x2 = 2.0d0*(c3acap(i)**2)*XC3B2*( dm1sq**BEX(i) )**2

        x3 = c3acap(i)*TH(i)*( 4.0d0/c3beta(i) )*BEX(i)
     $  *dm1sq**(BEX(i) - 1.0d0)

        xsum = x1 + x2 + x3
        WWdd(i) = (WWd(i)/dm1) + dm1sq*xsum

        HHd(i) = c3b(i)*HK(i)*WWd(i)

        HHdd(i) = c3b(i)*( HK(i)*WWdd(i) + (c3b(i) - 1.0d0)
     $  *HKK(i)*( WWd(i)**2 ) )

        HHt(i) = -2.0d0*TH(i)*c3b(i)*HK(i)

        HHtt(i) = 2.0d0*c3b(i)*HK(i) + 4.0d0*(TH(i)**2)*c3b(i)
     $  *(c3b(i) - 1.0d0)*HKK(i)

        x1 = 2.0d0/c3beta(i)
        HHdt(i) = -c3acap(i)*c3b(i)*x1*HK(i)*dm1*(dm1sq**BEX(i))
     $  - 2.0d0*TH(i)*c3b(i)*(c3b(i) - 1.0d0)*HKK(i)*WWd(i)

        PSd(i) = -2.0d0*c3ccap(i)*dm1*PS(i)
        PSdd(i) = ( 2.0d0*c3ccap(i)*dm1sq - 1.0d0 )
     $  *2.0d0*c3ccap(i)*PS(i)
        PSt(i) = -2.0d0*c3dcap(i)*tm1*PS(i)
        PStt(i) = ( 2.0d0*c3dcap(i)*tm1sq -1 )*2.0d0*c3dcap(i)*PS(i)
        PSdt(i) = 4.0d0*c3ccap(i)*c3dcap(i)*dm1*tm1*PS(i)

      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Now calculate phir and its partial derivatives.

c     Calculate phir.

      phir = 0.0d0

      DO i =1,7

c       deldi(i) = delta**c1d(i)
c       tauti(i) = tau**c1t(i)

        phir = phir + c1n(i)*deldi(i)*tauti(i)
      ENDDO

      DO i = 8,51

c       delci(i) = delta**c1c(i)
c       deldi(i) = delta**c1d(i)
c       tauti(i) = tau**c1t(i)
c       dmci(i) = c1d(i) - c1c(i)*delci(i)
c       edci(i) = DEXP(-delci(i))

        phir = phir + c1n(i)*deldi(i)*tauti(i)*edci(i)
      ENDDO

      DO i = 52,54

c       deldi(i) = delta**c2d(i)
c       tauti(i) = tau**c2t(i)
c       exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
c    $  - c2beta(i)*(tau - c2gamm(i))**2)

        phir = phir + c2n(i)*deldi(i)*tauti(i)*exxi(i)
      ENDDO

      DO i = 55,56

c       XC3B = 1.0d0/( 2.0d0*c3beta(i) )
c       TH(i) = (1.0d0 - tau) + c3acap(i)*dm1sq**XC3B
c       WW(i) = TH(i)**2 + c3bcap(i)*dm1sq**c3a(i)
c       HH(i) = WW(i)**c3b(i)
c       PS(i) = DEXP( -c3ccap(i)*dm1sq - c3dcap(i)*tm1**2 )

        phir = phir + c3n(i)*HH(i)*delta*PS(i)
      ENDDO

	if (deltax .ge. 1.0d0) then
	 do ix=1,2
	  phir = phir + nx(ix)*(deltax - 1.0d0)**dx(ix)*tau**tx(ix)
	 enddo
	end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phird.

      phird = 0.0d0

      DO i =1,7

c       deldm1(i) = delta**(c1d(i) - 1.0d0)
c       tauti(i) = tau**c1t(i)

        phird = phird + c1n(i)*c1d(i)*deldm1(i)*tauti(i)
      ENDDO

      DO i = 8,51

c       delci(i) = delta**c1c(i)
c       edci(i) = DEXP(-delci(i))
c       deldm1(i) = delta**(c1d(i) - 1.0d0)
c       tauti(i) = tau**c1t(i)
c       dmci(i) = c1d(i) - c1c(i)*delci(i)

        phird = phird + c1n(i)*edci(i)*deldm1(i)*tauti(i)*dmci(i)
      ENDDO

      DO i = 52,54

c       deldi(i) = delta**c2d(i)
c       tauti(i) = tau**c2t(i)
c       exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
c    $  - c2beta(i)*(tau - c2gamm(i))**2)
c       dade(i) = (c2d(i)/delta) -2.0d0*c2alph(i)*(delta - c2eps(i))

        phird = phird + c2n(i)*deldi(i)*tauti(i)*exxi(i)*dade(i)
      ENDDO

      DO i = 55,56

c       XC3B = 1.0d0/( 2.0d0*c3beta(i) )
c       BEX(i) = XC3B - 1.0d0
c       TH(i) = (1.0d0 - tau) + c3acap(i)*dm1sq**XC3B
c       WW(i) = TH(i)**2 + c3bcap(i)*dm1sq**c3a(i)
c       HH(i) = WW(i)**c3b(i)
c       PS(i) = DEXP( -c3ccap(i)*dm1sq - c3dcap(i)*tm1**2 )
c       HK(i) = WW(i)**(c3b(i) - 1.0d0)
c       PSd(i) = -2.0d0*c3ccap(i)*dm1*PS(i)
c       WWd(i) = dm1*( c3acap(i)*TH(i)*(2.0d0/c3beta(i))*dm1sq**BEX(i)
c    $  + 2.0d0*c3bcap(i)*c3a(i)*( dm1sq**(c3a(i) - 1.0d0) ) )
c       HHd(i) = c3b(i)*HK(i)*WWd(i)

        phird = phird + c3n(i)*( HH(i)*(PS(i) + delta*PSd(i))
     $  + HHd(i)*delta*PS(i) )
      ENDDO

	if (deltax .ge. 1.0d0) then
	 do ix=1,2
	  phird = phird + nx(ix)*dx(ix)*(deltax - 1.0d0)**(dx(ix) - 1.0d0)*rhocr/rhox*tau**tx(ix)
	 enddo
	end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirdd.

      phirdd = 0.0d0

      DO i =1,7
        phirdd = phirdd + c1n(i)*c1d(i)*(c1d(i) - 1.0d0)
     $  *deldm2(i)*tauti(i)
      ENDDO

      DO i = 8,51

c       delci(i) = delta**c1c(i)
c       edci(i) = DEXP(-delci(i))
c       dmci(i) = c1d(i) - c1c(i)*delci(i)

        xxt = dmci(i)*(dmci(i) - 1.0d0) - (c1c(i)**2)*delci(i)
        phirdd = phirdd + c1n(i)*edci(i)*deldm2(i)*tauti(i)*xxt
      ENDDO

      DO i = 52,54

c       exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
c    $  - c2beta(i)*(tau - c2gamm(i))**2)

        x1 = -2.0d0*c2alph(i)*deldi(i)
        x2 = 4.0d0*(c2alph(i)**2)*deldi(i)*(delta - c2eps(i))**2
        x3 = -4.0d0*c2d(i)*c2alph(i)*deldm1(i)*(delta - c2eps(i))
        x4 = c2d(i)*(c2d(i) - 1.0d0)*deldm2(i)
        xsum = x1 + x2 + x3 + x4
        phirdd = phirdd + c2n(i)*tauti(i)*exxi(i)*xsum
      ENDDO

      DO i = 55,56
        x1 = HH(i)*( 2.0d0*PSd(i) + delta*PSdd(i) )
        x2 = 2.0d0*HHd(i)*( PS(i) + delta*PSd(i) )
        x3 = HHdd(i)*delta*PS(i)
        xsum = x1 + x2 + x3
        phirdd = phirdd + c3n(i)*xsum
      ENDDO

	if (deltax .ge. 1.0d0) then
	 do ix=1,2
	  phirdd = phirdd + nx(ix)*dx(ix)*(dx(ix) - 1.0d0)*(deltax - 1.0d0)**(dx(ix) - 2.0d0)*(rhocr/rhox)**2*tau**tx(ix)
	 enddo
	end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirt.

      phirt = 0.0d0

      DO i =1,7
        phirt = phirt + c1n(i)*c1t(i)*deldi(i)*tautm1(i)
      ENDDO

      DO i = 8,51
        phirt = phirt + c1n(i)*c1t(i)*deldi(i)*tautm1(i)*edci(i)
      ENDDO

      DO i = 52,54
        phirt = phirt + c2n(i)*deldi(i)*tauti(i)*exxi(i)*tbtg(i)
      ENDDO

      DO i = 55,56
        xsum = HHt(i)*PS(i) + HH(i)*PSt(i)
        phirt = phirt + c3n(i)*delta*xsum
      ENDDO

	if (deltax .ge. 1.0d0) then
	 do ix=1,2
	  phirt = phirt + nx(ix)*(deltax - 1.0d0)**dx(ix)*tx(ix)*tau**(tx(ix) - 1.0d0)
	 enddo
	end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirtt.

      phirtt = 0.0d0

      DO i =1,7
        phirtt = phirtt + c1n(i)*c1t(i)*(c1t(i) - 1.0d0)
     $  *deldi(i)*tautm2(i)
      ENDDO

      DO i = 8,51
        phirtt = phirtt + c1n(i)*c1t(i)*(c1t(i) - 1.0d0)*deldi(i)
     $  *tautm2(i)*edci(i)
      ENDDO

      DO i = 52,54
        x1 = tbtg(i)**2
        x2 = -( c2t(i)/tausq ) - 2.0d0*c2beta(i)
        xsum = x1 + x2
        phirtt = phirtt + c2n(i)*deldi(i)*tauti(i)*exxi(i)*xsum
      ENDDO

      DO i = 55,56
        xsum = HHtt(i)*PS(i) + 2.0d0*HHt(i)*PSt(i) + HH(i)*PStt(i)
        phirtt = phirtt + c3n(i)*delta*xsum
      ENDDO

	if (deltax .ge. 1.0d0) then
	 do ix=1,2
	  phirtt = phirtt + nx(ix)*(deltax - 1.0d0)**dx(ix)*tx(ix)*(tx(ix) - 1.0d0)*tau**(tx(ix) - 2.0d0)
	 enddo
	end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirdt.

      phirdt = 0.0d0

      DO i =1,7
        phirdt = phirdt + c1n(i)*c1d(i)*c1t(i)*deldm1(i)*tautm1(i)
      ENDDO

      DO i = 8,51
        phirdt = phirdt + c1n(i)*c1t(i)*deldm1(i)*tautm1(i)
     $  *dmci(i)*edci(i)
      ENDDO

      DO i = 52,54
        xsum = dade(i)*tbtg(i)
        phirdt = phirdt + c2n(i)*deldi(i)*tauti(i)*exxi(i)*xsum
      ENDDO

      DO i = 55,56
        x1 = HH(i)*( PSt(i) + delta*PSdt(i) )
        x2 = delta*HHd(i)*PSt(i)
        x3 = HHt(i)*( PS(i) + delta*PSd(i) )
        x4 = HHdt(i)*delta*PS(i)
        xsum = x1 + x2 + x3 + x4
        phirdt = phirdt + c3n(i)*xsum
      ENDDO

	if (deltax .ge. 1.0d0) then
	 do ix=1,2
	  phirdt = phirdt + nx(ix)*dx(ix)*(deltax - 1.0d0)**(dx(ix) - 1.0d0)*rhocr/rhox*tx(ix)*tau**(tx(ix) - 1.0d0)
	 enddo
	end if

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate thermodynamic functions.

c     Helmholtz energy. The value is in kJ/kg-K.

      ax = rtx*( phi0 + phir )

c     Pressure. The value here is in kPa.

      px = rho*rtx*( 1.0d0 + delta*phird )

c     Internal energy.

      ux = rtx*tau*( phi0t + phirt )

c     Entropy.

      sx = rcnstw*( tau*(phi0t + phirt) - phi0 - phir )

c     Enthalpy.

      hx = rtx*( 1.0d0 + tau*(phi0t + phirt) + delta*phird )

c     Gibbs energy.

      gx = rtx*( 1.0d0 + phi0 + phir + delta*phird )

c     Alternate formulas for the Gibbs energy.

c     gx = hx - tempk*sx
c     gx = ax + hx - ux

c     Volume.

      vx = 1.0d0/rho

c     Isochoric heat capacity.

      cvx = -rcnstw*tausq*( phi0tt + phirtt )

c     Isobaric heat capacity.

      x1 = ( 1.0d0 + delta*phird - delta*tau*phirdt )**2
      x2 = 1.0d0 + 2.0d0*delta*phird + delsq*phirdd
      IF (DABS(x2) .GT. 1.0d-15) THEN
        cpx = cvx + rcnstw*(x1/x2)
      ELSE
        cpx = 1.0d+100
      ENDIF

c     Speed of sound.

      x1 = ( 1.0d0 + delta*phird - delta*tau*phirdt )**2
      x2 = tausq*( phi0tt + phirtt )
      x3 = x1/x2
      xxt = rtx*( 1.0d0 + 2.0d0*delta*phird
     $ + delsq*phirdd - x3 )
      IF (xxt .GT. 0.0d0) THEN
        wx = sqrt(xxt)
      ELSE
        wx = 0.0d0
      ENDIF

c     Joule-Thomsen coefficient.

      x1 = delta*phird + delsq*phirdd + delta*tau*phirdt
      x2 = ( 1.0d0 + delta*phird - delta*tau*phirdt )**2
      x3 = ( phi0tt + phirtt )*( 1.0d0 + 2.0d0
     $ *delta*phird + delsq*phirdd )
      mux = ( - x1/( x2 - tausq*x3 ) )/(rcnstw*rho)

c     Isothermal throttling coefficient.

      x1 = 1.0d0 + delta*phird - delta*tau*phirdt
      x2 = 1.0d0 + 2.0d0*delta*phird + delsq*phirdd
      dtx = ( 1.0d0 - ( x1/x2 ) )/rho

c     Isentropic temperature-pressure coefficient.

      x1 = 1.0d0 + delta*phird - delta*tau*phirdt
      x2 = x1**2
      x3 = ( phi0tt + phirtt )*( 1.0d0 + 2.0d0*delta*phird
     $ + delsq*phirdd )
      bsx = ( x1/( x2 - tausq*x3 ) )/(rcnstw*rho)

c     Derivative of pressure with respect to delta
c     (needed to perform Newton-Raphson iteration
c     to matched desired pressure). Here the value
c     is in kPa. Recall that:
c     px = rho*rtx*( 1.0d0 + delta*phird )

      pdx = ( px/delta ) + delta*rhocr*rtx*( phird + delta*phirdd )

c     Derivative of pressure with respect to tau (needed to calculate
c     the thermal expansion coefficient). Here the value is again
c     in kPa.

      ptx = ( -px/tau ) + px*delta*phirdt/(1.0d0 + delta*phird)

c     Compressibility. Here the value is in /Kpa.

      ktx = 1.0d0/(delta*pdx)

c     Calculate ztx (needed to calculate viscosity). Note: pcr here
c     is in MPa, but pdx is still in /kPa. Hence there is the need
c     to include the factor of 1000 kPa/MPa. Note that ztx itself
c     is dimensionless, so pressure unit correction needs to be made
c     here.

      ztx = 1000.0d0*pcr/pdx

c     An alternative formula is: ztx = 1000.0d0*delta*pcr*ktx
c     This can be useful for getting zeta from a calculator that gives
c     the compressibility (but not zeta).

c     Thermal expansion coefficient (thermal expansivity).
c     This calculation is based on the Maxwell relation:
c     (del P/del T) at const V = alpha/kappa

      avx = ktx*ptx*( -tau/tempk )

c     Parts needed to calculate residuals and Jacobian elements
c     if refining saturation properties at fixed temperature.

c     Helmholtz energy.

c     ax = rtx*( phi0 + phir )

      adx = rtx*( phi0d + phird )

c     Gibbs energy.

c     gx = rtx*( 1.0d0 + phi0 + phir + delta*phird )

      gdx = rtx*( phi0d + 2.0d0*phird + delta*phirdd )

      IF (qwrphi) THEN
        IF (.NOT.qsilent) THEN

          WRITE(noutpt,1120) tempk, rho, delta, tau

 1120     FORMAT(/7x,'tempk = ',e16.9,3x,'rho   = ',e16.9,
     $    /7x,'delta = ',e16.9,3x,'tau   = ',e16.9,//)

          WRITE(noutpt,1130) phi0, phir, phi0d, phird, phi0dd, phirdd,
     $    phi0t, phirt, phi0tt, phirtt, phi0dt, phirdt

 1130     FORMAT(/12x,'phi functions',/
     $    /5x,'phi0   = ',e16.9,5x,'phir   = ',e16.9,
     $    /5x,'phi0d  = ',e16.9,5x,'phird  = ',e16.9,
     $    /5x,'phi0dd = ',e16.9,5x,'phirdd = ',e16.9,
     $    /5x,'phi0t  = ',e16.9,5x,'phirt  = ',e16.9,
     $    /5x,'phi0tt = ',e16.9,5x,'phirtt = ',e16.9,
     $    /5x,'phi0dt = ',e16.9,5x,'phirdt = ',e16.9,//)

        ENDIF
      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Unit conversions. Direct evaluation of the IAPWS-95
c     model equations gives some of the results in in odd
c     units. Conversions are therefore desirable for most
c     usage. Prior to conversion, the units are determined
c     by the units adopted by the model for temperature (K),
c     density (kg/m3), and the gas constant (kJ/kg/K).
c     See Table 3.

c       Note that density is input as kg/m3.
c       Volume is therefore in m3/kg.

c       Pressure comes out in kPA. Divide by 1000 to obtain
c       pressure in MPa. One could divide instead by 100 to
c       obtain pressure in bars. Here we will use MPa.
c       The pdx (the partial derivative of pressure with respect
c       to delta) and ktx (the isothermal compressibility)
c       must also be corrected to be consistent with
c       pressure in MPa.

c       Internal energy, enthalpy, and Helmholtz energy come
c       out in kJ/kg.

c       Entropy and the two heat capacity functions come out
c       in kJ/kg/K.

c       Gibbs energy is therefore in kJ/kg.

c       Speed of sound comes out in sqrt(kJ/kg). Multiply by
c       sqrt(1000) to obtain the speed of sound in m/s.

c       Joule-Thomson coefficient comes out in K-m3/kJ.
c       These units are equivalent to the usual K/MPa.

c       Isothermal throttling coefficient comes out in m3/kg.
c       Divide by 1000 to obtain the result in the usual
c       kJ/kg/MPa.

c       Isentropic temperature-pressure coefficient comes out
c       K-m3/kJ (the same units as for the Joule-Thomson coefficient).
c       These units are equivalent to the usual K/MPa.

      px = 0.001d0*px
      pdx = 0.001d0*pdx
      ptx = 0.001d0*ptx
      ktx = 1000.0d0*ktx

      wx = wx*sqrt(1000.0d0)
      dtx = 0.001d0*dtx

      END
