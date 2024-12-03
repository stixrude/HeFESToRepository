      SUBROUTINE CALSCT(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, iter,
     $ betamx, bettol, btxtol, tempk, tau, press, rhosv, rhosl,
     $ deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $ dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $ ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $ dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $ ktxl, ztxl)

c     This routine calculates the saturation properties as a function
c     of specified temperature. This is done using Newton-Raphson
c     iteration to refine values of pressure, vapor density, and
c     liquid density, starting with results obtained using approximate
c     equations included by Wagner and Pruss (2002) in their description
c     of the IAPWS-95 model.

      IMPLICIT NONE

c     CALLING SEQUENCE VARIABLES.

      LOGICAL qerr, qfail, qprnt1, qprnt2, qprnt3, qsilent, qwrphi

      INTEGER nttyo, noutpt

      INTEGER iter

      REAL(8) pcr, rcnstw, rhocr, tempk, tau, tcr

      REAL(8) betamx, bettol, btxtol

      REAL(8) press, rhosv, rhosl, deltsv, deltsl

      REAL(8) pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv,
     $ bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv, ktxv, ztxv

      REAL(8) pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl,
     $ bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl, ktxl, ztxl

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     LOCAL VARIABLES.

      INTEGER k_par

      PARAMETER (k_par = 3)

      LOGICAL qrestl, qxiter, qsilnx

      INTEGER ILNOBL

      INTEGER i, itermx, irelax, j, j2, kdim, kmax, icutv,
     $ icutl, icutvq, icutlq, ncut

      CHARACTER(LEN=8) ux8, unam8

      REAL(8) delta, rho, rhosv0, rhosl0, psat, psat0, psatt

      REAL(8) delta0, dltsvq, dltslq, dltsv0, dltsl0, pxm, px0, pdiff

      REAL(8) aamatr(k_par,k_par), alpha(k_par), deltas(k_par),
     $ beta(k_par), betmx0

      REAL(8) xx, xxx

      REAL(8) arelax, dix, dltx

      REAL(8) bettl1, bettl2, bettl3

      DATA qsilnx / .TRUE. /

      DATA bettl1 / 1.0d-8 /
      DATA bettl2 / 1.0d-7 /
      DATA bettl3 / 1.0d-6 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (.NOT.qsilent) THEN

        WRITE (nttyo,1010) tempk, tau
        WRITE (noutpt,1010) tempk, tau

 1010   FORMAT(/3x,'CALSCT: Temp(K) = ',f9.4,13x,'tau    = ',e16.9,/)

      ENDIF

      kmax = k_par
      kdim = kmax

      qerr = .FALSE.
      qfail = .FALSE.
      qrestl = .FALSE.
      qxiter = .FALSE.
      itermx = 35
      btxtol = bettol

      IF (tempk.LE.298.15d0) THEN
        qrestl = .TRUE.
        btxtol = bettl1
        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1020) tempk, bettl1
          WRITE (noutpt,1020) tempk, bettl1

 1020     FORMAT(6x,'Temp(K) = ',f10.5,', LIES BETWEEN THE TRIPLE',
     $    ' POINT and 298.15K.',/6x,'APPLYING LOOSER CONVERGENCE',
     $    ' TOLERANCE bettl1 = ',1pe10.3,'.',/)

        ENDIF
      ENDIF

      IF (tempk.GT.647.090d0 .AND. tempk.LT.tcr) THEN
        IF (.NOT. qsilent) THEN

          WRITE (nttyo,1040)
          WRITE (noutpt,1040)

 1040     FORMAT(6x,'AS PRESENTLY TUNED, SATURATION CURVE ITERATION',
     $    /6x,'WILL NOT CONVERGE NORMALLY FOR TEMPERATURES GREATER',
     $    /6x,'THAN 647.090 K AND LESS THAN THE CRITICAL TEMPERATURE',
     $    /6x,'OF 647.096 K. A MAXIMUM OF FIVE ITERATIONS WILL BE',
     $    ' DONE.',/)

        ENDIF
        qxiter = .TRUE.
      ENDIF

      IF (tempk.GT.647.082d0) THEN
        qrestl = .TRUE.
        btxtol = bettl2
        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1050) tempk, bettl2
          WRITE (noutpt,1050) tempk, bettl2

 1050     FORMAT(6x,'Temp(K) = ',f10.5,', VERY CLOSE TO THE CRITICAL',
     $    ' TEMPERATURE.',/6x,'APPLYING LOOSER CONVERGENCE TOLERANCE',
     $    ' bettl1 = ',1pe10.3,'.',/)

        ENDIF
      ENDIF

      arelax = 1.0d0

      rhosv = 0.0d0
      rhosl = 0.0d0
      press = 0.0d0

c     Calculate approximate saturation pressure and
c     corresponding densities of liquid and vapor.
c     These results are not those of the IAPWS-95 model
c     itself, but can serve as starting estimates.

      CALL APXSCT(noutpt, nttyo, qprnt1, qsilent, rhocr, rhosv,
     $ rhosl, tempk, tcr, pcr, psat, psatt, qerr)

      IF (qerr) THEN

c       Error, specified temperature is not in the allowed range.

        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1210) tempk
          WRITE (noutpt,1210) tempk

 1210     FORMAT(/6x,'ERROR (CALSCT): Temp = ',f9.4,'K. This is out',
     $    /9x,'of range for the saturation curve.')

        ENDIF

        qfail = .TRUE.
        GO TO 999
      ENDIF

      icutv = 0
      icutl = 0
      icutvq = 0
      icutlq = 0

      deltsv = rhosv/rhocr
      deltsl = rhosl/rhocr

c     Save the values from the approximation.

      dltsvq = deltsv
      dltslq = deltsl

  120 CONTINUE

      IF (qprnt1 .AND. .NOT.qsilent) THEN

        WRITE (nttyo,1240) tempk, psat, tau, deltsv, rhosv,
     $  deltsl, rhosl
        WRITE (noutpt,1240) tempk, psat, tau, deltsv, rhosv,
     $  deltsl, rhosl

 1240   FORMAT(/6x,'CALSCT: Starting values for Temp(K) = ',f9.4,
     $  //9x,'psat   = ',e16.9,' MPa',6x,'tau    = ',e16.9,
     $  /9x,'deltsv = ',e16.9,10x,'rhosv  = ',e16.9,
     $  /9x,'deltsl = ',e16.9,10x,'rhosl  = ',e16.9,/)

      ENDIF

      press = psat
      rhosv0 = rhosv
      rhosl0 = rhosl

      iter = 0
      irelax = 0
      betmx0 = 1.0d+100

      psat0 = psat
      dltsv0 = deltsv
      dltsl0 = deltsl

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Below is the return point to refine the saturation
c     curve properties.

  250 CONTINUE

c     First calculate the vapor properties by calling EVAI95
c     with the vapor delta.

      unam8 = 'Vapor'
      j2 = ILNOBL(unam8)

      IF (qprnt2) THEN

        WRITE (noutpt,1250) unam8(1:j2), deltsv, rhosv

 1250   FORMAT(6x,a,/9x,'delta  = ',e16.9,10x,'rho   = ',e16.9,/)

      ENDIF

c     Calling sequence substitutions:

c       deltsv for delta
c       rhosv sfor rho
c       pxv for px
c       axv for ax
c       uxv for ux
c       sxv for sx
c       hxv for hx
c       cvxv for cvx
c       cpxv for cpx
c       wxv for wx
c       muxv for mux
c       dtxv for dtx
c       bsxv for bsx
c       gxv for gx
c       vxv for vx
c       pdxv for pdx
c       ptxv for ptx
c       adxv for adx
c       gdxv for gdx
c       avxv for avx
c       ktxv for ktx
c       ztxv for ztx
c       qsilnx for qsilent

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, deltsv, rhosv, pxv, axv, uxv, sxv, hxv,
     $ cvxv, cpxv, wxv, muxv, dtxv, bsxv, gxv, vxv, pdxv, ptxv, adxv,
     $ gdxv, avxv, ktxv, ztxv, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilnx)

      IF (iter .LE. 0) THEN
        IF (.NOT.qsilent .AND.qprnt1) THEN

          WRITE (nttyo,1270) pxv, pdxv
          WRITE (noutpt,1270) pxv, pdxv

 1270     FORMAT(/9x,'pxv    = ',e16.9,' MPa',6x,'pdxv   = ',e16.9)

        ENDIF
      ENDIF

c     Now calculate the liquid properties by calling EVAI95
c     with the liquid delta.

      unam8 = 'Liquid'
      j2 = ILNOBL(unam8)

      IF (qprnt2) THEN
        WRITE (noutpt,1250) unam8(1:j2), deltsl, rhosl
      ENDIF

c     Calling sequence substitutions:

c       deltsl for delta
c       rhosl sfor rho
c       pxl for px
c       axl for ax
c       uxl for ux
c       sxl for sx
c       hxl for hx
c       cvxl for cvx
c       cpxl for cpx
c       wxl for wx
c       muxl for mux
c       dtxl for dtx
c       bsxl for bsx
c       gxl for gx
c       vxl for vx
c       pdxl for pdx
c       ptxl for ptx
c       adxl for adx
c       gdxl for gdx
c       avxl for avx
c       ktxl for ktx
c       ztxl for ztx
c       qsilnx for qsilent

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, deltsl, rhosl, pxl, axl, uxl, sxl, hxl,
     $ cvxl, cpxl, wxl, muxl, dtxl, bsxl, gxl, vxl, pdxl, ptxl, adxl,
     $ gdxl, avxl, ktxl, ztxl, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilnx)

      IF (iter .LE. 0) THEN
        IF (.NOT.qsilent .AND.qprnt1) THEN

          WRITE (nttyo,1292) pxl, pdxl
          WRITE (noutpt,1292) pxl, pdxl

 1292     FORMAT(9x,'pxl    = ',e16.9,' MPa',6x,'pdxl   = ',e16.9)

        ENDIF
      ENDIF

      IF (qprnt2) THEN
        xx = (axv - axl)

        WRITE (noutpt,2420) axv, axl, xx
        WRITE (nttyo,2420) axv, axl, xx

 2420   FORMAT(/9x,'axv    = ',e16.9,' kJ/kg',4x,'axl   = ',e16.9,
     $  /9x,'adif   = ',e16.9,' kJ/kg')

        WRITE (noutpt,2430) adxv, adxl
        WRITE (nttyo,2430) adxv, adxl

 2430   FORMAT(9x,'adxv   = ',e16.9,' kJ/kg',4x,'adxl  = ',e16.9)

      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c       The pdx for vapor cannot be negative.
c       Under-relax to prevent this.

        IF (pdxv. lt. 0.0d0) THEN

          IF (iter .LE. 0) THEN
            IF (icutv .ge. 30) GO TO 410
            icutv = icutv + 1
            rhosv = 0.995d0*rhosv
            deltsv = rhosv/rhocr
            IF (qprnt1) THEN

              WRITE (nttyo,1140)
              WRITE (noutpt,1140)

 1140         FORMAT(/3x,'CALSCT: DECREASING INITIAL VAPOR DENSITY',
     $        ' TO AVOID NEGATIVE PDX',/)

            ENDIF
            GO TO 120
          ELSE
            IF (icutv .ge. 30) GO TO 410
            icutv = icutv + 1
            IF (qprnt1) THEN

              WRITE (nttyo,1150)
              WRITE (noutpt,1150)

 1150         FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID NEGATIVE',
     $        ' PDX FOR VAPOR',/)

            ENDIF
            arelax = 0.25d0
            GO TO 150
          ENDIF
        ENDIF
        icutv = 0

c       The pdx for liquid cannot be negative.
c       Under-relax to prevent this.

        IF (pdxl. lt. 0.0d0) THEN

          IF (iter .LE. 0) THEN
            IF (icutl .ge. 30) GO TO 410
            icutl = icutl + 1
            rhosl = 1.001d0*rhosl
            deltsl = rhosl/rhocr
            IF (qprnt1) THEN

              WRITE (nttyo,1160)
              WRITE (noutpt,1160)

 1160         FORMAT(/3x,'CALSCT: INCREASING INITIAL LIQUID DENSITY',
     $        ' TO AVOID NEGATIVE PDX',/)

            ENDIF
            GO TO 120
          ELSE
            IF (icutl .ge. 30) GO TO 410
            icutl = icutl + 1
            IF (qprnt1) THEN

              WRITE (nttyo,1170)
              WRITE (noutpt,1170)

 1170         FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID NEGATIVE',
     $        ' PDX FOR LIQUID',/)

            ENDIF
            arelax = 0.25d0
            GO TO 150
          ENDIF
        ENDIF
        icutl = 0

c       The revised delta for vapor cannot be less than a good
c       fraction of the value obtained from the intial approximation.

        IF (deltsv. lt. 0.90d0*dltsvq) THEN

          IF (icutvq .ge. 30) GO TO 410
          icutvq = icutvq + 1
          arelax = 0.25d0

          IF (qprnt1) THEN

            WRITE (nttyo,1180)
            WRITE (noutpt,1180)

 1180       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID VAPOR DELTA',
     $      /3x,'TOO MUCH LESS THAN THE APPROXIMATION VALUE.',/)

          ENDIF

          GO TO 150

        ENDIF
        icutvq = 0

c       The revised delta for liquid cannot be greater than a small
c       fraction of the value obtained from the intial approximation.

        IF (deltsl. gt. 1.10d0*dltslq) THEN

          IF (icutlq .ge. 30) GO TO 410
          icutlq = icutlq + 1
          IF (qprnt1) THEN

            WRITE (nttyo,1190)
            WRITE (noutpt,1190)

 1190       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID LIQUID DELTA',
     $      /3x,'TOO MUCH GREATER THAN THE APPROXIMATION VALUE.',/)

          ENDIF
          arelax = 0.25d0
          GO TO 150
        ENDIF
        icutlq = 0

        GO TO 170

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  150 CONTINUE

      IF (qprnt2) THEN
        IF (arelax .NE. 1.0d0) THEN
          WRITE (noutpt,1530) arelax
        ENDIF
      ENDIF

      xxx = 0.0d0

      DO j = 1,3
        deltas(j) = arelax*deltas(j)
        xx = DABS(deltas(j))
        xxx = DMAX1(xxx,xx)
      ENDDO

      psat = psat0 + deltas(3)
      deltsv = dltsv0 + deltas(1)
      deltsl = dltsl0 + deltas(2)
      rhosv = rhocr*deltsv
      rhosl = rhocr*deltsl
      GO TO 250

  170 CONTINUE

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Have obtained valid (outside the unstable zone) vapor and
c     liquid properties for the current iteration. Improve the
c     calculated saturation properties by solving three equations
c     in three unknowns. The equations are all in terms of pressure.
c     The unknowns to be found are psat, deltsv, and deltsl.

c     Calculate the Maxwell criterion pressure
c     (Gibbs energy equality expressed through the
c     Helmholtz energies and the pressure)

      dix = (1.0d0/deltsl) - (1.0d0/deltsv)
      dltx = deltsl - deltsv
      IF (.NOT.qsilent .AND. qprnt2) THEN

        WRITE (noutpt,1295) dix, dltx

 1295   FORMAT(9x,'dix    = ',e16.9,'    ',6x,'dltx   = ',e16.9)

      ENDIF

      IF (DABS(dix) .GT. 1.0d-15) THEN
        IF (DABS(axv - axl) .GT. 1.0d-15) THEN

c         Normal calculation, result in kPa.

          pxm = rhocr*( axv - axl )/dix

c         Convert from kPa to MPa.

          pxm = 0.001d0*pxm
        ELSE

c         There is no difference in the Helmholtz
c         energies of the vapor and the liquid.

          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1300)
            WRITE (noutpt,1300)

 1300       FORMAT(/3x,'CALSCT: THE TWO HELMHOLTZ ENERGIES ARE EQUAL.',
     $      /3x,'THE MAXWELL EQUATION CANNOT BE EVALUATED.')

          ENDIF

          GO TO 410
        ENDIF
      ELSE

c       Exception intended for the critical point.

        IF (DABS(tempk - tcr).LE.1.0d-10 .AND.
     $    DABS(deltsv - 1.0d0).LE.1.0d-10 .AND.
     $    DABS(deltsl - 1.0d0).LE.1.0d-10) THEN

c         Am at the critical point.

          pxm = pxv
        ELSE

c         Not at the critical point, but the vapor
c         and liquid densities have converged.

          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1310)
            WRITE (noutpt,1310)

 1310       FORMAT(/3x,'CALSCT: THE TWO DELTA VALUES ARE EQUAL.',
     $      /3x,'THE MAXWELL EQUATION CANNOT BE EVALUATED.')

          ENDIF

          GO TO 410
        ENDIF
      ENDIF

      IF (iter .LE. 0) THEN
        IF (.NOT.qsilent .AND. qprnt1) THEN

          WRITE (nttyo,1320) pxm
          WRITE (noutpt,1320) pxm

 1320     FORMAT(/6x,'Maxwell',/9x,'pxm    = ',e16.9,' MPa',/)

        ENDIF
      ENDIF

c     Calculate residual functions.

      alpha(1) = pxm - psat
      alpha(2) = pxv - psat
      alpha(3) = pxl - psat

      beta(1) = DABS(alpha(1)/psat)
      beta(2) = DABS(alpha(2)/psat)
      beta(3) = DABS(alpha(3)/psat)

      betamx = DMAX1(beta(1), beta(2), beta(3))

      IF (qprnt2) THEN

        WRITE (noutpt,1330) (alpha(i), i = 1,kdim)

 1330   FORMAT(/6x,'Residuals',
     $  /9x,'alpha(1) = ',e16.9,1x,'(Maxwell)',
     $  /9x,'alpha(2) = ',e16.9,1x,'(Vapor)',
     $  /9x,'alpha(3) = ',e16.9,1x,'(Liquid)')

        WRITE (noutpt,1332) (beta(i), i = 1,kdim)

 1332   FORMAT(/6x,'Relative Residuals',
     $  /10x,'beta(1) = ',e16.9,1x,'(Maxwell)',
     $  /10x,'beta(2) = ',e16.9,1x,'(Vapor)',
     $  /10x,'beta(3) = ',e16.9,1x,'(Liquid)')

      ENDIF

      IF (.NOT.qsilent) THEN

        WRITE (noutpt,1335) iter, psat, betamx
        WRITE (nttyo,1335) iter, psat, betamx

 1335   FORMAT(6x,'CALSCT: iter= ',i3,', psat= ',1pe16.9,', betamx= ',
     $  e12.5)

      ENDIF

c     Note: using a convergence tolerance below 1.0d-11
c     may lead to non-convergence due to the limitations
c     of 64-bit arithmetic.

      IF (betamx .LE. btxtol) THEN

c       Iteration has converged.

        press = psat

        IF (.NOT.qsilent .AND. qprnt1) THEN

          WRITE (nttyo,1360) press, betamx, btxtol
          WRITE (noutpt,1360) press, betamx, btxtol

 1360     FORMAT(/3x,'CALSCT: ITERATION CONVERGED TO press= ',1pe16.9,
     $    ' (MPa).',/3x,'MAX NORM betamx = ',e10.3,
     $    ', TOLERANCE btxtol = ',1pe10.3)

        ENDIF

        IF (deltsv .LE. 1.0d-12) THEN

          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1365)
            WRITE (noutpt,1365)

 1365       FORMAT(/3x,'CALSCT: THE SOLUTION IS A FALSE ONE. THE VAPOR',
     $      /3x,'DENSITY AND PRESSURE ARE NEAR-ZERO.',/)

          ENDIF

          GO TO 410
        ENDIF

        GO TO 420

      ELSEIF (iter .GE. itermx) THEN

c       Have done the maximum number of iterations.

        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1370)
          WRITE (noutpt,1370)

 1370     FORMAT(/3x,'CALSCT: HAVE DONE THE MAXIMUM NUMBER OF',
     $    ' ITERATIONS.')

        ENDIF

        GO TO 410

      ELSEIF (qxiter .AND. iter.GE.5) THEN

c       Have done the maximum number of iterations.
c       Report calculated results with warning.

        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1375) betamx, btxtol
          WRITE (noutpt,1375) betamx, btxtol

 1375     FORMAT(/3x,'CALSCT: CALCULATED RESULTS ARE APPROXIMATE.',
     $    /3x,'MAX NORM betamx = ',1pe10.3,', TOLERANCE btxtol = ',
     $    e10.3)

        ENDIF

        GO TO 420

      ELSE

c       Make corrections and do another iteration.

        iter = iter + 1

c       The Jacobian matrix J here is aamatr(kdim,kdim).

        DO i = 1,3
          DO j = 1,3
            aamatr(i,j) = 0.30
          ENDDO
        ENDDO

        aamatr(1,1) = -(1.0d0/(dix*deltsv**2)) + (adxv/(axv - axl))
        aamatr(1,2) =  (1.0d0/(dix*deltsl**2)) - (adxl/(axv - axl))
        aamatr(1,1) = pxm*aamatr(1,1)
        aamatr(1,2) = pxm*aamatr(1,2)
        aamatr(1,3) = -1.0d0

        aamatr(2,1) = pdxv
        aamatr(2,2) = 0.0d0
        aamatr(2,3) = -1.0d0

        aamatr(3,1) = 0.0d0
        aamatr(3,2) = pdxl
        aamatr(3,3) = -1.0d0

        IF (qprnt3) THEN

          WRITE (noutpt,1410)

 1410     FORMAT(/3x,'Starting matrix')

          WRITE (noutpt,1420) ((aamatr(i,j), j = 1,kdim),
     $    alpha(i), i = 1,kdim)

 1420     FORMAT(/6x,'Jacobian and residual',
     $    /9x,1pe12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $    /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $    /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5)

        ENDIF

c       Since this matrix is only 3x3, solve using Gaussian
c       elimination with partial pivoting. There is no need
c       to use the "Jordan" part of the Gauss-Jordan algorithm
c       for a matrix of this size.

c       The simultaneous equations have the form:

c         aamatr(i,1)*deltas(1) + aamatr(i,2)*deltas(2)
c           aamatr(i,3)*deltas(3) = -alpha(i), i = 1,3

        CALL GAUSSE(noutpt, nttyo, qerr, qprnt3, kdim,
     $  kmax, aamatr, alpha, deltas)

        IF (qerr) THEN
          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1430)
            WRITE (noutpt,1430)

 1430       FORMAT(/6x,'ERROR (CALSCT): Gaussian elimination failed.')

          ENDIF

          GO TO 410
        ENDIF

c       Reverse the sign of deltas to be consistent
c       with definition of correction.

        DO j = 1,kdim
          deltas(j) = -deltas(j)
        ENDDO

c       IF (qprnt2) THEN
c
c         WRITE (noutpt,1510)
c
c1510     FORMAT(/6x,'Corrections')
c
c         DO j = 1,kdim
c
c           WRITE (noutpt,1520) j, deltas(j)
c
c1520       FORMAT(13x,i2,3x,1pe16.9)
c
c         ENDDO
c       ENDIF

        ncut = 0
        arelax = 1.0d0

        IF (qprnt2) THEN
          IF (arelax .NE. 1.0d0) THEN

            WRITE (noutpt,1530) arelax

 1530       FORMAT(/9x,'arelax = ',e16.9)

          ENDIF
        ENDIF

c       Save current values.

        psat0 = psat
        dltsv0 = deltsv
        dltsl0 = deltsl

c       Make corrections.

  300   CONTINUE
        IF (qprnt2) THEN

          WRITE (noutpt,1510)

 1510     FORMAT(/6x,'Corrections')

          DO j = 1,kdim

            WRITE (noutpt,1520) j, deltas(j)

 1520       FORMAT(13x,i2,3x,1pe16.9)

          ENDDO
        ENDIF

        psat = psat0 + deltas(3)
        deltsv = dltsv0 + deltas(1)
        deltsl = dltsl0 + deltas(2)

c       The delta for liquid cannot be less than the delta for vapor.
c       Under-relax to prevent this.

        IF (deltsl. lt. deltsv) THEN
          IF (ncut .ge. 20) GO TO 410
          ncut = ncut + 1
          IF (qprnt1) THEN

            WRITE (nttyo,1540)
            WRITE (noutpt,1540)

 1540       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID VAPOR-LIQUID',
     $      ' DENSITY INVERSION',/)

          ENDIF
          arelax = 0.25d0
          GO TO 320
        ENDIF

c       Corrected delta values must be positive to avoid
c       a singularity in the equation-of-state model equations.

        IF (deltsv.LE.0.0d0 .OR. deltsl.LE.0.0d0) THEN
          IF (ncut .ge. 20) GO TO 410
          ncut = ncut + 1
          IF (qprnt1) THEN

            WRITE (nttyo,1580)
            WRITE (noutpt,1580)

 1580       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID NEGATIVE',
     $      ' DELTA.',/)

          ENDIF
          arelax = 0.25d0
          GO TO 320
        ENDIF

        GO TO 330

  320   CONTINUE
        xxx = 0.0d0
        DO j = 1,3
          deltas(j) = arelax*deltas(j)
          xx = DABS(deltas(j))
          xxx = DMAX1(xxx,xx)
        ENDDO
        GO TO 300

  330   rhosv = rhocr*deltsv
        rhosl = rhocr*deltsl

c       Test for vapor-liquid density inversion.
c       That is, is the vapor more dense than the liquid?

        IF (qprnt2) THEN

          WRITE (noutpt,1600) psat, deltsv, deltsl

 1600     FORMAT(/6x,'Corrected variables',
     $    /9x,'psat   = ',e16.9,' MPa',
     $    /9x,'deltsv = ',e16.9,/9x,'deltsl = ',e16.9)

        ENDIF

        GO TO 250

      ENDIF

      GO TO 420

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  410 CONTINUE

c     Iteration failed.

      qfail = .TRUE.
      IF (.NOT.qsilent) THEN

        WRITE (nttyo,1690)
        WRITE (noutpt,1690)

 1690   FORMAT(/3x,'CALSCT: ITERATION FAILED.',/)

      ENDIF

      GO TO 999

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  420 CONTINUE

c     Calculations are complete for the last line read
c     from the input file. Write the results.

      IF (.NOT.qsilent) THEN

        WRITE (nttyo,2400)
        WRITE (noutpt,2400)

 2400   FORMAT(/3x,'CALSCT key results:',//8x,'Temp(K)',10x,
     $  'press(MPa)')

        WRITE (nttyo,2410) tempk,psat
        WRITE (noutpt,2410) tempk,psat

 2410   FORMAT(6x,f9.4,7x,e16.9)

        WRITE (noutpt,2415) tempk, psat, tau, deltsv, rhosv,
     $   deltsl, rhosl

 2415   FORMAT(/6x,'Temp(K) = ',f9.4
     $  //9x,'psat   = ',e16.9,' MPa',6x,'tau    = ',e16.9,
     $  /9x,'deltsv = ',e16.9,10x,'rhosv  = ',e16.9,
     $  /9x,'deltsl = ',e16.9,10x,'rhosl  = ',e16.9)

      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 999  CONTINUE

      END
