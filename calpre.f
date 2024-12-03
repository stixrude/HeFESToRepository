      SUBROUTINE CALPRE(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qwrphi, qrhog, rhog, pcr, rcnstw, rhocr, tcr, iter,
     $ betamx, bettol, rhotol, btxtol, udescr, press, tempk, tau,
     $ delta, rho, psat, px, ux, sx, hx, cvx, cpx, wx, mux, dtx,
     $ bsx, ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx)

c     Find the thermodynamic properties of water at given temperature
c     and pressure. The problem reduces to finding the value of reduced
c     density (delta) that is consistent with the desired pressure.
c     The Newton-Raphson method is employed. Under-relaxation techniques
c     are required to assure that delta or other requisite parameters
c     do not take on out-of-bounds values in the iteration process.
c     Small negative values of calculated pressure are okay. Zero or
c     negative values for calculated "pdx" (pressure derivative with
c     respect to delta) imply the unstable zone and must be avoided.

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qerr, qfail, qprnt1, qprnt2, qprnt3, qwrphi, qrhog

      CHARACTER(LEN=24) udescr

      INTEGER nttyo, noutpt

      INTEGER iter

      REAL(8) press, tempk, tau, delta, rho, rhog, psat

      REAL(8) pcr, rcnstw, rhocr, tcr

      REAL(8) betamx, bettol, btxtol, rhotol

      REAL(8) px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      LOGICAL qsilent

      INTEGER ILNOBL

      INTEGER icdlta, itermx, irelax, irlxmx, j2

      CHARACTER(LEN=8) ux8

      REAL(8) rhoidg, rhosup, rhocpa, rhosv, rhosl, rxx, psatt

      REAL(8) arelax, delta0, deltas, delts0, deltol, dltamx, dltsx,
     $ pdiff, px0, rho0, sbetmx

      DATA qsilent / .TRUE. /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      itermx = 60
      IF (qrhog) itermx = 80
      btxtol = bettol
      deltol = 1.0d-12
      icdlta = 0
      qerr = .FALSE.
      qfail = .FALSE.

      if (.not. qsilent) WRITE (nttyo,1050) tempk, tau, press
      if (.not. qsilent) WRITE (noutpt,1050) tempk, tau, press

 1050 FORMAT(/3x,'CALPRE: Temp(K) = ',f9.4,13x,'tau    = ',e16.9,
     $ /11x,'press   = ',e16.9,' MPa',/)

      IF (qrhog) THEN

        if (.not. qsilent) WRITE (nttyo,1060) rhog
        if (.not. qsilent) WRITE (noutpt,1060) rhog

 1060   FORMAT(11x,'rhog   = ',e16.9,' kg/m3 (starting value)',/)

      ENDIF

      iter = 0
      rho = 0.0d0

c     Obtain a description (udescr) of the H2O fluid.
c     Ignore the qerr error flag.

      CALL FDESCR(nttyo, noutpt, qprnt1, qprnt2, qprnt3,
     $ qwrphi, pcr, rcnstw, rhocr, tcr, bettol, btxtol, psat,
     $ rho, rhosv, rhosl, rhotol, tempk, tau, press, udescr)

      IF (udescr(1:8) .EQ. 'vapor   ') THEN

c       Vapor: assume ideal gas behavior.

c       Ideal gas.
        rho = 1000.0d0*press/(tempk*rcnstw)

c       Alternate: ideal gas correction to saturation
c       density (slightly different).
c       rho = (press/psat)*rhosv

      ELSEIF (udescr(1:8) .EQ. 'liquid  ') THEN

c       Liquid: use a liquid-like density.

c       The liquid density on the saturation curve.
        rho = rhosl

c       Something slightly greater than the saturation
c       density.
c       rho = 1.0002d0*rhosl

      ELSEIF (udescr(1:8) .EQ. 'compress') THEN

c       Estimate the density of compressed liquid.

c       Ideal gas correction to critical point density.
c       This produces one of two false solutions at 500K
c       and 25 MPa. Do not use this.
c       rho = (press/pcr)*(tcr/tempk)*rhocr

c       Twice the ideal gas correction to the critical
c       point density.
c       rho = 2.0d0*(press/pcr)*(tcr/tempk)*rhocr

c       The saturated liquid density is a minimum value.

        rho = 1.10d0*rhosl

c       Close to the upper limit for this field
c       (T near the triple point, 1000 MPa).
c       For higher pressure, a higher value might
c       be needed.
c       rho = 1250.0d0

        rho = DMAX1(rho, rhosl)
        rho = DMIN1(rho, 1400.0d0)

      ELSEIF (udescr(1:8).EQ.'supercri') THEN

c       Estimate the density of supercritical fluid.

c       Ideal gas.
c       rho = 1000.0d0*press/(tempk*rcnstw)

c       SUPCRT92 estimate, about 15% higher than ideal gas.
c       rho = 2500.0d0*press/tempk

c       Ideal gas correction to critical point density.
c       rho = (press/pcr)*(tcr/tempk)*rhocr

c       Twice the ideal gas correction to the critical
c       point density.
        rho = 2.0d0*(press/pcr)*(tcr/tempk)*rhocr

c       Close to the upper limit for this P, T field.
c       (T near the critical point, 1000 MPa).
c       Calculated value of 1068.7 kg/m3 is rounded up.
c       For higher pressure, a higher value might
c       be needed.
c       rho = 1100.0d0

        rho = DMIN1(rho, 1100.0d0)

      ELSEIF (udescr(1:8).EQ.'hot vapo') THEN

c       Estimate the density of hot vapor.

c       Ideal gas.
        rhoidg = 1000.0d0*press/(tempk*rcnstw)

c       SUPCRT92 estimate, about 15% higher than ideal gas.
        rhosup = 2500.0d0*press/tempk

c       Ideal gas correction to critical point density.
        rhocpa = (press/pcr)*(tcr/tempk)*rhocr

c       The upper limit for this field, the critical pressure
c       (rhocr), 22.064 MPa.
c       rho = rhocr

        IF (press .LE. 1.00d0) THEN
          rho = rhoidg
        ELSEIF (press .LE. 18.0d0) THEN
          rho = rhosup
        ELSE
          rho = rhocpa
        ENDIF

        rho = DMIN1(rho, rhocr)

      ELSE

c       The H2O fluid type could not be determined.

        if (.not. qsilent) WRITE (nttyo,1070)
        if (.not. qsilent) WRITE (noutpt,1070)

 1070   FORMAT(/6x,'WARNING (CALPRE): The H2O fluid type could not be',
     $  ' determined.',/9x,'A good starting estimate of density',
     $  ' could not be established.')

c       Will try four times the critical density.

        rho = 4*rhocr

      ENDIF

      IF (qrhog) THEN

c       Use the user-specified starting value for rho.

        rho = rhog
      ENDIF

      delta = rho/rhocr

      irelax = 0
      irlxmx = 20
      arelax = 0.25d0
      sbetmx = 0.0d0
      betamx = 0.0d0
      j2 = ILNOBL(udescr)

      if (.not. qsilent) WRITE (nttyo,1080) udescr(1:j2)
      if (.not. qsilent) WRITE (noutpt,1080) udescr(1:j2)

 1080 FORMAT(6x,'This appears to be ',a,'.',/)

c     Below is the return point for iterating on delta to
c     obtain a desired pressure.

  250 CONTINUE

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, delta, rho, px, ax, ux, sx, hx,
     $ cvx, cpx, wx, mux, dtx, bsx, gx, vx, pdx, ptx, adx,
     $ gdx, avx, ktx, ztx, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilent)

c     Test to see if the pressure has converged to the desired value.

      pdiff = px - press
      sbetmx = pdiff/press
      betamx = DABS(sbetmx)

      if (.not. qsilent) WRITE (noutpt,1250) iter, px, sbetmx
      if (.not. qsilent) WRITE (nttyo,1250) iter, px, sbetmx

 1250 FORMAT(6x,'CALPRE: iter= ',i3,', px= ',1pe16.9,', sbetmx= ',
     $ e12.5)

      IF (qprnt1) THEN

        if (.not. qsilent) WRITE (noutpt,1260) delta, rho
        if (.not. qsilent) WRITE (nttyo,1260) delta, rho

 1260   FORMAT(/9x,'delta  = ',e16.9,10x,'rho    = ',e16.9)

        if (.not. qsilent) WRITE (noutpt,1290) px, pdx
        if (.not. qsilent) WRITE (nttyo,1290) px, pdx

 1290   FORMAT(9x,'px     = ',e16.9,' MPa',6x,'pdx    = ',e16.9)

      ENDIF

      IF (betamx .LE. btxtol) THEN

c       Have converged.

        IF (qprnt1) THEN

          if (.not. qsilent) WRITE (nttyo,1670) betamx, btxtol
          if (.not. qsilent) WRITE (noutpt,1670) betamx, btxtol

 1670     FORMAT(/3x,'CALPRE: ITERATION CONVERGED.',
     $    /3x,'MAX NORM betamx =',1pe10.3,', TOLERANCE btxtol = ',
     $    e10.3,/)

        ENDIF

        GO TO 420
      ENDIF

      IF (iter .GE. itermx) THEN

c       Have done the maximum number of iterations.

        if (.not. qsilent) WRITE (nttyo,1370)
        if (.not. qsilent) WRITE (noutpt,1370)

 1370   FORMAT(/3x,'CALPRE: HAVE DONE THE MAXIMUM NUMBER OF',
     $  ' ITERATIONS.')

        GO TO 410
      ENDIF

c     Estimate a new value for delta (and rho).

      iter = iter + 1
      irelax = 0
      delta0 = delta
      rho0 = rho
      px0 = px
      deltas = -pdiff/pdx
      delts0 = deltas
      dltamx = DABS(deltas)

      IF (qprnt1) THEN

        if (.not. qsilent) WRITE (noutpt,1372) dltamx
        if (.not. qsilent) WRITE (nttyo,1372) dltamx

 1372   FORMAT(9x,'dltamx = ',e16.9,/)

      ENDIF

      IF (iter .GE. 5) THEN
        IF (dltamx .LE. deltol) icdlta = icdlta + 1
      ENDIF

      IF (icdlta .GE. 5) THEN

        IF (betamx .LE. 1.0d-6) THEN

c         Have pseudo-converged.

          if (.not. qsilent) WRITE (nttyo,1375) betamx, btxtol, dltamx, deltol
          if (.not. qsilent) WRITE (noutpt,1375) betamx, btxtol, dltamx, deltol

 1375     FORMAT(/3x,'CALPRE: ITERATION PSEUDO-CONVERGED.',
     $    /6x,'MAX NORM betamx = ',1pe10.3,', TOLERANCE btxtol = ',
     $    e10.3,/6x,'MAX NORM dltamx = ',e10.3,
     $    ', TOLERANCE deltol = ',1pe10.3)

          GO TO 420

        ENDIF

      ENDIF

      IF (icdlta .GE. 10) THEN

c       Iteration is failing to result in any improvement.

        if (.not. qsilent) WRITE (nttyo,1377) betamx, btxtol, dltamx, deltol
        if (.not. qsilent) WRITE (noutpt,1377) betamx, btxtol, dltamx, deltol

 1377   FORMAT(/3x,'CALPRE: ITERATION IS NOT LEADING TO IMPROVEMENT.',
     $  /6x,'MAX NORM betamx = ',1pe10.3,', TOLERANCE btxtol = ',
     $  e10.3,/6x,'MAX NORM dltamx = ',e10.3,
     $  ', TOLERANCE deltol = ',1pe10.3)

        GO TO 410

      ENDIF

      dltsx = 0.10d0*delta0
      IF (deltas .GT. dltsx) THEN

c       Clamp deltas to keep delta from large fractional increase.

        IF (qprnt1) THEN

          if (.not. qsilent) WRITE (nttyo,1277)
          if (.not. qsilent) WRITE (noutpt,1277)

 1277     FORMAT(/3x,'CLAMPING CHANGE IN DELTA TO AVOID LARGE',
     $    ' FRACTIONAL INCREASE.')

        ENDIF

        deltas = dltsx
      ENDIF

      dltsx = -0.40d0*delta0
      IF (deltas .LT. dltsx) THEN

c       Clamp deltas to keep delta from large fractional decrease.

        IF (qprnt1) THEN

          if (.not. qsilent) WRITE (nttyo,1278)
          if (.not. qsilent) WRITE (noutpt,1278)

 1278     FORMAT(/3x,'CLAMPING CHANGE IN DELTA TO AVOID LARGE',
     $    ' FRACTIONAL DECREASE.')

        ENDIF

        deltas = dltsx
      ENDIF

c     Newton-Raphson correction.

  220 delta = delta0 + deltas

      IF (delta .LT. 1.0d-15) then

c       Under-relax to keep delta from being too small. Only
c       positive values are physical.

        IF (irelax .GE. irlxmx) GO TO 405
          IF (qprnt1) THEN

          if (.not. qsilent) WRITE (nttyo,1330)
          if (.not. qsilent) WRITE (noutpt,1330)

 1330     FORMAT(/3x,'UNDER-RELAXING TO AVOID NEGATIVE DELTA.')

        ENDIF

        irelax = irelax + 1
        deltas = arelax*deltas
        GO TO 220
      ENDIF

      rho = rhocr*delta

      GO TO 250

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  405 CONTINUE

c     Have hit the limit for under-relaxation steps.

      if (.not. qsilent) WRITE (ux8,'(i8)') irlxmx
      j2 = ILNOBL(ux8)
      CALL LEJUST(ux8)

      if (.not. qsilent) WRITE (nttyo,1350) ux8(1:j2)
      if (.not. qsilent) WRITE (noutpt,1350) ux8(1:j2)

 1350 FORMAT(/3x,'CALPRE: HAVE HIT THE UNDER-RELAXATION LIMIT OF ',a,
     $ /3x,' TIMES PER ITERATION.')

c     Continue iteration.

      rho = rhocr*delta

      GO TO 250

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  410 CONTINUE

c     Iteration failed.

      qfail = .TRUE.

      WRITE (nttyo,1690)
      WRITE (noutpt,1690)

 1690 FORMAT(/3x,'CALPRE: ITERATION FAILED.',/)

      GO TO 999

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  420 CONTINUE

c     Calculations are complete for the last line read
c     from the input file. Write the results.

      if (.not. qsilent) WRITE (nttyo,2400)
      if (.not. qsilent) WRITE (noutpt,2400)

 2400 FORMAT(/8x,'Temp(K)',10x,'press(MPa)')

      if (.not. qsilent) WRITE (nttyo,2410) tempk,px
      if (.not. qsilent) WRITE (noutpt,2410) tempk,px

 2410 FORMAT(6x,f9.4,7x,e16.9)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  999 CONTINUE

      END
