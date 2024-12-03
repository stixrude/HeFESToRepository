      SUBROUTINE FDESCR(nttyo, noutpt, qprnt1, qprnt2, qprnt3,
     $ qwrphi, pcr, rcnstw, rhocr, tcr, bettol, btxtol, psat,
     $ rho, rhosv, rhosl, rhotol, tempk, tau, press, udescr)

c     Given the temperature and pressure, find the appropriate
c     description of the H2O fluid. A problem may occur if the
c     pressure is equal or nearly equal to the saturation pressure.
c     Here comparing the pressure with the saturation pressure
c     pressure may lead to the wrong description, as vapor and
c     liquid coexist at the saturation pressure. It then becomes
c     neccesary to compare the fluid density with the saturated
c     vapor and saturated liquid densities. If the density is
c     known (as following a CALDLT calculation), it will be used.
c     If it is not known (as in a CALPRE calculation), it will
c     not be used. In that case, the results obtained here will
c     determine the starting density estimate, thus in essence
c     choosing "vapor" or "liquid" for pressures close to the
c     saturation pressure.

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     CALLING SEQUENCE VARIABLES.

      LOGICAL qprnt1, qprnt2, qprnt3, qwrphi

      CHARACTER(LEN=24) udescr

      INTEGER nttyo, noutpt

      REAL(8) tempk, press

      REAL(8) pcr, rcnstw, rhocr, tau, tcr

      REAL(8) psat, rho, rhosv, rhosl, rhotol

      REAL(8) bettol, btxtol

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     LOCAL VARIABLES.

      LOGICAL qsilent, qerrx, qfailx

      INTEGER iter

      REAL(8) betamx, btxsav, deltsv, deltsl

      REAL(8) pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv,
     $ bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv, ktxv, ztxv

      REAL(8) pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl,
     $ bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl, ktxl, ztxl

      REAL(8) ptest, btest, rtestl, rtestv

      DATA qsilent / .TRUE. /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      udescr = 'unknown'
      qerrx = .FALSE.
      qfailx = .FALSE.
      psat = 0.0d0
      rhosv = 0.0d0
      rhosl = 0.0d0

      IF (tempk.LT.273.15d0) THEN

c       Note that the allowed temperature range has been extended a bit
c       on the low end to include 0C.

        WRITE (nttyo,1010)
        WRITE (noutpt,1010)

 1010   FORMAT(/6x,'WARNING (FDESCR): Temperature is too low to',
     $  ' assign',/9x,'a description to the H2O fluid.')

        GO TO 999
      ENDIF

      IF (tempk .LE. tcr) THEN

        btxsav = btxtol

c       Calculate the saturation curve properties.

c       Calling sequence substitutions:

c         psat for press
c         qerrx for qerr
c         qfailx for qfail

        CALL CALSCT(nttyo, noutpt, qerrx, qfailx, qprnt1, qprnt2,
     $  qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, iter,
     $  betamx, bettol, btxtol, tempk, tau, psat, rhosv, rhosl,
     $  deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $  dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $  ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $  dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $  ktxl, ztxl)

        IF (qerrx) THEN

          WRITE (nttyo,1020)
          WRITE (noutpt,1020)

 1020     FORMAT(/6x,'WARNING (FDESCR): The call to CALSCT to',
     $    ' calculate',/9x,'the saturation pressure failed.')

          rhosv = 0.0d0
          rhosl = 0.0d0
        ENDIF

c       Reset the convergence tolerance.

        btxtol = btxsav

        IF (press .GT. pcr) THEN
          udescr = 'compressed liquid'

          IF (qerrx) THEN

            WRITE (nttyo,1030)
            WRITE (noutpt,1030)

 1030       FORMAT(/6x,'WARNING (FDESCR): Because the call to CALSCT',
     $      /9x,'failed, an arbitrary liquid-like density will be',
     $      /9x,'assigned as a starting value for compressed liquid.')

            rhosl = 1.05d0
          ENDIF

        ELSE

          IF (qerrx) THEN

            WRITE (nttyo,1040)
            WRITE (noutpt,1040)

 1040       FORMAT(/6x,'WARNING (FDESCR): Because the call to CALSCT',
     $      /9x,'failed, the vapor and liquid states cannot be',
     $      /9x,'distinguished from one another. Liquid is assigned',
     $      /9x,'arbitrarily.')

            udescr = 'liquid'
            rhosl = 1.05d0
          ELSE

            IF (press .GE. psat) THEN
              udescr = 'liquid'
            ELSE
              udescr = 'vapor'
            ENDIF

c           Use density (rho) if available and pressure is close
c           to psat.

            IF (rho .GT. 0.0d0) THEN

              ptest = (press - psat)/psat
              btest = 10*btxtol

              IF (DABS(ptest) .LE. btest) THEN

c               Here press is very close to psat.
c               Use rho to determine vapor or liquid.

                rtestl = (rho - rhosl)/rhosl
                rtestv = (rho - rhosv)/rhosv

                IF (DABS(rtestl) .LE. rhotol) THEN
                  udescr = 'liquid'
                ELSEIF (DABS(rtestv) .LE. rhotol) THEN
                  udescr = 'vapor'
                ELSE

                  WRITE (nttyo,1050)
                  WRITE (noutpt,1050)

 1050             FORMAT(/6x,'WARNING (FDESCR): Could not use density',
     $            ' rho to discriminate',/6x,'vapor from liquid.')

                ENDIF

              ENDIF

            ENDIF

          ENDIF

        ENDIF

      ELSE

        IF (press .GT. pcr) THEN
          udescr = 'supercritical fluid'
        ELSE
          udescr = 'hot vapor'
        ENDIF

      ENDIF

  999 CONTINUE

      END
