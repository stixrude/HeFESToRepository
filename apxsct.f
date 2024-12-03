      SUBROUTINE APXSCT(noutpt, nttyo, qprnt1, qsilent, rhocr, rhosv,
     $ rhosl, tempk, tcr, pcr, psat, psatt, qerr)

c     Evaluate the pressure (psat) as a function of temperature along
c     the vapor-liquid equilibrium curve, using equation 2.5 of Wagner
c     and Pruss (2002). Also calculate the densities of the liquid and
c     vapor phases using equations 2.6 and 2.7 from the same source.
c     Results are not fully consistent with IAPWS-95 to high precision,
c     but may serve as close approximations or starting values for
c     refinement.

c       psat = saturation pressure (MPa)
c       rhosv = density of vapor (kg/m3)
c       rhosl = density of liquid (kg/m3)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qprnt1, qsilent, qerr

      INTEGER nttyo, noutpt

      REAL(8) delta, rhosv, rhosl, rhocr, tempk, tcr, pcr

      REAL(8) psat, psatt

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

c     Model coefficients.

      REAL(8) c1a(1:6), c1b(1:6), c1c(1:6)
      REAL(8) e1b(1:6), e1c(1:6)

      DATA c1a(1:6) / -7.85951783d0, 1.84408259d0, -11.7866497d0,
     $ 22.6807411d0, -15.9618719d0, 1.80122502d0 /

      DATA c1b(1:6) / 1.99274064d0, 1.09965342d0, -0.510839303d0,
     $ -1.75493479d0, -45.5170352d0, -6.74694450d+05 /

      DATA c1c(1:6) / -2.03150240d0, -2.68302940d0, -5.38626492d0,
     $ -17.2991605d0, -44.7586581d0, -63.9201063d0 /

      LOGICAL qfirst

      SAVE qfirst, e1b, e1c

      DATA qfirst / .TRUE. /

      REAL(8) varth, x1, x2

      INTEGER i

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (qprnt1 .AND. .NOT.qsilent) THEN

        WRITE(nttyo,1100) tempk
        WRITE(noutpt,1100) tempk

 1100   FORMAT(/3x,'APXSCT: Input tempk= ',f8.4' (K)',/)

      ENDIF

c     If this is the first call, calculate the exponents for the
c     density equations.

      IF (qfirst) THEN

        e1b(1) = 1.0d0/3.0d0
        e1b(2) = 2.0d0/3.0d0
        e1b(3) = 5.0d0/3.0d0
        e1b(4) = 16.0d0/3.0d0
        e1b(5) = 43.0d0/3.0d0
        e1b(6) = 110.0d0/3.0d0

        e1c(1) = 2.0d0/6.0d0
        e1c(2) = 4.0d0/6.0d0
        e1c(3) = 8.0d0/6.0d0
        e1c(4) = 18.0d0/6.0d0
        e1c(5) = 37.0d0/6.0d0
        e1c(6) = 71.0d0/6.0d0

        qfirst = .FALSE.

      ENDIF

c     Zero parameters to be calculated.

      psat = 0.0d0
      rhosv = 0.0d0
      rhosl = 0.0d0

c     Check to see that the temperature is in the allowed range.

      qerr = .FALSE.
      IF (tempk.LT.273.15d0 .OR. tempk.GT.tcr) THEN
        qerr = .TRUE.
        GO TO 999
      ENDIF

      varth = 1.0d0 - (tempk/tcr)

c     Saturation pressure.

      x1 = c1a(1)*varth + c1a(2)*varth**1.5d0 + c1a(3)*varth**3.0d0
     $ + c1a(4)*varth**3.5d0 + c1a(5)*varth**4.0d0
     $ + c1a(6)*varth**7.5d0

      x2 = ( tcr/tempk )*x1

      psat = pcr*DEXP(x2)

c     Derivative of saturation pressure with respect to temperature.

      x1 = c1a(1) + 1.5d0*c1a(2)*varth**0.5d0
     $ + 3.0d0*c1a(3)*varth**2.0d0 + 3.5d0*c1a(4)*varth**2.5d0
     $ + 4.0d0*c1a(5)*varth**3.0d0 + 7.5d0*c1a(6)*varth**6.5d0

      x2 = DLOG( psat/pcr ) + x1

      psatt = -( psat/tempk )*x2

c     Density of liquid.

      x1 = 1.0d0 + c1b(1)*varth**(1.0d0/3.0d0)
     $ + c1b(2)*varth**(2.0d0/3.0d0) + c1b(3)*varth**(5.0d0/3.0d0)
     $ + c1b(4)*varth**(16.0d0/3.0d0) + c1b(5)*varth**(43.0d0/3.0d0)
     $ + c1b(6)*varth**(110.0d0/3.0d0)

      rhosl = rhocr*x1

c     Density of vapor.

      x1 = c1c(1)*varth**(2.0d0/6.0d0) + c1c(2)*varth**(4.0d0/6.0d0)
     $ + c1c(3)*varth**(8.0d0/6.0d0) + c1c(4)*varth**(18.0d0/6.0d0)
     $ + c1c(5)*varth**(37.0d0/6.0d0) + c1c(6)*varth**(71.0d0/6.0d0)

      rhosv = rhocr*DEXP(x1)

  999 CONTINUE

      END
