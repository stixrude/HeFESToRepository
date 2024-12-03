      SUBROUTINE GAUSSE(noutpt, nttyo, qerr, qprnt3, kdim,
     $ kmax, aamatr, alpha, deltas)

c     Solve the matrix equation aamatr*deltas = alpha using Gaussian
c     elimination with partial pivoting.

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qerr, qprnt3

      INTEGER nttyo, noutpt

      INTEGER kdim, kmax

      REAL(8) aamatr(kmax,kmax), alpha(kmax), deltas(kmax)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      INTEGER i, j, k, kk, kopt

      REAL(8) aamx, abmx, amx, axx

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      qerr = .FALSE.

      DO k = 1,kdim - 1

c       Find the optimal pivot row in rows k through kdim
c       Unless row k happens to be the optimal row, exchange
c       it for the optimal row.

        kopt = k
        aamx = DABS(aamatr(k,k))

        DO kk = k + 1, kdim
          abmx = DABS(aamatr(kk,k))
          IF (abmx .GT. aamx) THEN
            kopt = kk
            aamx = abmx
          ENDIF
        ENDDO

        IF (aamx .LE. 1.0d-14) THEN

          WRITE(nttyo,1120)
          WRITE(noutpt,1120)

 1120     FORMAT(/6x,'ERROR (GAUSSE): Have a singular matrix',/)

          qerr = .FALSE.
          GO TO 999
        ENDIF

        IF (kopt .GT. k) THEN

c         Exchange rows k and kopt in the extended matrix
c         (don't forget the right-hand-side part).

          DO j =  k,kdim
            amx = aamatr(k,j)
            aamatr(k,j) = aamatr(kopt,j)
            aamatr(kopt,j) = amx
          ENDDO
          axx = alpha(k)
          alpha(k) = alpha(kopt)
          alpha(kopt) = axx

          IF (qprnt3) THEN

            WRITE (noutpt,1130)

 1130       FORMAT(/3x,'Pivot (row exchange)')

            WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $      alpha(i), i = 1,kdim)

 1140       FORMAT(/6x,'Jacobian and residual',
     $      /9x,1pe12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $      /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $      /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5)

          ENDIF

        ENDIF

        amx = aamatr(k,k)
        aamatr(k,k) = 1.0d0
        DO j =  k + 1,kdim
          aamatr(k,j) = aamatr(k,j)/amx
        ENDDO
        alpha(k) = alpha(k)/amx

        IF (qprnt3) THEN

          WRITE (noutpt,1150)

 1150     FORMAT(/3x,'Scale the pivot row for unit diagonal',
     $    ' entry')

          WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $    alpha(i), i = 1,kdim)

        ENDIF

        DO kk = k + 1,kdim

          axx = aamatr(kk,k)

          DO j = k,kdim
            aamatr(kk,j) = aamatr(kk,j) - axx*aamatr(k,j)
          ENDDO
          alpha(kk) = alpha(kk) - axx*alpha(k)

        ENDDO

        IF (qprnt3) THEN

          WRITE (noutpt,1160)

 1160     FORMAT(/3x,'Zero lower row entries in the pivot',
     $    ' column')

          WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $    alpha(i), i = 1,kdim)

        ENDIF

      ENDDO

c     Scale the last row for unit diagonal entry.

      amx = aamatr(k,k)
      aamatr(k,k) = 1.0d0
      DO j =  k + 1,kdim
        aamatr(k,j) = aamatr(k,j)/amx
      ENDDO
      alpha(k) = alpha(k)/amx

      IF (qprnt3) THEN
        WRITE (noutpt,1150)
        WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $  alpha(i), i = 1,kdim)
      ENDIF

c     Find the solution vector.

      DO k = kdim,1,-1
        deltas(k) = alpha(k)
        DO kk = k + 1,kdim
          deltas(k) = deltas(k) - aamatr(k,kk)*deltas(kk)
        ENDDO
      ENDDO

c     deltas(3) = alpha(3)
c     deltas(2) = ( alpha(2) - aamatr(2,3)*deltas(3) )
c     deltas(1) = ( alpha(1) - aamatr(1,2)*deltas(2)
c    $ - aamatr(1,3)*deltas(3) )

  999 CONTINUE

      END
