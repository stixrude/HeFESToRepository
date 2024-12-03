      SUBROUTINE LEJUST(ustr)

c     This subroutine left-justifies the non-blank portion of the string
c     ustr.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       ustr   = the output string variable

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j,jj,jbl,j1,nchars

      INTEGER IFNOBL

c-----------------------------------------------------------------------

c     Get the length of the string variable.

      nchars = len(ustr)

c     Get the position of the first non-blank character and the number
c     of blanks on the left-hand-side.

      j1 = IFNOBL(ustr)
      jbl = j1 - 1

      IF (jbl .GT. 0) THEN
        DO jj = j1,nchars
          j = jj - jbl
          ustr(j:j) = ustr(jj:jj)
        ENDDO
        DO j = nchars - jbl + 1,nchars
          ustr(j:j) = ' '
        ENDDO
      ENDIF

      END
