      INTEGER FUNCTION IFNOBL(ustr)

c     This subroutine finds the position of the first non-blank
c     character in the string ustr.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       IFNOBL = the position of the first non-blank character

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j,nchars

c-----------------------------------------------------------------------

c     Get the length of the string variable.

      nchars = len(ustr)

c     Find the first non-blank character.

      IFNOBL = 0
      DO j = 1,nchars
        IF (ustr(j:j) .NE. ' ') THEN
          IFNOBL = j
          GO TO 999
        ENDIF
      ENDDO

  999 CONTINUE

      END
