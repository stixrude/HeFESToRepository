      INTEGER FUNCTION ILNOBL(ustr)

c     This subroutine finds the position of the last non-blank character
c     in the string ustr.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       ILNOBL = the position of the first non-blank character

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

      ILNOBL = 0
      DO j = nchars,1,-1
        IF (ustr(j:j) .NE. ' ') THEN
          ILNOBL = j
          GO TO 999
        ENDIF
      ENDDO

  999 CONTINUE

      END
