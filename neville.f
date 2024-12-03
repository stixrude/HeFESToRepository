       SUBROUTINE Neville(A, B, xm, xt, res, err)
C  Standard method for polynomial interpolation 
C  Given arrays of x,y values A,B find the y-value (res) at the given x-value (xt) and estimate the uncertainty (err)
C  xm is the order of the interpolation
       INTEGER, intent ( in ) :: xm
       double precision, dimension(xm), intent( in ) :: A, B
       double precision, intent( in ) :: xt
       double precision, intent( out ) :: res, err
       double precision, dimension(xm) :: P
       INTEGER j,k
       P(:) = B(:)
       DO j=1, xm
          DO k=1, xm - j
               P(k)=((xt - A(k+j)) * P(k) + (A(k)-xt) * P(k+1))/(A(k)-A(k+j))
            END DO
       END DO
       res=P(1)
       err=P(1)-P(2)
       return
       END SUBROUTINE Neville
