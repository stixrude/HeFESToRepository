        double precision function Helm(x,d,idos)

        include 'P1'
        include 'dos.inc'
        include 'const.inc'

	integer idos,jlo,kl,klo,ku,mint
	double precision x,d,c2,c4,dx,dy,fein,sfac,u,v,xl,xu,y,yl,yu
C  Empirical low-T behavior of Sin
        double precision, parameter :: asin=23.594, bsin=6.1
C  non-Empirical high-T behavior of F of Sin
C  parameter comes from difference between theta(0) and theta_i for the 
C  sin dispersion relation (see notes) 8/30/01
        double precision, parameter :: csin=0.188256
C  Einstein + Taylor Series in d*x limit of the width of the Optic Continuum
        double precision, parameter :: dmin = 0.4
        if (x .eq. 0.) then
         Helm = 0.
         return
        end if
        if (idos .eq. 5) then
         Helm = log(x) - 1./3.
         return
        end if
        mint = 4
        sfac = (2./pirad)**3
        if (idos .ne. 4) then
         call bserch(xa,na,x,jlo)
         klo = min(max(jlo-(mint-1)/2,1),na+1-mint)
        end if

        if (idos .eq. 1) go to 10
        if (idos .eq. 2) go to 20
        if (idos .eq. 3) go to 30
        if (idos .eq. 4) go to 40

10      continue                                  !  Debye
        call neville(xa(klo),deb3(klo),mint,x,y,dy)
        Helm = y + log(1. - exp(-x))
        if (x .gt. xamax) Helm = -(pirad*pirad*pirad*pirad/15.)/(x*x*x)
        if (x .lt. xamin) Helm = log(x) - 4./3.
     &                         + 1.0 - (1.0 - deb2(imin))/xamin*x
        return

20      continue                                  !  Einstein
        Helm = log(1. - exp(-x))
        if (x .gt. xamax) Helm = - exp(-x)
        return

30      continue                                  !  Sin
        call neville(xa(klo),sin3(klo),mint,x,y,dy)
        Helm = y + log(1. - exp(-x))
        if (x .gt. xamax) Helm = -(pirad*pirad*pirad*pirad/15.)/(x*x*x)*sfac
     &                         - asin*x**(1. - bsin)
        if (x .lt. xamin) Helm = log(x) - 4./3.
     &                         + 1.0 - (1.0 - sin2(imin))/xamin*x
     &                         + csin
        return

40      continue                                  !  Optic Continuum
        if (d .le. dmin) then
         dx = x*d/2.
         u = exp(-x)
         v = u/(1. - u)
         fein = log(1. - u)
         c2 = -2.*v - 2.*v*v
         c4 = -2.*v - 14.*v*v - 24.*v*v*v - 12.*v*v*v*v
         Helm = fein + c2*dx*dx/(2.*6.) + c4*dx*dx*dx*dx/(2.*120.)
         return
        end if
        xu = x*(1. + d/2.)
        xl = x*(1. - d/2.)
        call bserch(xa,na,xu,jlo)
        ku = min(max(jlo-(mint-1)/2,1),na+1-mint)
        call bserch(xa,na,xl,jlo)
        kl = min(max(jlo-(mint-1)/2,1),na+1-mint)
        call neville(xa(ku),opt3(ku),mint,xu,yu,dy)
        call neville(xa(kl),opt3(kl),mint,xl,yl,dy)
        if (xu .ne. 0.) yu = yu + log(1. - exp(-xu))
        if (xl .ne. 0.) yl = yl + log(1. - exp(-xl))
        if (xu .lt. xamin) yu = log(1. - exp(-xu)) 
     &   - (1.0 - (1.0 - opt2(imin))/xamin*xu)
        if (xl .lt. xamin) yl = log(1. - exp(-xl)) 
     &   - (1.0 - (1.0 - opt2(imin))/xamin*xl)
        if (xu .gt. xamax) yu = -(pirad*pirad/6.)/(xu)
        if (xl .gt. xamax) yl = -(pirad*pirad/6.)/(xl)
        if (xl .eq. 0.) yl = 0.
        if (xu .eq. 0.) yu = 0.
        Helm = (xu*yu - xl*yl)/(xu - xl) 
c       print*, xu,yu,xl,yl
        
        return
        end
