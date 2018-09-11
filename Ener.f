        double precision function Ener(x,d,idos)

        include 'P1'
        include 'dos.inc'
        include 'const.inc'

	integer idos,jlo,kl,klo,ku,mint
	double precision x,d,dy,sfac,xl,xu,y,yl,yu
C  Empirical low-T behavior of Sin
        double precision, parameter :: asin=23.594, bsin=6.1
C  Einstein limit of the width of the Optic Continuum
        double precision, parameter :: dmin = 0.01
        if (x .eq. 0.) then
         Ener = 0.
         return
        endif
        if (idos .eq. 5) then
         Ener = 1.0
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
        call neville(xa(klo),deb2(klo),mint,x,y,dy)
        Ener = y
        if (x .gt. xamax) Ener = (pirad*pirad*pirad*pirad/5.)/(x*x*x)
        if (x .lt. xamin) Ener = 1.0 - (1.0 - deb2(imin))/xamin*x

        return

20      continue                                  !  Einstein
        call neville(xa(klo),ein2(klo),mint,x,y,dy)
        Ener = y
        if (x .gt. xamax) Ener = x*exp(-x)
        if (x .lt. xamin) Ener = x/(exp(x) - 1.)
        return

30      continue                                  !  Sin
        call neville(xa(klo),sin2(klo),mint,x,y,dy)
        Ener = y
        if (x .gt. xamax) Ener = (pirad*pirad*pirad*pirad/5.)/(x*x*x)*sfac
     &                         + asin*(bsin-1.)*x**(1. - bsin)
        if (x .lt. xamin) Ener = 1.0 - (1.0 - sin2(imin))/xamin*x
        return

40      continue                                  !  Optic Continuum
        if (d .lt. dmin) then
         call bserch(xa,na,x,jlo)
         klo = min(max(jlo-(mint-1)/2,1),na+1-mint)
         go to 20
        end if
        xu = x*(1. + d/2.)
        xl = x*(1. - d/2.)
        call bserch(xa,na,xu,jlo)
        ku = min(max(jlo-(mint-1)/2,1),na+1-mint)
        call bserch(xa,na,xl,jlo)
        kl = min(max(jlo-(mint-1)/2,1),na+1-mint)
        call neville(xa(ku),opt2(ku),mint,xu,yu,dy)
        call neville(xa(kl),opt2(kl),mint,xl,yl,dy)
        if (xu .lt. xamin) yu = 1.0 - (1.0 - opt2(imin))/xamin*xu
        if (xl .lt. xamin) yl = 1.0 - (1.0 - opt2(imin))/xamin*xl
        if (xu .gt. xamax) yu = (pirad*pirad/6.)/(xu)
        if (xl .gt. xamax) yl = (pirad*pirad/6.)/(xl)
        Ener = (xu*yu - xl*yl)/(xu - xl)
        if (Ener .lt. 0.) Ener = 0.
        
        return
        end
