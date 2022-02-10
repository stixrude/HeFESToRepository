        double precision function Heat(x,d,idos)

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
         Heat = 0.
         return
        end if
        if (idos .eq. 5) then
         Heat = 1.0
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
        call neville(xa(klo),deb1(klo),mint,x,y,dy)
        Heat = y
        if (x .gt. xamax) Heat = (12.*pirad*pirad*pirad*pirad/15.)/(x*x*x)
        if (x .lt. xamin) Heat = 1.0 - (1.0 - deb1(imin))/xamin**2*x*x
        return

20      continue                                  !  Einstein
        call neville(xa(klo),ein1(klo),mint,x,y,dy)
        Heat = y
        if (x .gt. xamax) Heat = x*x*exp(-x)
        if (x .lt. xamin) Heat = 1.0 - (1.0 - ein1(imin))/xamin**2*x*x
        return

30      continue                                  !  Sin
        if (x .gt. xamax) then
	 Heat = (12.*pirad*pirad*pirad*pirad/15.)/(x*x*x)*sfac 
     &                         + asin*bsin*(bsin-1)*x**(1-bsin)
	 return
	end if
        if (x .lt. xamin) then
	 Heat = 1.0 - (1.0 - sin1(imin))/xamin**2*x*x
	 return
	end if
        call neville(xa(klo),sin1(klo),mint,x,y,dy)
        Heat = y
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
        call neville(xa(ku),opt1(ku),mint,xu,yu,dy)
        call neville(xa(kl),opt1(kl),mint,xl,yl,dy)
        if (xu .lt. xamin) yu = 1.0 - (1.0 - opt1(imin))/xamin**2*xu*xu
        if (xl .lt. xamin) yl = 1.0 - (1.0 - opt1(imin))/xamin**2*xl*xl
        if (xu .gt. xamax) yu = (pirad*pirad/3.)/(xu)
        if (xl .gt. xamax) yl = (pirad*pirad/3.)/(xl)
        Heat = (xu*yu - xl*yl)/(xu - xl)
        if (Heat .lt. 0.) Heat = 0.

        return
        end
