        double precision function vcon(x,d,idos)

        include 'P1'
        include 'const.inc'

	integer idos
	double precision x,d,sfac

        if (x .eq. 0.) then
         vcon = 0.
         return
        endif
        sfac = (2./pirad)**3

        if (idos .eq. 1) go to 10
        if (idos .eq. 2) go to 20
        if (idos .eq. 3) go to 30
        if (idos .eq. 4) go to 40

10      continue                                  !  Debye
        vcon = 3.*x*x
        return

20      continue                                  !  Einstein
        vcon = 1./d
        return

30      continue                                  !  Sin
        vcon = 3.*sfac*(asin(x))**2/sqrt(1. - x*x)
        return

40      continue                                  !  Optic Continuum
        vcon = 1./d
        return

        end
