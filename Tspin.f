	double precision function Tspin(idum)

C  Optionally set upper bound to range of temperature
C  To remove the upper bound set Tspin = 0.
	
        include 'P1'
	integer idum
	double precision apar,Ti,Pi
        common /state/ apar(nspecp,nparp),Ti,Pi

	Tspin = 0.
C  Spinodal limits of inv010220 with the addition of artificial high temperature ol2 phase
c	Tspin = 5100. + (7500.-5100.)*Pi/10.
C  Spinodal limits of inv010123A with the addition of artificial high temperature ol2 phase
c	Tspin = 4900. + (7200.-4900.)*Pi/10.
C  Spinodal limits for MgSiO3 fluid
c	Tspin = 20000. + (10000000.-20000.)*Pi/10.

c	print*, 'End of Tspin',idum,Tspin

	return
	end
