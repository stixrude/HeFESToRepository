	double precision function Tspin(idum)

C  Optionally set upper bound to range of temperature
C  To remove the upper bound set Tspin = 0.
	
        include 'P1'
	integer idum
	double precision apar,Ti,Pi
        common /state/ apar(nspecp,nparp),Ti,Pi

	Tspin = 2300. + (4100.-2300.)*Pi/10.
	Tspin = 3400. + (5300.-3400.)*Pi/10.
C  Spinodal limits of inv010220 with the addition of artificial high temperature ol2 phase
	Tspin = 5100. + (7500.-5100.)*Pi/10.
c	Tspin = 10100. + (7500.-5100.)*Pi/10.

c	Tspin = 0.

	return
	end
