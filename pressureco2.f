	double precision function pressureco2(Vi)
       USE def_constants
       USE properties
	include 'P1'
	double precision apar,Ti,Pi
	double precision Vi,v
        common /state/ apar(nspecp,nparp),Ti,Pi

	v = Vi/W/1.d6				! m^3/kg
	call pressure_sw(Ti,v,pressureco2)	! Pa

	pressureco2 = pressureco2/1.d9 - Pi

c	print '(a11,99d16.5)', 'pressureco2',v,Vi,Ti,Pi,pressureco2

	return
	end
