	double precision function videal(fn,Pi,Ti)

	double precision fn,Pi,Ti

	include 'P1'
	include 'const.inc'

	videal = fn*Rgas*Ti/Pi/1000.

	return
	end
