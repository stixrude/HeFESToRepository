	double precision function vmurnaghan(ispec)

	include 'P1'
	include 'const.inc'
	integer ispec
	double precision apar,Ti,Pi,Vo,Ko,Kop,gamma,fn,Pth

        common /state/ apar(nspecp,nparp),Ti,Pi

	fn = apar(ispec,1)
        Vo = apar(ispec,6)
        Ko = apar(ispec,7)
        Kop = apar(ispec,8)
	gamma = apar(ispec,26)
	Pth = gamma*3.*fn*Rgas*Ti/Vo/1000.

        vmurnaghan = Vo*(1. + (Pi - Pth)*(Kop/Ko))**(-1./Kop)

	return
	end
