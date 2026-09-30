	double precision function gh2gasfunc(Ti)

C  Coefficient values from hollanpowell_11

	double precision Ti
	double precision ah2,bh2,ch2,dh2,Trh2,STrh2,hh2gas,sh2gas

        ah2 = 23.3
        bh2 = +4.627e-3
        ch2 = 0.
        dh2 = +76.3
        Trh2 = 298.15
        STrh2 = 130.680
        hh2gas = 0. + ah2*(Ti - Trh2) + 0.5*bh2*(Ti**2 - Trh2**2) - ch2*(1./Ti - 1./Trh2) + 2.*dh2*(sqrt(Ti) - sqrt(Trh2))
        sh2gas = STrh2 + ah2*log(Ti/Trh2) + bh2*(Ti - Trh2) - 0.5*ch2*(1./Ti**2 - 1./Trh2**2)
     &          - 2.*dh2*(1./sqrt(Ti) - 1./sqrt(Trh2))
        gh2gasfunc = (hh2gas - Ti*sh2gas + STrh2*Trh2)/1000.

        write(31,*) 'hydrogen gas enthalpy, entropy, gibbs (kJ/mol) = ',Ti,hh2gas/1000.,sh2gas,gh2gasfunc

	return
	end
