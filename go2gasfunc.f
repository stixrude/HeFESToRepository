	double precision function go2gasfunc(Ti)

C  Feb. 22, 2025.  I no longer recall the source of the coeffficient values (a-d) below.  
C  The reference entropy is exactly that of the JANAF tables.
C  A spot check shows good agreement with the JANAF table values of H and S at other temperatures.
C  To do: use the NIST formulation, which is more involved, being piecewise continuous in temperature.

	double precision Ti
	double precision ao2,bo2,co2,do2,Tro2,STro2,ho2gas,so2gas

        ao2 = 47.255
        bo2 = -4.550e-4
        co2 = 4.402e+5
        do2 = -393.5
        Tro2 = 298.15
        STro2 = 205.147
        ho2gas = 0. + ao2*(Ti - Tro2) + 0.5*bo2*(Ti**2 - Tro2**2) - co2*(1./Ti - 1./Tro2) + 2.*do2*(sqrt(Ti) - sqrt(Tro2))
        so2gas = STro2 + ao2*log(Ti/Tro2) + bo2*(Ti - Tro2) - 0.5*co2*(1./Ti**2 - 1./Tro2**2)
     &          - 2.*do2*(1./sqrt(Ti) - 1./sqrt(Tro2))
        go2gasfunc = (ho2gas - Ti*so2gas + STro2*Tro2)/1000.

        write(31,*) 'oxygen gas enthalpy, entropy, gibbs (kJ/mol) = ',Ti,ho2gas/1000.,so2gas,go2gasfunc

	return
	end
