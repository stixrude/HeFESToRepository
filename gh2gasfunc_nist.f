	double precision function gh2gasfunc(Ti)

C  NIST formulation, which is piecewise continuous in temperature.
C  https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1#Thermo-Gas

	integer iwindow
	double precision Ti,t
	double precision ah2,bh2,ch2,dh2,eh2,fh2,gh2,hh2,Trh2,STrh2,hh2gas,sh2gas
	double precision aa(3),ba(3),ca(3),da(3),ea(3),fa(3),ga(3),ha(3)
	data aa/33.066178,	18.563083,	43.413560/
	data ba/-11.363417,	12.257357,	-4.293079/
	data ca/11.432816,	-2.859786,	1.272428/
	data da/-2.772874,	0.268238,	-0.096876/
	data ea/-0.158558,	1.977990,	-20.533862/
	data fa/-9.980797,	-1.147438,	-38.515158/
	data ga/172.707974,	156.288133,	162.081354/
	data ha/0.0,	0.0,	0.0/

	t = Ti/1000.
	if (Ti .le. 1000.) iwindow = 1
	if (Ti .gt. 1000. .and. Ti. le. 2500.) iwindow = 2
	if (Ti .gt. 2500.) iwindow = 3
	ah2 = aa(iwindow)
	bh2 = ba(iwindow)
	ch2 = ca(iwindow)
	dh2 = da(iwindow)
	eh2 = ea(iwindow)
	fh2 = fa(iwindow)
	gh2 = ga(iwindow)
	hh2 = ha(iwindow)
	write(31,*) 'hydrogen coefficients',ah2,bh2,ch2,dh2,eh2,fh2,gh2,hh2

        Trh2 = 298.15
        STrh2 = 130.680
	hh2gas = ah2*t + bh2*t**2/2. + ch2*t**3/3. + dh2*t**4/4. - eh2/t + fh2 - hh2
	sh2gas = ah2*log(t) + bh2*t + ch2*t**2/2. + dh2*t**3/3. - eh2/(2.*t**2) + gh2
        gh2gasfunc = (1000.*hh2gas - Ti*sh2gas + STrh2*Trh2)/1000.

        write(31,*) 'hydrogen gas enthalpy, entropy, gibbs (kJ/mol) = ',Ti,hh2gas,sh2gas,gh2gasfunc

	return
	end
