        subroutine thermlel(ispec,Vi,Fel,Eel,Sel,Pel,Cvel,betael,Kel,Kelp,aktel,daktdvel)
	
	include 'P1'
        include 'theory.inc'
	include 'const.inc'

	integer ispec
	double precision eta,aktel,apar,d2tdv2,d2xdv2,d2zdv2,d3tdv3,d3xdv3,d3zdv3
	double precision daktdvel,dtdv,dxdv,dzdv,fel
        double precision Kel,Kelp,pel,Pi,sel,telo,Ti,x
	double precision Eel,dlnzdv,d2lnzdv2,dlntdv,d2lntdv2,d3lnzdv3,d3lntdv3,Tel,DTel,Tinf
	double precision zelo,Vo,xi,y,zel,bpar,Vi
	double precision betael,Cvel,zelr
        common /state/ apar(nspecp,nparp),Ti,Pi

        Vo =   apar(ispec,6)
	Telo = apar(ispec,13)
	eta = apar(ispec,14)
	Tinf = apar(ispec,15)
	zelo = apar(ispec,28)
	xi = apar(ispec,29)
	bpar = apar(ispec,30)
	zelo = zelo/1.e6
	Tel = Tinf + Telo*(Vi/Vo)**eta

        x = Telo/Ti
	DTel = Tel - Tinf
        dtdv = eta*DTel/Vi
        d2tdv2 = DTel/Vi**2*eta*(eta-1.)
	d3tdv3 = DTel/Vi**3*eta*(eta-1.)*(eta-2.)
	dxdv = dtdv/Ti
	d2xdv2 = d2tdv2/Ti
	d3xdv3 = d3tdv3/Ti
        dlntdv = dtdv/Tel
        d2lntdv2 = -dlntdv**2 + d2tdv2/Tel
        d3lntdv3 = 2.*dlntdv**3 - 3.*dlntdv*d2tdv2/Tel + d3tdv3/Tel
        y = Vi/Vo
        zel = zelo/(bpar*y**(-xi) + 1.)
        zelr = 1./(bpar*y**(-xi) + 1.)
        dzdv = xi*zel/Vi*(1. - zelr)
        d2zdv2 = dzdv/Vi*(xi - 2.*xi*zelr - 1.)
        d3zdv3 = d2zdv2/Vi*(xi - 2.*xi*zelr - 2.) - dzdv**2/Vi*2.*xi*zelr/zel
        dlnzdv = dzdv/zel
        d2lnzdv2 = -dlnzdv**2 + d2zdv2/zel
        d3lnzdv3 = 2.*dlnzdv**3 - 3.*dlnzdv*d2zdv2/zel + d3zdv3/zel

        Sel = zel*(Ti - Tel - Tel*log(Ti/Tel))
        Sel = 0.5*zel*(log(Ti/Tel))**2					! kJ/mol/K
        Cvel = zel*max((Ti - Tel),0.)
        Cvel = zel*log(Ti/Tel)						! kJ/mol/K
        Eel = 0.5*zel*(Ti - Tel)**2
        Eel = zel*(Ti*log(Ti/Tel) - (Ti - Tel))				! kJ/mol
        betael = 1000.*dzdv*(Ti - Tel) - 1000.*zel*dtdv
        betael = 1000.*(dlnzdv*Cvel - zel*dlntdv)			! J/mol/K/(cm^3/mol)
        if (Ti .le. Tel .or. Tel .eq. 0.) then
         Sel = 0.
         Cvel = 0.
         Eel = 0.
         betael = 0.
        end if

	Fel = Eel - Ti*Sel						! kJ/mol

        Pel = -dlnzdv*Fel + zel*dtdv*((Ti - Tel) - Ti*log(Ti/Tel))
	Pel = -dlnzdv*Fel - dlntdv*Eel					! GPa
	if (Ti .le. Tel .or. Tel .eq. 0.) Pel = 0.

        aktel = dlnzdv*Sel - zel*dtdv*log(Ti/Tel)
	aktel = dlnzdv*Sel - dlntdv*Cvel				! GPa/K
	if (Ti .le. Tel .or. Tel .eq. 0.) aktel = 0.

	daktdvel = d2lnzdv2*Sel + dlnzdv*aktel - dzdv*dtdv*log(Ti/Tel) - zel*d2tdv2*log(Ti/Tel) + zel*dtdv/Tel*dtdv
	daktdvel = d2lnzdv2*Sel + dlnzdv*aktel - d2lntdv2*Cvel - 0.001*dlntdv*betael	! GPa/K/(cm^3/mol)
	if (Ti .le. Tel .or. Tel .eq. 0.) daktdvel = 0.

	Kel = -Vi*(-d2lnzdv2*Fel + dlnzdv*Pel + (dzdv*dtdv + zel*d2tdv2)*((Ti - Tel) - Ti*log(Ti/Tel))
     &      + zel*dtdv*(-dtdv + Ti/Tel*dtdv))
	Kel = -Vi*((dlnzdv + dlntdv)*Pel - d2lnzdv2*Fel - d2lntdv2*Eel - dlntdv*aktel*Ti)
	if (Ti .le. Tel .or. Tel .eq. 0.) Kel = 0.

 	Kelp = -1. +Vi**2/Kel*Ti**2*(
     &         d3zdv3*(0.5*(1. - x**2) + x*log(x)) 
     &  +      (3.*d2zdv2*dxdv + 3.*dzdv*d2xdv2 + zel*d3xdv3)*(1. - x + log(x))
     &  +      (3.*dzdv*dxdv**2 + 3.*zel*dxdv*d2xdv2)*(1./x - 1.)
     &  -      zel*dxdv**3/x**2)
	Kelp = -1. + Vi**2/Kel*(2.*(d2lnzdv2 + d2lntdv2)*Pel - (dzdv + dtdv)*Kel/Vi
     &       - d3lnzdv3*Fel - d3lntdv3*Eel - 2*d2lntdv2*aktel*Ti - dlntdv*daktdvel*Ti)
	if (Ti .le. Tel .or. Tel .eq. 0. .or. Kel .eq. 0.) Kelp = 0.

	return
	end
