        subroutine landau(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)

C  As in SLB11

        include 'P1'
        include 'const.inc'
        include 'theory.inc'

	integer ispec
	double precision vi,qorder,tc,glan,cplan,alplan,betlan,slan,vlan,apar,Pi,smax,tco,Ti,vmax
        common /state/ apar(nspecp,nparp),Ti,Pi

        glan = 0.
        slan = 0.
        vlan = 0.
        cplan = 0.
        alplan = 0.
        betlan = 0.
        qorder = 0.
        Tc = 0.
        Tco = apar(ispec,38)
        smax = apar(ispec,39)
        vmax = apar(ispec,40)
        if (Tco .le. 0.) return
        Tc = Tco + vmax/(0.001*smax)*Pi
        if (Ti .ge. Tc) return

        qorder = ((Tc - Ti)/Tc)**(1./4.)

        glan = smax*((Ti - Tc)*qorder**2 + 1./3.*Tco*qorder**6)
        slan = -smax*qorder**2*(1.5 - 0.5*Tco/Tc)
        vlan = -vmax*qorder**2*(1. + 0.5*Ti/Tc*(1. - Tco/Tc))
        cplan = smax*Ti/(2.*Tc*qorder**2)*(1.5 - 0.5*Tco/Tc)
        alplan = 1./2.*vmax/Vi*Ti/(Tc**2*qorder**2)
        alplan = 1./2.*vmax/Vi/(Tc*qorder**2)
        alplan = 0.5*vmax/Vi/Tc*(1./qorder**2*(1. + 0.5*Ti/Tc*(1. - Tco/Tc))
     &         - qorder**2*(1. - Tco/Tc))
        betlan = 1./(2.*qorder**2)*vmax/Vi*vmax/(0.001*Tc*smax)*Ti/Tc
        betlan = vmax/Vi*vmax/(0.001*Tc*smax)*Ti/Tc*
     &         (  0.5/qorder**2*(1. + 0.5*Ti/Tc*(1. - Tco/Tc))
     &         +  qorder**2*(Tco/Tc - 0.5))

        return
        end
