        subroutine landaucr(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)

C  As in SLB11 except:
C  Coefficient of Q^6 is 1/3T_C (rather than 1/3T_C0)
C  Reference state is ORDERED.

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
c        if (Ti .ge. Tc) return

	qorder = 0.
	if (Ti .le. Tc) then
         qorder = ((Tc - Ti)/Tc)**(1./4.)
         cplan = smax*Ti/(2.*Tc*qorder**2)
         alplan = 1./2.*vmax/Vi*Ti/(Tc**2*qorder**2)
         betlan = 1./(2.*qorder**2)*vmax/Vi*vmax/(0.001*Tc*smax)*(Ti/Tc)**2
	end if

        glan = smax*((Ti - Tc)*qorder**2 + 1./3.*Tc*qorder**6 + 2./3.*Tc - Ti) 
        slan = -smax*(qorder**2 - 1.)
        vlan = -vmax*(qorder**2 - 1./3.*qorder**6 - 2./3.)

        return
        end
