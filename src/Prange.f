	subroutine Prange(Psmall,Plarge)

C  Must set a range of pressure over which to search
C  Psmall can be a negative quantity.  Psmall should be greater than the pressure of the spinodal

        include 'P1'
        include 'chem.inc'
        include 'absent.inc'
	include 'const.inc'
	logical isochor
	integer ispec,ires
	double precision apar,Ti,Pi,Psmall,Plarge,Vo,Ko,Kop,fspin,Pspin,ph,v1,v2,P1,P2,f1,f2
	double precision vlan,x1,x2,xx,fret,htl
	double precision, parameter :: tol=1.e-6
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ ispec,isochor

	Psmall = -1.e15
	Plarge = +1.e15
	do 1 ispec=1,nspec
	 if (absents(ispec)) go to 1
         vlan = apar(ispec,40)
	 Vo = apar(ispec,6)
	 Ko = apar(ispec,7) 
	 Kop = apar(ispec,8) 
	 htl = apar(ispec,31)
	 ph = 0.001*apar(ispec,26)/Vo*3.*apar(ispec,1)*Rgas*(Ti - apar(ispec,4))
C  Set the lower limit of the pressure range to the spinodal.
C  Estimate the strain, fspin, and the pressure, Pspin at spinodal from second order finite strain theory about the current point
	 fspin = -(5./3.*Pi + Ko)/(7.*Ko)
	 Pspin = 3.*Ko*f1*(1. + 2.*f1)**(5./2.) + ph
C  Or estimate the strain from the volume of the spinodal at T=T_0
	 v1 = max(apar(ispec,51),apar(ispec,53))
	 v2 = min(apar(ispec,52),apar(ispec,54))
	 f1 = 0.5*((v1/Vo)**(-2./3.) - 1.)
	 f2 = 0.5*((v2/Vo)**(-2./3.) - 1.)
	 P1 = 3.*Ko*f1*(1. + 2.*f1)**(5./2.)*(1. + 1.5*(Kop - 4.)*f1)
	 P2 = 3.*Ko*f2*(1. + 2.*f2)**(5./2.)*(1. + 1.5*(Kop - 4.)*f2) + ph
	 Plarge = min(Plarge,P1)
	 Psmall = max(Psmall,P2)
C  Or find the volume of the spinodal at T=Ti by minimizing the pressure.  Assume that Vo<V_spinodal<Vsp(T_0)
         x1 = Vo - vlan
         x2 = v2
         xx = x2*(1. - tol)
	 if (htl .eq. 0.) then
          call nlmin_V(xx,x1,x2,fret,ires)
	 else
	  call nlmin_VL(xx,x1,x2,fret,ires)
	 end if
	 P2 = fret + Pi
	 Psmall = max(Psmall,P2)
c	 write(31,*) 'Prange prelim',ispec,v1,v2,f1,f2,P1,P2,Plarge,Psmall,xx,fret,x1,x2,ires,Pi,htl
1	continue

	write(31,*) 'Prange',Pi,Ti,Psmall,Plarge

	return
	end
