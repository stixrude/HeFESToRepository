        subroutine landauqr(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)

C  As in SLB11 except:
C  Q^4 = (T_C - T)/T_C0
C  Reference state is ORDERED.
C  Note that in this formulation, qorder may exceed unity.  This may be problematic at very high pressure where e.g. quartz may become more stable than stishovite.
C  For now I have fixed this by limiting the value of qorder to qmax (1 July, 2016)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'
        include 'chem.inc'
	include 'lag.inc'

	integer ispec,iph,i
	double precision vi,qorder,tc,glan,cplan,alplan,betlan,slan,vlan,apar,Pi,smax,tco,Ti,vmax
	double precision t,gfunc,dgdt,d2gdt2,afe,bfe,cfe,ntot,x
        common /state/ apar(nspecp,nparp),Ti,Pi
	double precision, parameter :: qmax=2.0

	ntot = 0.
	iph = lphase(ispec)
	do 1 i=1,nspec
1	ntot = ntot + f(iph,i)*n(i)
	x = n(ispec)/ntot

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
C  Set magnetic entropy of bcc iron to zero at high pressure
C  Expression of jacobsschmidfetzer_10 Eq. 19-20
c	print*, 'Test for bcc',smax,Tco,abs(Tco-1043.01)
	if (abs(smax-9.46028) .le. 1e-5 .and. abs(Tco-1043.01) .le. 1e-5) then
c	if (abs(smax-9.46028) .le. 1e-6) then
	 if (Vi .le. 6.0) return
c	 print*, 'found bcc iron in landauqr'
	 t = Ti/Tco
	 alplan = 0.0
	 betlan = 0.0
	 afe = 0.9053
	 bfe = 0.91805
	 cfe = 0.64173
	 if (t .le. 1.0) then
	  gfunc = 1. - (afe/t + bfe*(t**3/6. + t**9/135. + t**15/600.))
	  dgdt = afe/t**2 - bfe*(t**2/2. + t**8/15. + t**14/40.)
	  d2gdt2 = -2.*afe/t**3 - bfe*(t + 8.*t**7/15. + 14.*t**13/40.)
	 else
	  gfunc = -cfe*(t**(-5)/10. + t**(-15)/315. + t**(-25)/1500.)
	  dgdt = cfe*(t**(-6)/2. + t**(-16)/21 + t**(-26)/60.)
	  d2gdt2 = -cfe*(3*t**(-7) + 16.*t**(-17)/21. + 26.*t**(-27)/60.)
	 end if
	 glan = smax*Ti*(gfunc - 1.)
	 slan = -smax*(gfunc - 1.) -smax*t*dgdt
	 cplan = -smax*(2.*t*dgdt + t**2*d2gdt2)
	 return
	end if
C  magnetite.  Assume same expression as for Fe.  For generalization, see Hillert's book, pg. 523.
c	print*, 'Test for mag',smax,Tco,abs(Tco-1043.01)
	if (abs(smax-43.1758) .le. 1.0 .and. abs(Tco-845.5) .le. 10.) then
c	if (abs(smax-43.1758) .le. 1.0 .and. abs(Tco-1845.5) .le. 10.) then
	 t = Ti/Tco
	 alplan = 0.0
	 betlan = 0.0
	 afe = 0.9053
	 bfe = 0.91805
	 cfe = 0.64173
	 if (t .le. 1.0) then
	  gfunc = 1. - (afe/t + bfe*(t**3/6. + t**9/135. + t**15/600.))
	  dgdt = afe/t**2 - bfe*(t**2/2. + t**8/15. + t**14/40.)
	  d2gdt2 = -2.*afe/t**3 - bfe*(t + 8.*t**7/15. + 14.*t**13/40.)
	 else
	  gfunc = -cfe*(t**(-5)/10. + t**(-15)/315. + t**(-25)/1500.)
	  dgdt = cfe*(t**(-6)/2. + t**(-16)/21 + t**(-26)/60.)
	  d2gdt2 = -cfe*(3*t**(-7) + 16.*t**(-17)/21. + 26.*t**(-27)/60.)
	 end if
	 glan = smax*Ti*(gfunc - 1.)
	 slan = -smax*(gfunc - 1.) -smax*t*dgdt
	 cplan = -smax*(2.*t*dgdt + t**2*d2gdt2)
	 return
	end if
        if (Tco .le. 0.) return
c        Tc = Tco + vmax/(0.001*smax)*Pi - Tco*(1. - x)
        Tc = Tco + vmax/(0.001*smax)*Pi
c        if (Ti .ge. Tc) return

	qorder = 0.
	if (Ti .le. Tc) then
         qorder = ((Tc - Ti)/Tco)**(1./4.)
	 if (qorder .gt. qmax) qorder = qmax
         cplan = smax*Ti/(2.*Tco*qorder**2)
         alplan = vmax/Vi/(2.*Tco*qorder**2)
         betlan = vmax/Vi*vmax/(0.001*smax)/(2.*Tco*qorder**2)
	end if
c	write(71,*) 'landauqr',ispec,Tco,vmax,smax,Pi,qorder,Tc,x

        glan = smax*((Ti - Tc)*(qorder**2 - 1.) + 1./3.*Tco*(qorder**6 - 1.))
        slan = -smax*(qorder**2 - 1.)
        vlan = -vmax*(qorder**2 - 1.)

        return
        end
