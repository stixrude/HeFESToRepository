        subroutine landauqr(ispec,Vi,qorder,Tc,glan,cplan,alplan,betlan,slan,vlan)

C  As in SLB11 except:
C  Q^4 = (T_C - T)/T_C0
C  Reference state is ORDERED.
C  Note that in this formulation, qorder may exceed unity.  This may be problematic at very high pressure where e.g. quartz may become more stable than stishovite.
C  For now I have fixed this by limiting the value of qorder to qmax (1 July, 2016)
C  Reduce qmax from 2 to 1.5 (2 April, 2022).  The reason is that with qmax=2, neph is stabilized in pyrolite at P>400 GPa along a 1600 K adiabat (conditions relevant for super-Earths).
C  Limit Tc instead of q.  The reason is that Tc increases without bound with increasing P, so glan, which goes like -Tc, decreases without bound.  Limit Tc/Tco<=10. (July, 2024)

        include 'P1'
        include 'const.inc'
        include 'theory.inc'
        include 'chem.inc'
	include 'lag.inc'

	integer ispec,iph,i,icfe
	double precision vi,qorder,tc,glan,cplan,alplan,betlan,slan,vlan,apar,Pi,smax,tco,Ti,vmax
	double precision t,gfunc,dgdt,d2gdt2,afe,bfe,cfe,ntot,x,p,glansave
        common /state/ apar(nspecp,nparp),Ti,Pi
	common /mag/ icfe
	double precision, parameter :: qmax=1.5
	double precision, parameter :: Tcmax=10.

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
        if (Tco .le. 0.) return
c	print*, 'Test for bcc',smax,Tco,abs(Tco-1043.01)
	if (abs(smax-9.46028) .le. 1e-5 .and. abs(Tco-1043.01) .le. 1e-5) then
c	if (icfe .ne. 0 .and. s(icfe,ispec) .ne. 0.) then
	 alplan = 0.0
	 betlan = 0.0
	 p = 0.40
	 call hillert(Ti,Tco,p,smax,glan,slan,cplan)
	 return
	end if
C  magnetite.  Assume same expression as for Fe.  For generalization, see Hillert's book, pg. 523.
c	if (abs(smax-43.1758) .le. 1.0 .and. abs(Tco-845.5) .le. 10.) then
	if (abs(Tco-845.5) .le. 1.) then
	 alplan = 0.0
	 betlan = 0.0
	 p = 0.40
	 call hillert(Ti,Tco,p,smax,glan,slan,cplan)
	 return
	end if
C  wustite.  Assume same expression as for Fe.  For generalization, see Hillert's book, pg. 523.
c	if (abs(smax-53.52540) .le. 1.0 .and. abs(Tco-191.0) .le. 10.) then
c	 alplan = 0.0
c	 betlan = 0.0
c	 p = 0.60
c	 call hillert(Ti,Tco,p,smax,glan,slan,cplan)
c	 return
c	end if

        Tc = Tco + vmax/(0.001*smax)*Pi
	if (Tc .gt. Tcmax*Tco) Tc = Tcmax*Tco

	qorder = 0.
	if (Ti .le. Tc) then
         qorder = ((Tc - Ti)/Tco)**(1./4.)
c	 if (qorder .gt. qmax) qorder = qmax
         cplan = smax*Ti/(2.*Tco*qorder**2)
         alplan = vmax/Vi/(2.*Tco*qorder**2)
         betlan = vmax/Vi*vmax/(0.001*smax)/(2.*Tco*qorder**2)
	end if

        glan = smax*((Ti - Tc)*(qorder**2 - 1.) + 1./3.*Tco*(qorder**6 - 1.))
        slan = -smax*(qorder**2 - 1.)
        vlan = -vmax*(qorder**2 - 1.)

c	write(71,'(a8,i5,99f14.5)') 'landauqr',ispec,Pi,Ti,Tco,vmax,smax,qorder,Tc,x,glan/1000.,cplan

        return
        end
