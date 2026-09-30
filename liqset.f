        subroutine liqset(ispec,apar,aliq,fn,wm,To,Fo,Vo,Telo,eta,Telinf,zelo,xi,bpar,fmth)

C  Assign liquid state parameter including those in apar and the a_ij in aliqc to scalar variables and to aliq
C  Also compute bounds in volume search

        include 'P1'
	include 'const.inc'

	integer ispec,i,j,nbm,nobm,noth,mfit,maxij,lineart,noln
	double precision fn,wm,to,fo,vo,telo,eta,telinf,zelo,xi,bpar
        double precision Ko,Kop,Kopp,So,CVo,betao,akto,daktdvo,fmth,htl
	double precision flow,fupp,par,parold,swap,vlow,vsp,vupp,vsplow,vspupp,detsp,asp,bsp
	double precision f1,f2,v1,v2,c,fextremum,vextremum
        double precision apar(nspecp,nparp)
        double precision aliq(nparp,nparp),aliqc
	double precision, parameter :: vsquaredminimum = 0.1
	double precision, parameter :: Tsmall=1.e-5
        double precision, parameter :: vsmall = 1.d-6    !  Used for setting valid volume domain
        common /liqc/ aliqc(nspecp,nparp,nparp),mfit,nobm,noth,maxij,lineart,nbm,noln

        do 21 i=1,nparp
         do 21 j=1,nparp
          aliq(i,j) = aliqc(ispec,i,j)
21      continue

        fn =      apar(ispec,1)		!	Atoms in formula unit
        wm =      apar(ispec,3)		!	Formula mass (g/mol)
        To =      apar(ispec,4)		!	T_0 (K)
        Fo =      apar(ispec,5)		!	F_0 (kJ/mol)
        Vo =      apar(ispec,6)		!	V_0 (cm^3/mol)
        Ko =      apar(ispec,7)		!	K_0 (GPa)
        Kop =     apar(ispec,8)		!	K_0 prime (-)
        Kopp =    apar(ispec,9)		!	K_0K_0 prime prime (zero for third order)
        So =      apar(ispec,10)	!	S_0 (J/mol/K)
        CVo =     apar(ispec,11)	!	C_V0 (J/mol/K)
        betao =   apar(ispec,12)	!	beta_0 (J/(mol K cm^3/mol))
	Telo =    apar(ispec,13)	!	T_el0-T_elinf (K)
	eta =     apar(ispec,14)	!	-eta (-)
        Telinf =  apar(ispec,15)	!	T_elinf (K)
	akto =    apar(ispec,26)	!	(alpha K_T)_0 (MPa/K)
	daktdvo = apar(ispec,27)	!	(d alpha K_T / d V)_0 (MPa/(K cm^3/mol))
        zelo =    apar(ispec,28)	!	z_el0 (mJ/mol/K)
	bpar =    apar(ispec,29)	!	Volume dependence of z_el now assumed to be zero
	xi =      apar(ispec,30)	!	Volume dependence of z_el now assumed to be zero
	htl =     apar(ispec,31)	!	Material index (1=liquid, 4=vapor)
	fmth =    apar(ispec,33)	!	Rosenfeld-Tarazona exponent (-)

	if (To .lt. 0.) then
	 To = -Tsmall
	end if

C  Computed quantities

C  Volume limits set by real vibrational frequency
C  Find roots of Eq. 41 SLB05
C  Identify domain of positive v^2
	apar(ispec,51) = 0.
	apar(ispec,52) = 0.
	if (htl .ne. 0.) then
C  Not a solid.  Eq. 41 does not apply
	 go to 11
	end if
11	continue

C  Spinodal instabilities at T=T_0
C  Find roots of Eq. 32 SLB05
C  Identify domain of positive bulk modulus
	c = 1.0
        bsp = (3.*Kop - 5.)
        asp = 27./2.*(Kop - 4.)
	detsp = bsp*bsp - 4.*asp*c
	vspupp = 1.d+15
	vsplow = 1.d-15
	if (detsp .lt. 0.) then
C  No roots: domain of positive bulk modulus is unbounded (does not occur for BM3).
	 go to 20
	end if
	if (asp .eq. 0.) then
C  Only one root, which has a negative value (f1=-1/7), and therefore corresponds to an upper bound on the volume (K_0'=4)
	 f1 = -c/bsp
	 vspupp = Vo*(2.*f1 + 1.)**(-3./2.)
	 go to 20
	end if
C  Two roots 
	f1 = (-bsp - sqrt(detsp))/(2.*asp)
	f2 = (-bsp + sqrt(detsp))/(2.*asp)
	if (max(f1,f2) .lt. 0.) then
C  Only an upper bound (Kop>4)
	 vspupp = Vo*(2.*max(f1,f2) + 1.)**(-3./2.)
	 go to 20
	end if
	if (min(f1,f2) .gt. 0.) then
C  Only a lower bound (does not occur for BM3).
	 vsplow = Vo*(2.*min(f1,f2) + 1.)**(-3./2.)
	 go to 20
	end if
C  One positive root and one negative root: upper and lower bounds (Kop<4)
	vspupp = Vo*(2.*min(f1,f2) + 1.)**(-3./2.)
	vsplow = Vo*(2.*max(f1,f2) + 1.)**(-3./2.)
20	continue
	apar(ispec,53) = max(vsplow,Vo/10.) + vsmall
	apar(ispec,54) = min(vspupp,Vo*10.) - vsmall
C  Relax upper bound on volume for liquids
c	if (htl .eq. 1) apar(ispec,54) = Vo*1000.
c	print*, 'spinodal limits: f1,f2,vsplow,vspupp,apar(ispec,51),apar(ispec,52),fextremum,vextremum'
c     &   ,f1,f2,vlow,vupp,apar(ispec,51),apar(ispec,52),fextremum,vextremum
c	write(31,'(a34,i5,99f12.5)') 'V bounds: vibrational and spinodal'
c     &   ,ispec,apar(ispec,51),apar(ispec,52),apar(ispec,53),apar(ispec,54),f1,f2
cc     &   ,a,b,det,vlow,vupp,asp,bsp,detsp,vsplow,vspupp

        return
100     format(a19,4f13.5)
        end
