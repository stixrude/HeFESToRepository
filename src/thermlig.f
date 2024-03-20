        subroutine thermlig(ispec,Vi,Fig,Eig,Sig,Pig,Cvig,betaig,Kig,Kigp,aktig,daktdvig)
	
	include 'P1'
        include 'theory.inc'
	include 'const.inc'
	include 'elem.inc'

	integer ispec,ielem,jelem
	double precision fn,apar
	double precision dbl,ee,fac,fig
        double precision Pi,Pig,sig,Ti,vratio
	double precision Vi,aktig,daktdvig,Kig,Kigp,Eig
	double precision betaig,Cvig
        common /state/ apar(nspecp,nparp),Ti,Pi
        fac = hplanck/sqrt(2.*pirad*boltzk)
        ee = exp(1.)

	fn =   apar(ispec,1)

        Sig = 0.
        do 12 ielem=1,nelem
         if (sat(ielem,ispec) .eq. 0.) go to 12
         dbl = fac/sqrt(Ti*wat(ielem)/(1000.*avn))
         vratio = (Vi/avn*1.e-6)/(sat(ielem,ispec)*dbl**3)
         Sig = Sig + sat(ielem,ispec)*(log(vratio) + 2.5)
12      continue
        Sig = Rgas*Sig					! J/mol/K

        Cvig = fn*1.5*Rgas				! J/mol/K

	Eig = Cvig*Ti					! J/mol

        Fig = 0.
        do 11 ielem=1,nelem
         if (sat(ielem,ispec) .eq. 0.) go to 11
         dbl = fac/sqrt(Ti*wat(ielem)/(1000.*avn))
         vratio = (Vi/avn*1.e-6)/(sat(ielem,ispec)*dbl**3)
         Fig = Fig + sat(ielem,ispec)*log(vratio*ee)
11      continue
        Fig = -Rgas*Ti*Fig				! J/mol

        betaig = 0.0

	Pig = 0.001*fn*Rgas*Ti/Vi			! GPa

	aktig = 0.001*fn*Rgas/Vi			! GPa/K
	daktdvig = -0.001*fn*Rgas/(Vi*Vi)		! GPa/K/(cm^3/mol)
	Kig = Pig
	Kigp = 1.0

	return
	end
