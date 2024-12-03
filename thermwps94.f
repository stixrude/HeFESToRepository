	subroutine thermw(ispec,Vi,volnl,CP,CV,gamma,K,KS,alp,Ftot,ent,pressure,E)

C  Equation of state of Pitzer and Sterner (1994) JCP
C  Ti in units of Kelvin
C  rho in units of mol/cm^3
C  pressure in units of GPa

        include 'P1'
        include 'const.inc'
	include 'water.inc'
	include 'theory.inc'

	double precision, parameter :: Sconst = 86.808 -  19.103896651846917	!	J/mol/K to recover JANAF value at 1 bar 373 K
        integer i,ispec
	double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
	double precision d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
	double precision e1,e2,e3,e4,e5,e6,e7,e8,e9,e10
	double precision rho,Ftot,fac,dbl,ee,Fig,vratio,K,den,deriv1,deriv2,Tr,rhor,sum
	double precision Sig,Fxs,S,ent,CVig,CV,Sxs,f,ft,ftt,Kig,akt,pt,p,dent,pressure,gamma
	double precision CP,KS,alp,E,cs
        double precision Pi,Ti,volnl
        double precision Vi,apar
        logical isochor
        common /state/ apar(nspecp,nparp),Ti,Pi
        fac = hplanck/sqrt(2.*pirad*boltzk)
        ee = exp(1.)
        dbl = fac/sqrt(Ti*wmol/(1000.*avn))
        vratio = (Vi/avn*1.e-6)/dbl**3
        Fig = -log(vratio*ee)
	Sig = log(vratio) + 2.5
c	print*, 'De Broglie thermal wavelength (A)',dbl*1e10,fac,Ti,wmol,avn,hplanck,pirad,boltzk
c	print*, 'Volume ratio',vratio
c	print*, 'Ideal free energy, entropy',Fig,Sig

	rho = 1./Vi

C  Ideal Part
        rhor = wmol*rho*1000./rhocritical
        Tr = Tcritical/Ti
        sum = 0.
        do 1 i=4,8
         sum = sum + a(i)*log(1. - exp(-g(i)*Tr))
1       continue
        Fig = log(rhor) + a(1) + a(2)*Tr + a(3)*log(Tr) + sum
c	print*, 'Fig reduced',Fig
        Fig = Rgas*Ti*Fig
c	print*, 'Fig',Fig

	sum = 0.
	do 2 i=4,8
	 sum = sum + a(i)*g(i)/(exp(g(i)*Tr) - 1.)
2	continue
	Sig = a(2) + a(3)/Tr + sum
c	print*, 'Sig reduced',Sig
	Sig = Rgas*Tr*Sig - Fig/Ti
c	print*, 'Sig',Sig

	sum = 0.
	do 3 i=4,8
	 sum = sum + a(i)*g(i)**2*exp(-g(i)*Tr)/(1. - exp(-g(i)*Tr))**2
3	continue
	CVig = -a(3)/Tr**2 - sum
c	print*, 'CVig reduced',CVig
	CVig = -Tr**2*Rgas*CVig
c	print*, 'CVig',CVig

C  Equation 4
	c1  = c1a(1)*Ti**(-4)  + c1a(2)*Ti**(-2)  + c1a(3)*Ti**(-1)  + c1a(4)  + c1a(5)*Ti  + c1a(6)*Ti**2
	c2  = c2a(1)*Ti**(-4)  + c2a(2)*Ti**(-2)  + c2a(3)*Ti**(-1)  + c2a(4)  + c2a(5)*Ti  + c2a(6)*Ti**2
	c3  = c3a(1)*Ti**(-4)  + c3a(2)*Ti**(-2)  + c3a(3)*Ti**(-1)  + c3a(4)  + c3a(5)*Ti  + c3a(6)*Ti**2
	c4  = c4a(1)*Ti**(-4)  + c4a(2)*Ti**(-2)  + c4a(3)*Ti**(-1)  + c4a(4)  + c4a(5)*Ti  + c4a(6)*Ti**2
	c5  = c5a(1)*Ti**(-4)  + c5a(2)*Ti**(-2)  + c5a(3)*Ti**(-1)  + c5a(4)  + c5a(5)*Ti  + c5a(6)*Ti**2
	c6  = c6a(1)*Ti**(-4)  + c6a(2)*Ti**(-2)  + c6a(3)*Ti**(-1)  + c6a(4)  + c6a(5)*Ti  + c6a(6)*Ti**2
	c7  = c7a(1)*Ti**(-4)  + c7a(2)*Ti**(-2)  + c7a(3)*Ti**(-1)  + c7a(4)  + c7a(5)*Ti  + c7a(6)*Ti**2
	c8  = c8a(1)*Ti**(-4)  + c8a(2)*Ti**(-2)  + c8a(3)*Ti**(-1)  + c8a(4)  + c8a(5)*Ti  + c8a(6)*Ti**2
	c9  = c9a(1)*Ti**(-4)  + c9a(2)*Ti**(-2)  + c9a(3)*Ti**(-1)  + c9a(4)  + c9a(5)*Ti  + c9a(6)*Ti**2
	c10 = c10a(1)*Ti**(-4) + c10a(2)*Ti**(-2) + c10a(3)*Ti**(-1) + c10a(4) + c10a(5)*Ti + c10a(6)*Ti**2

C  Temperature derivative of Equation 4
	d1  = -4.*c1a(1)*Ti**(-5)  + -2.*c1a(2)*Ti**(-3)  + -c1a(3)*Ti**(-2)            + c1a(5)  + 2.*c1a(6)*Ti
	d2  = -4.*c2a(1)*Ti**(-5)  + -2.*c2a(2)*Ti**(-3)  + -c2a(3)*Ti**(-2)            + c2a(5)  + 2.*c2a(6)*Ti
	d3  = -4.*c3a(1)*Ti**(-5)  + -2.*c3a(2)*Ti**(-3)  + -c3a(3)*Ti**(-2)            + c3a(5)  + 2.*c3a(6)*Ti
	d4  = -4.*c4a(1)*Ti**(-5)  + -2.*c4a(2)*Ti**(-3)  + -c4a(3)*Ti**(-2)            + c4a(5)  + 2.*c4a(6)*Ti
	d5  = -4.*c5a(1)*Ti**(-5)  + -2.*c5a(2)*Ti**(-3)  + -c5a(3)*Ti**(-2)            + c5a(5)  + 2.*c5a(6)*Ti
	d6  = -4.*c6a(1)*Ti**(-5)  + -2.*c6a(2)*Ti**(-3)  + -c6a(3)*Ti**(-2)            + c6a(5)  + 2.*c6a(6)*Ti
	d7  = -4.*c7a(1)*Ti**(-5)  + -2.*c7a(2)*Ti**(-3)  + -c7a(3)*Ti**(-2)            + c7a(5)  + 2.*c7a(6)*Ti
	d8  = -4.*c8a(1)*Ti**(-5)  + -2.*c8a(2)*Ti**(-3)  + -c8a(3)*Ti**(-2)            + c8a(5)  + 2.*c8a(6)*Ti
	d9  = -4.*c9a(1)*Ti**(-5)  + -2.*c9a(2)*Ti**(-3)  + -c9a(3)*Ti**(-2)            + c9a(5)  + 2.*c9a(6)*Ti
	d10 = -4.*c10a(1)*Ti**(-5) + -2.*c10a(2)*Ti**(-3) + -c10a(3)*Ti**(-2)           + c10a(5) + 2.*c10a(6)*Ti

C  Second temperature derivative of Equation 4
	e1  = +20.*c1a(1)*Ti**(-6)  + +6.*c1a(2)*Ti**(-4)  + +2.*c1a(3)*Ti**(-3)                  + 2.*c1a(6)
	e2  = +20.*c2a(1)*Ti**(-6)  + +6.*c2a(2)*Ti**(-4)  + +2.*c2a(3)*Ti**(-3)                  + 2.*c2a(6)
	e3  = +20.*c3a(1)*Ti**(-6)  + +6.*c3a(2)*Ti**(-4)  + +2.*c3a(3)*Ti**(-3)                  + 2.*c3a(6)
	e4  = +20.*c4a(1)*Ti**(-6)  + +6.*c4a(2)*Ti**(-4)  + +2.*c4a(3)*Ti**(-3)                  + 2.*c4a(6)
	e5  = +20.*c5a(1)*Ti**(-6)  + +6.*c5a(2)*Ti**(-4)  + +2.*c5a(3)*Ti**(-3)                  + 2.*c5a(6)
	e6  = +20.*c6a(1)*Ti**(-6)  + +6.*c6a(2)*Ti**(-4)  + +2.*c6a(3)*Ti**(-3)                  + 2.*c6a(6)
	e7  = +20.*c7a(1)*Ti**(-6)  + +6.*c7a(2)*Ti**(-4)  + +2.*c7a(3)*Ti**(-3)                  + 2.*c7a(6)
	e8  = +20.*c8a(1)*Ti**(-6)  + +6.*c8a(2)*Ti**(-4)  + +2.*c8a(3)*Ti**(-3)                  + 2.*c8a(6)
	e9  = +20.*c9a(1)*Ti**(-6)  + +6.*c9a(2)*Ti**(-4)  + +2.*c9a(3)*Ti**(-3)                  + 2.*c9a(6)
	e10 = +20.*c10a(1)*Ti**(-6) + +6.*c10a(2)*Ti**(-4) + +2.*c10a(3)*Ti**(-3)                 + 2.*c10a(6)

	den = c2 + c3*rho + c4*rho**2 + c5*rho**3 + c6*rho**4
	deriv1 = c3 + 2.*c4*rho + 3.*c5*rho**2 + 4.*c6*rho**3
	deriv2 = 2.*c4 + 6.*c5*rho + 12.*c6*rho**2
	dent = d2 + d3*rho + d4*rho**2 + d5*rho**3 + d6*rho**4

C  Helmholtz free energy Equation 1.  f = F/RT.  F in J/mol
        f = c1*rho + 1./den - 1./c2
     &  - (c7/c8)*(exp(-c8*rho) - 1.) - (c9/c10)*(exp(-c10*rho) - 1.)
	Ftot = Rgas*Ti*f + Fig - Ti*Sconst
c	print*, 'f=F/RT',f

C  Entropy.  ent in J/mol/K
	ft = d1*rho - ((d2 + d3*rho + d4*rho**2 + d5*rho**3 + d6*rho**4)/den**2 - d2/c2**2) 
     &  + (c7*d8/c8*rho   + c7*d8/c8**2   - d7/c8)*exp(-c8*rho)
     &  + (c9*d10/c10*rho + c9*d10/c10**2 - d9/c10)*exp(-c10*rho)
     &  + d7/c8 - c7*d8/c8**2 + d9/c10 - c9*d10/c10**2
c	print*, 'ft=df/dT',ft,-f*Ti/Tr
	S = -Rgas*Ti*ft - Rgas*f + Sig
	ent = S + Sconst

C  Heat Capacity.  CV in J/mol/K
	ftt = e1*rho - (e2 + e3*rho + e4*rho**2 + e5*rho**3 + e6*rho**4)/den**2
     &  + 2.*(d2 + d3*rho + d4*rho**2 + d5*rho**3 + d6*rho**4)**2/den**3
     &  + e2/c2**2 - 2.*d2**2/c2**3
     &  + (d7*d8/c8*rho + c7*e8/c8*rho - c7*d8**2/c8**2*rho + d7*d8/c8**2 + c7*e8/c8**2 - 2.*c7*d8**2/c8**3 
     &  - e7/c8 + d7*d8/c8**2)*exp(-c8*rho)
     &  - (c7*d8/c8*rho   + c7*d8/c8**2   - d7/c8)*d8*rho*exp(-c8*rho)
     &  + (d9*d10/c10*rho + c9*e10/c10*rho - c9*d10**2/c10**2*rho + d9*d10/c10**2 + c9*e10/c10**2 - 2.*c9*d10**2/c10**3 
     &  - e9/c10 + d9*d10/c10**2)*exp(-c10*rho)
     &  - (c9*d10/c10*rho   + c9*d10/c10**2   - d9/c10)*d10*rho*exp(-c10*rho)
     &  + e7/c8 - d7*d8/c8**2 - d7*d8/c8**2 - c7*e8/c8**2 + 2.*c7*d8**2/c8**3
     &  + e9/c10 - d9*d10/c10**2 - d9*d10/c10**2 - c9*e10/c10**2 + 2.*c9*d10**2/c10**3
c	print*, 'ftt=d2f/dT2',ftt
	CV = -Rgas*Ti**2*ftt - 2.*Rgas*Ti*ft + CVig

C  Bulk modulus.  K in GPa
	K = 1. + 2.*c1*rho - 2.*rho*(deriv1/den**2) - rho**2*(deriv2/den**2) + 2.*rho**2*(deriv1**2/den**3) 
     &  + 2.*c7*rho*exp(-c8*rho) - c7*c8*rho**2*exp(-c8*rho) + 2.*c9*rho*exp(-c10*rho) - c9*c10*rho**2*exp(-c10*rho)
	K = K*Rgas*Ti/Vi/1000.
 
C  Pressure Equation 2. p=P/RT.  pressure in GPa
        p = rho + c1*rho**2
     &  - rho**2*((c3 + 2.*c4*rho + 3.*c5*rho**2 + 4.*c6*rho**3)/den**2)
     &  + c7*rho**2*exp(-c8*rho) + c9*rho**2*exp(-c10*rho)
	pressure = Rgas*Ti*p/1000.

C  Thermal pressure coefficient alpha*K_T.  akt in MPa/K
	pt = d1*rho**2 - rho**2*((d3 + 2.*d4*rho + 3.*d5*rho**2 + 4.*d6*rho**3)/den**2 
     &  - 2.*(c3 + 2.*c4*rho + 3.*c5*rho**2 + 4.*c6*rho**3)*dent/den**3)
     &  + d7*rho**2*exp(-c8*rho) - c7*d8*rho**3*exp(-c8*rho)
     &  + d9*rho**2*exp(-c10*rho) - c9*d10*rho**3*exp(-c10*rho)
c	print*, 'p,pt',p,pt
	akt = Rgas*Ti*pt + 1000.*pressure/Ti

C  Internal energy
	E = F + Ti*S

C  Thermal expasivity
	alp = 0.001*akt/K

C  Gruneisen parameter
	gamma = akt*Vi/CV

C  Isobaric heat capacity
	CP = CV*(1. + alp*gamma*Ti)

C  Adiabatic bulk modulus
	KS = K*(1. + alp*gamma*Ti)

C  Sound speed
	cs = sqrt(KS*Vi/wmol)

	return
	end
