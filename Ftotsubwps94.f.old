	subroutine Ftotsubw(ispec,Vi,Ftot)

C  Equation of state of Pitzer and Sterner (1994) JCP
C  Ti in units of Kelvin
C  rho in units of mol/cm^3
C  pressure in units of GPa
C  Includes ideal part of Saul and Wagner (1989)

        include 'P1'
        include 'const.inc'
	include 'water.inc'

        integer i,ispec
	double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
	double precision d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
	double precision e1,e2,e3,e4,e5,e6,e7,e8,e9,e10
	double precision rho,Ftot,fac,dbl,ee,Fig,vratio,Fxs,rhor,Tr,sum
	double precision den,deriv1,deriv2,dent
        double precision Pi,Ti
        double precision Vi,apar
	double precision volve,entve,cpve,bkve,CP,ent,K
        logical isochor
        common /volent/ volve,entve,cpve,bkve
        common /state/ apar(nspecp,nparp),Ti,Pi
        fac = hplanck/sqrt(2.*pirad*boltzk)
        ee = exp(1.)
        dbl = fac/sqrt(Ti*wmol/(1000.*avn))
        vratio = (Vi/avn*1.e-6)/dbl**3
	Fig = -log(vratio*ee)

	rho = 1./Vi

C  Ideal Part
	rhor = wmol*rho/rhocritical
	Tr = Tcritical/Ti
	sum = 0.
	do 1 i=4,8
	 sum = sum + a(i)*log(1. - exp(-g(i)*Tr))
1	continue
	Fig = log(rhor) + a(1) + a(2)*Tr + a(3)*log(Tr) + sum
	Fig = Rgas*Ti*Fig

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

C  Equation 1
	Fxs = c1*rho + 1./(c2 + c3*rho + c4*rho**2 + c5*rho**3 + c6*rho**4) - 1./c2 
     &  - (c7/c8)*(exp(-c8*rho) - 1.) - (c9/c10)*(exp(-c10*rho) - 1.)
	Ftot = Rgas*Ti*Fxs + Fig			! J/mol
c	Ftot = Fig

        volve = Vi
        entve = ent
        cpve = Cp
        bkve = K

	return
	end
