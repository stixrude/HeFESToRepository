	double precision function pressurew(Vi)

C  Equation of state of Pitzer and Sterner (1994) JCP
C  Ti in units of Kelvin
C  rho in units of mol/cm^3
C  pressure in units of GPa

        include 'P1'
        include 'const.inc'
	include 'water.inc'

        integer i,ispec
	double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
	double precision rho,pressure
        double precision Pi,Ti
        double precision Vi,apar
        logical isochor
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chor/ ispec,isochor

	rho = 1./Vi

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

C  Equation 2
	pressure = rho + c1*rho**2 
     &  - rho**2*((c3 + 2.*c4*rho + 3.*c5*rho**2 + 4.*c6*rho**3)/(c2 + c3*rho + c4*rho**2 + c5*rho**3 + c6*rho**4)**2)
     &  + c7*rho**2*exp(-c8*rho) + c9*rho**2*exp(-c10*rho)

C  Units: GPa
	pressure = Rgas*Ti*pressure/1000.

	pressurew = pressure - Pi

c	print*, Ti,Pi,Vi,pressure,pressurew

	return
	end
