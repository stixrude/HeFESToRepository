	subroutine icebcc(Vi,un,pn,kn,kpn)
C  Following frenchredmer_15 add zero-point nuclear contribution to VII-VII'-X transitions 
C  that captures anomalies in pressure and bulk modulus.
C  Use a simpler functional form that reproduces their result closely:
C  un(rho) = u7(rho) + f(rho)*(u10(rho) - u7(rho))
C  where u7 and u10 are the internal energies of VII and X respectively
C  and the Fermi function f(rho) = [1. + exp(-(rho-rho0)/delta))]^-1
C  Assume that u7(rho) = 0.  Then
C  un(rho) = f(rho)*u10(rho)
C  Further assume a linear form for u10(rho) = a + b*rho
	double precision Vi, rho, un, pn, kn, kpn, xmcapw
	double precision u, up, upp, uppp, f, fp, fpp, fppp
	data xmcapw / 0.018015268d0 /
c	double precision, parameter :: rho0 = 3.0, delta = 0.20, ua = -1.0, ub = 0. ! rho, delta in g/cm^3, ua in kJ/g, ub in kJ*cm^3/g^2 = GPa/g^2
	double precision, parameter :: rho0 = 2.9, delta = 0.15, ua = -0.8, ub = 0. ! rho, delta in g/cm^3, ua in kJ/g, ub in kJ*cm^3/g^2 = GPa/g^2
c	double precision, parameter :: rho0 = 2.9, delta = 0.15, ua = -0.0, ub = 0. ! rho, delta in g/cm^3, ua in kJ/g, ub in kJ*cm^3/g^2 = GPa/g^2

	rho = 1000.*xmcapw/Vi
C  Internal energy of Ice X
	u = ua + ub*(rho - rho0)
c	u = ua + ub*log(rho/rho0)
C  Fermi function
	f = 1./(1. + exp(-(rho - rho0)/delta))
C  Derivatives
	up = ub 
	upp = 0.
	uppp = 0.
c	up = ub/rho
c	upp = -ub/rho**2
c	uppp = 2.*ub/rho**3
	fp = f*(1. - f)/delta
	fpp = f*(1. - 2.*f)*(1. - f)/delta**2
	fppp = f*(1. - f)/delta**3*(1. - 6.*f + 6.*f**2)
C  Thermodynamic quantities
	un = f*u
	pn = rho**2*(fp*u + f*up)
	kn = 2.*pn + rho**3*(fpp*u + 2.*fp*up + f*upp)
	kpn = 5. - 6.*pn/kn + rho**4/kn*(fppp*u + 3.*fpp*up + 3.*fp*upp + f*uppp)

        un = un*1000.*1000.*xmcapw                      ! J/mol

	return
	end
