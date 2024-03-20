	function tlindeman(vol,wmagg,fnagg,thet)

C  Temperature via the Lindemann law 
C  Eq. (1) of Wolf and Jeanloz (1984) JGR

        include 'P1'
        include 'const.inc'

	double precision vol,wmagg,fnagg,thet,amass,amu,aspace,fac,hbar,tlindeman
	double precision, parameter :: flin = 0.1533, Angstrom = 1.e10

	aspace = (vol/fnagg/avn)**(1./3.)*0.01*Angstrom
	amass = wmagg/fnagg
	amu = 1.e-3/avn
	hbar = hplanck/2/pirad
	fac = (hbar)**2/boltzk/amu*Angstrom**2
	tlindeman = amass/fac*thet**2*aspace**2/9.*flin**2

	print*, fac,amass,aspace,thet,flin,tlindeman
	return
	end
