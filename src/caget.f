	subroutine caget(func,a,b,ires)

C  Slight modification and re-write of routine zbrac from Numerical Recipes (Press et al., 1992)

	integer i,iter
	double precision fa,fb
	double precision func,a,b
	integer ires
	external func
	double precision, parameter :: gold=1.618
	integer, parameter :: itmax=100

	if (a .eq. b) stop 'No initial range in cage'
	ires = 1
	iter = 0
	fa = func(a)
	fb = func(b)
	do 1 i=1,itmax
	 if (fa*fb .lt. 0.) return
	 if (abs(fa) .lt. abs(fb)) then
	  a = a + gold*(a - b)
	  fa = func(a)
	 else
	  b = b + gold*(b - a)
	  fb = func(b)
	 end if
1	continue
	ires = 0

	return
	end
	 
