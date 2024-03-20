	subroutine cage(func,a,b,clow,cupp,ires)

C  Modification of routine zbrac from Numerical Recipes (Press et al., 1992)
C  Add absolute lower and upper limits on brackets

	integer i,iter
	double precision clow,cupp,at,bt,fa,fb
	double precision func,a,b
	integer ires
	external func
	double precision, parameter :: gold=1.618
	integer, parameter :: itmax=100

	if (a .ge. b) stop 'No initial range in cage'
	if (a .lt. clow) stop 'initial guess violates lower bound'
	if (b .gt. cupp) stop 'initial guess violates upper bound'
	ires = 1
	iter = 0
	fa = func(a)
	fb = func(b)
c	write(*,*) a,b,fa,fb
	do 1 i=1,itmax
	 if (fa*fb .lt. 0.) return
	 if (abs(fa) .lt. abs(fb)) then
	  at = a + gold*(a - b)
	  if (at .lt. clow) then
	   ires = -1
	   return
	  end if
	  a = at
	  fa = func(a)
	 else
	  bt = b + gold*(b - a)
	  if (bt .gt. cupp) then
	   ires = -2
	   return
	  end if
	  b = bt
	  fb = func(b)
	 end if
c	 write(31,*) i,a,b,fa,fb
1	continue
	ires = 0

	return
	end
	 
