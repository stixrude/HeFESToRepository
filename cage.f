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

	ires = 0
c	if (a .ge. b) stop 'No initial range in cage'
c	if (a .lt. clow) stop 'initial guess violates lower bound'
c	if (b .gt. cupp) stop 'initial guess violates upper bound'
	if (a .ge. b) then
c	 print*, 'No initial range in cage'
	 return
	end if
	if (a .lt. clow) then
c	 print*, 'initial guess violates lower bound'
	 return
	end if
	if (b .gt. cupp) then
c	 print*, 'initial guess violates upper bound'
	 return
	end if
	ires = 1
	iter = 0
c	print*, 'In cage',a,b,clow,cupp,ires
	fa = func(a)
	fb = func(b)
c	write(*,*) a,b,fa,fb
	at = a
	bt = b
	do 1 i=1,itmax
	 if (fa*fb .lt. 0.) return
	 if (abs(fa) .lt. abs(fb)) then
c	  at = a + gold*(a - b)
	  if (at .le. clow) then
	   ires = -1
	   return
	  end if
	  at = max(a + gold*(a - b),clow)
	  a = at
	  fa = func(a)
	 else
c	  bt = b + gold*(b - a)
	  if (bt .ge. cupp) then
	   ires = -2
	   return
	  end if
	  bt = min(b + gold*(b - a),cupp)
	  b = bt
	  fb = func(b)
	 end if
c	 write(31,*) i,a,b,fa,fb
c	 print*, 'cage',i,a,b,clow,cupp,fa,fb,ires
1	continue
	ires = 0

	return
	end
	 
