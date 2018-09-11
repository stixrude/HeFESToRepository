	double precision function dfac(n)

C  compute the factorial of n
C  Modified after routine factrl (Numerical Recipes, Press et al., 1992)
C  to use intrinsic gammln function
	
	double precision x
	integer i,n,nmax
	integer, parameter :: ndim=33
	double precision s(ndim)
	save nmax,s
	data nmax,s(1)/0,1./

	if (n .lt. 0) then
	 print*, 'Error in computation of factorial dfac: negative argument'
	end if

	x = real(n)

	if (n .le. nmax) then
	 dfac = s(n+1)
	else if (n .le. ndim-1) then
	 do 1 i=nmax+1,n
	  s(i+1) = i*s(i)
1	 continue
	 nmax = n
	 dfac = s(n+1)
	else
	 dfac = exp(log_gamma(x+1.))
	end if

	return
	end

