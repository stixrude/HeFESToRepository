	subroutine cages(func,a,b,n,ac,bc,nc)

C  Slight modification and re-write of routine zbrak from Numerical Recipes (Press et al., 1992)

	double precision dx,fdn,fup,x
	integer n,nc,i,ncc
	double precision func,a,b,ac(nc),bc(nc)
	external func

	ncc = 0
	x = a
	dx = (b - a)/float(n)
	fdn = func(x)
	do 1 i=1,n
	 x = x + dx
	 fup = func(x)
	 if (fdn*fup .lt. 0.) then
	  ncc = ncc + 1
	  ac(ncc) = x - dx
	  bc(ncc) = x
	  if (ncc .eq. nc) go to 10
	 end if
	 fdn = fup
1	continue
10	continue
	nc = ncc
	return
	end

	
