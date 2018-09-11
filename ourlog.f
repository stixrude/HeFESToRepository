	double precision function ourlog(y)

	double precision y,x
c	double precision, parameter :: small=1.d-10
	double precision, parameter :: small=3.d-308

  	x = abs(y)
	if (x .le. small) then
  	 ourlog = log(small)
c	 ourlog = log(small) - 0.5 + 0.5*(x/small)**2
c	 ourlog = log(small) - 0.75 + (x/small)**2 - 0.25*(x/small)**4
c	 ourlog = log(small) - 1./small*(small - x)
	else
	 ourlog = log(x)
	end if

	if (y .le. -small) then
c	 write(31,*) 'WARNING: Negative Argument of logarithm',y,ourlog
	else
c	 if (y .le. small) write(31,*) 'WARNING: Small Argument of logarithm',y,ourlog
	end if

	return
	end
