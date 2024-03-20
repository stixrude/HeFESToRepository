	subroutine bserch(xx,n,x,i)

C  Search the array xx for the value x.  Return the index i of the value xx(i) that occurs just before x.
C  Assume that the values in the array xx vary monotonically.

	integer i,j,k,n
	double precision x,xx(n)

	i = 0
	j = n + 1
	do 
	 k = (i + j)/2
	 if ((xx(n) .ge. xx(1)) .eqv. (xx(k) .ge. x)) then
	  j = k
	 else
	  i = k
	 end if
	 if (i + 1 .ge. j) exit
	end do
	if (x .eq. xx(1)) then
	 i = 1
	end if
	if (x .eq. xx(n)) then
	 i = n
	end if

	return
	end
