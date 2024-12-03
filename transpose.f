	subroutine transpose(m,n,mp,np,a,at)

	integer m,n,mp,np,i,j,mt,nt
	double precision a(mp,np),at(np,mp)

	mt = n
	nt = m
	do 1 i=1,mt
	 do 1 j=1,nt
	  at(i,j) = a(j,i)
1	continue

	return
	end
	  
