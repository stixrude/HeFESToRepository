           subroutine myconfeas(val, ndim, x, grad, need_gradient, d)
	    include 'P1'
	    include 'chem.inc'

	integer ndim,i
            integer need_gradient
            double precision val, x(ndim), grad(ndim), d(nspecp)

	    val = d(nspec+1)
	    do 1 i=1,ndim
1	    val = val + d(i)*x(i)
            if (need_gradient.ne.0) then
	     do 2 i=1,nspec
	      grad(i) = d(i)
2	     continue
            endif
c	write(31,'(a8,99e12.5)') 'myconstr',val,(x(i),i=1,ndim),grad
	    return
            end
