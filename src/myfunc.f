         subroutine myfunc(val, ndim, x, grad, need_gradient, f_data)

ctest	include 'P1'									! For testing purposes
ctest	include 'chem.inc'								! For testing purposes
	double precision func,f_data
          double precision val, x(ndim), grad(ndim)
          integer ndim, need_gradient,ncall,i
	  common /funcom/ ncall
	  data ncall/0/
	  ncall = ncall + 1

	  val = func(x)
          if (need_gradient.ne.0) then
	   call dfunc(x,grad)
          endif
ctest	call nform(x,n,n1,q2,nspec,ndim)						! For testing purposes
ctest	write(31,'(a8,i10,999f22.12)') 'myfunc',ncall,val,(n(i),i=1,nspec),x,grad	! For testing purposes
	 return
         end
