         subroutine myfunc(val, ndim, x, grad, need_gradient, f_data)

	double precision func,f_data
          double precision val, x(ndim), grad(ndim)
          integer ndim, need_gradient,ncall
	  common /funcom/ ncall
	  data ncall/0/
	  ncall = ncall + 1

c	write(95,*) 'MyFunc'
	  val = func(x)
          if (need_gradient.ne.0) then
	   call dfunc(x,grad)
          endif
c	write(31,'(a8,i10,99f22.12)') 'myfunc',ncall,val,x,grad
c	print '(i5,99f12.5)', ncall,x,val,grad
	 return
         end
