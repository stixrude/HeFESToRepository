         subroutine myvol(val, ndim, x, grad, need_gradient, f_data)

	double precision f_data,pressure
          double precision val, x(ndim), grad(ndim)
          integer ndim, need_gradient
	  val = pressure(x)
          if (need_gradient.ne.0) then
	   call dfunc(x,grad)
          endif
c	write(31,'(a5,99f12.5)') 'myvol',x,val,grad
	 return
         end
