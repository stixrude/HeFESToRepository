         subroutine myvoll(val, ndim, x, grad, need_gradient, f_data)

	double precision f_data,pressurel
          double precision val, x(ndim), grad(ndim)
          integer ndim, need_gradient
	  val = pressurel(x)
          if (need_gradient.ne.0) then
	   call dfunc(x,grad)
          endif
c	write(31,'(99f12.5)') x,val,grad
	 return
         end
