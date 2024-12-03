         subroutine myvolw(val, ndim, x, grad, need_gradient, f_data)

	double precision f_data,pressurew
          double precision val, x(ndim), grad(ndim)
          integer ndim, need_gradient
	  val = pressurew(x)
          if (need_gradient.ne.0) then
c	   call dfunc(x,grad)
          endif
c	write(31,'(a5,99f12.5)') 'myvol',x,val,grad
	 return
         end
