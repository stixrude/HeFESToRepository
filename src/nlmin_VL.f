            subroutine nlmin_VL(x,xlow,xupp,minf)
	include 'P1'
	include 'numpar.inc'

	double precision xlow,xupp
            external myvoll
             double precision  x(nspecp), minf
             integer ires
             integer*8 opt
             integer, parameter :: maxeval=500
             include 'nlopt.f'

c             call nlo_create(opt, NLOPT_GN_ISRES, 1)
             call nlo_create(opt, NLOPT_LN_COBYLA, 1)
             call nlo_set_min_objective(ires, opt, myvoll, 0)

c	call nlo_set_lower_bounds1(ires, opt, 1.e-15)
c	call nlo_set_upper_bounds1(ires, opt, 1.e+5)
	call nlo_set_lower_bounds1(ires, opt, xlow)
	call nlo_set_upper_bounds1(ires, opt, xupp)
     
c             call nlo_set_xtol_rel(ires, opt, 1.D-4)
c             call nlo_set_xtol_abs(ires, opt, ssmall)
             call nlo_set_ftol_abs(ires, opt, 1.E-10)
             call nlo_set_maxeval(ires, opt, maxeval)
     
             call nlo_optimize(ires, opt, x, minf)
             if (ires.lt.0) then
c               write(*,*) 'nlopt V failed!',ires
             else
c                write(*,*) 'found min at ', (x(i),i=1,nnull)
c                write(*,*) 'min val = ', minf
             endif

        if (ires .eq. 5) write(31,*) 'WARNING: Reached maximum evaluations in nlmin_V',maxeval
     
             call nlo_destroy(opt)
     
     	end
