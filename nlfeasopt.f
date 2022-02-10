            subroutine nlfeasopt(x,nnew,ndim,minf,icase,itry,ires)
	include 'P1'
	include 'numpar.inc'
	include 'chem.inc'
	include 'absent.inc'
	include 'lag.inc'

	integer ndim,icase,i,ic,j,jconstr,jspec,nconstr,need_gradient,iph,ispec,jph,kph,itry
	double precision grad,val,vsum
            external myfeas, myconfeas
             double precision da(nspecp,nspecp)
             double precision x(nspecp), minf
	     double precision nnew(nspecp)
	     logical valid,validc
             integer ires
             integer*8 opt
             include 'nlopt.f'
	     integer, parameter :: itmax = 100
	double precision, parameter :: ftol=1.e-5, xtol=1.e-5
c	print*, 'nlmin',nconstr,ndim

             call nlo_create(opt, NLOPT_LD_SLSQP, ndim)
             call nlo_set_min_objective(ires, opt, myfeas, 0)

	        do 1 i=1,nspecp
                 do 1 j=1,nspecp
1                da(i,j) = 0.
     
        nconstr = 0
C---> Constrain the amount of each present species to be greater than zero
c        do 22 ispec=1,nspec
c         if (absents(ispec)) go to 22
c         nconstr = nconstr + 1
c         da(ispec,nconstr) = -1.0
c22      continue
c        go to 23
C<---
C  Set inequality constraints on N_jk (j components, k sites)
            valid = validc(vsum,x,da,nconstr)

        do 21 ispec=1,nspec
         if (absents(ispec)) go to 21
         iph = lphase(ispec)
         if (ispec  .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph)-1) then
          nconstr = nconstr+1
          da(ispec,nconstr) = -1.0
         end if
21      continue
23	continue
c       write(31,*) 'nlmin nconstr after additions = ',nconstr

c	write(31,*) 'After validc',vsum,nconstr
c	write(31,'(71f9.5)') (da(1,i),i=1,nconstr)
c	write(31,'(71f9.5)') (x(i),i=1,nspec)
       do 2 jconstr=1,nconstr
	 da(nspec+1,jconstr) = ssmall*10.
         call nlo_add_inequality_constraint(ires, opt,
     &     myconfeas, da(1,jconstr), 1.D-8)
2      continue

C  Set equality constraints for composition and absent n_i
	do 3 ic=1,nc
	 jconstr = jconstr + 1
	 do 4 jspec=1,nspec
	  da(jspec,jconstr) = s(ic,jspec)
4	 continue
	 da(nspec+1,jconstr) = -b(ic)
         call nlo_add_equality_constraint(ires, opt,
     &     myconfeas, da(1,jconstr), 1.D-8)
3	continue

c             call nlo_set_xtol_rel(ires, opt, 1.D-4)
c             call nlo_set_xtol_abs(ires, opt, ssmall)
             call nlo_set_ftol_abs(ires, opt, 1.E-10)
             call nlo_set_maxeval(ires, opt, 100)
     
             call nlo_optimize(ires, opt, x, minf)
             if (ires.lt.0) then
               write(31,*) 'WARNING: nlopt failed in nlfeas!',ires,(x(i),i=1,ndim)
             else
c                write(31,*) 'found min at ', (x(i),i=1,ndim)
c                write(31,*) 'min val = ', minf
             endif

	icase = 0
       do 5 jconstr=1,nconstr
	call myconfeas(val, ndim, x, grad, need_gradient, da(1,jconstr))
	if (val .gt. 1.e-8) then
	 write(31,*) 'WARNING: Constraint is violated',jconstr,val
	 icase = -1
	end if
5      continue

        write(31,*) 'Leaving nlfeasopt with n-vector:',minf,ires,icase
        write(31,'(71f9.5)') (n(i),i=1,nspec)

c	call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
     
             call nlo_destroy(opt)
     
	return
     	end
