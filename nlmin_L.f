            subroutine nlmin_L(x,ndim,minf,iter,ires)
	include 'P1'
	include 'numpar.inc'
	include 'chem.inc'
	include 'lag.inc'
	include 'absent.inc'

	integer ndim,i,j,jconstr,nconstr,iph,ispec
	double precision dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,phugo,vtarg,starg,tlast,vhugo,vsum,wmagg
            external myfunc, myconstraint
             double precision da(nspecp,nspecp)
             double precision x(nspecp), minf
	     double precision lb(nspecp), ub(nspecp),palb(nspecp),paub(nspecp)
	     double precision Tspin,Psmall,Plarge,apar,Pi,Ti
	     logical valid,validc
	     logical chcalc,adcalc,hucalc,pabounds
             integer ires,maxeval,ncall,ncallold,iter,nvet,nvep
             integer*8 opt,local_opt
	     parameter (maxeval=1000)
	     double precision, parameter :: Tsmall = 1.e-5, eps=1.e-8
             common /funcom/ ncall
        common /nlbounds/ palb,paub,pabounds
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        common /state/ apar(nspecp,nparp),Ti,Pi
	     data ncallold/0/
             include 'nlopt.f'

	opt = 0
	local_opt = 0
             call nlo_create(opt, NLOPT_LD_SLSQP, ndim)
c             call nlo_create(opt, NLOPT_LD_AUGLAG, ndim)
c             call nlo_create(local_opt, NLOPT_LD_SLSQP, ndim)
c	     call nlo_set_local_optimizer(ires, opt, local_opt)
             call nlo_set_min_objective(ires, opt, myfunc, 0)
c             if (adcalc) call nlo_set_max_objective(ires, opt, myfunc, 0)
c        write(31,*) 'nlmin ndim =',ndim

C  Bounds
	 if (adcalc .or. chcalc) then
             call nlo_get_lower_bounds(ires, opt, lb)
             call nlo_get_upper_bounds(ires, opt, ub)
	     call Prange(Psmall,Plarge)
	     if (adcalc) lb(nnull+nvet) = Tsmall
	     if (adcalc .and. Tspin(ncall) .ne. 0.) ub(nnull+nvet) = Tspin(ncall)
	     if (chcalc) lb(nnull+nvep) = Psmall
	     if (chcalc) ub(nnull+nvep) = Plarge
             call nlo_set_lower_bounds(ires, opt, lb)
             call nlo_set_upper_bounds(ires, opt, ub)
	end if
	if (pabounds .and. ndim .gt. 0) then
         call nlo_get_lower_bounds(ires, opt, lb)
         call nlo_get_upper_bounds(ires, opt, ub)
	 do 11 i=1,ndim
	  lb(i) = palb(i) - eps
	  ub(i) = paub(i) + eps
c	  write(31,*) 'pabounds',i,lb(i),ub(i)
11	 continue
         call nlo_set_lower_bounds(ires, opt, lb)
         call nlo_set_upper_bounds(ires, opt, ub)
	end if
	  
c	write(31,*) 'T Bounds',lb(ndim),ub(ndim),Tspin(ncall),ndim,nvet,adcalc

	        do 1 i=1,nspecp
                 do 1 j=1,nspecp
1                da(i,j) = 0.
     
	nconstr = 0
C---> Constrain the amount of each present species to be greater than zero
	do 22 ispec=1,nspec
	 if (absents(ispec)) go to 22
	 nconstr = nconstr + 1
	 da(ispec,nconstr) = -1.0
22	continue
	go to 23
C<---
C--->  Constrain the total amount of each component on each site to be greater than zero
            valid = validc(vsum,x,da,nconstr)
c	write(31,*) 'nlmin nconstr after validc = ',nconstr

C  Add positive definitiveness constraint to each member of pairs of species with identical formulae in the same phase
C  This is necessary to prevent a large negative amount of one member of the pair.
	do 21 ispec=1,nspec
	 if (absents(ispec)) go to 21
	 iph = lphase(ispec)
         if (ispec  .ge. iophase(iph) .and. ispec .le. iophase(iph)+mophase(iph)-1) then
	  nconstr = nconstr+1
	  da(ispec,nconstr) = -1.0
	 end if
21	continue
C<---
23	continue
c	write(31,*) 'nlmin nconstr after additions = ',nconstr

       do 2 jconstr=1,nconstr
         call nlo_add_inequality_constraint(ires, opt,
     &     myconstraint, da(1,jconstr), 1.D-12)
2      continue

	    if (adcalc) then
	     da(nspec+nvet,nconstr+nvet) = 1.0
         call nlo_add_equality_constraint(ires, opt,
     &     myconstraint, da(1,nconstr+nvet), 1.D-5)
	    end if

	    if (chcalc) then
	     da(nspec+nvep,nconstr+nvep) = 1.0
         call nlo_add_equality_constraint(ires, opt,
     &     myconstraint, da(1,nconstr+nvep), 1.D-10)
	    end if

c             call nlo_set_xtol_rel(ires, opt, 1.D-4)
c             call nlo_set_xtol_abs(ires, opt, ssmall)
             call nlo_set_ftol_abs(ires, opt, 1.D-15)
             call nlo_set_maxeval(ires, opt, maxeval)
c             call nlo_set_ftol_abs(ires, local_opt, 1.D-15)
c             call nlo_set_maxeval(ires, local_opt, maxeval)
     
             call nlo_optimize(ires, opt, x, minf)
             if (ires.lt.0) then
c               write(31,'(a22,2i5,99f12.5)') 'WARNING: nlopt failed!',ires,iter,Pi,Ti,vtarg,starg
c     &          ,(x(i),i=1,ndim),(lb(i),i=1,ndim),(ub(i),i=1,ndim)
               write(31,*) 'WARNING: nlopt failed!',ires,iter,Pi,Ti,vtarg,starg
     &          ,(x(i),i=1,ndim),(lb(i),i=1,ndim),(ub(i),i=1,ndim)
             else
c                write(31,*) 'nlopt found min at ', (x(i),i=1,ndim)
c                write(31,*) 'nlopt min val = ', minf
             endif

	iter = ncall - ncallold
	if (iter .eq. maxeval) write(31,*) 'WARNING: Reached maximum evaluations in nlmin_L',maxeval,iter
c	write(31,*) 'myfunc in nlmin_L',ncall,ncallold,iter
c	write(95,*) 'Optimized',iter,minf
	ncallold = ncall
     
             call nlo_destroy(opt)
c             call nlo_destroy(local_opt)
     
	return
     	end
