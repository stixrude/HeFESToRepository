        subroutine restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)

        include 'P1'
	include 'numpar.inc'

	integer nspec,nph,ispec,lph
        double precision f(nphasep,nspecp),dn(nspecp)
	double precision lags(nspecp)
        logical absents(nspecp),absent(nphasep)
        logical allow(nphasep),reslog,resall
        logical spinod(nspecp),spinph(nphasep)
        common /spinc/ spinod,spinph
	reslog = .false.

C  Restore absent species in present phases
        do 81 lph=1,nph
         if (absent(lph)) go to 81
	 if (.not. allow(lph)) go to 81
         do 82 ispec=1,nspec
          if (f(lph,ispec) .eq. 0.) go to 82
	  if (spinod(ispec)) go to 82
	  if (resall) go to 83
	  if (lags(ispec) .gt. 0.) go to 82
83	  continue
          if (absents(ispec)) then
           absents(ispec) = .false.
           dn(ispec) = tradd
	   print*, 'Restoring species',ispec
	   write(31,*) 'Restoring species',ispec
	   reslog = .true.
          end if
82       continue
81      continue

        return
        end
