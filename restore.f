        subroutine restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)

        include 'P1'
	include 'numpar.inc'

C  Add absent species in present phases:
C  resall = false.  Restore only those species that have negative (favorable) chemical affinities (lags)
C  resall = true.   Restore also those species that have positive chemical affinities up to lagsmax.

	integer nspec,nph,ispec,lph
	double precision, parameter :: lagsmax = 10.
        double precision f(nphasep,nspecp),dn(nspecp)
	double precision lags(nspecp)
        logical absents(nspecp),absent(nphasep)
        logical allow(nphasep),reslog,resall
        logical spinod(nspecp),spinph(nphasep)
        character*80 phname(nphasep),sname(nspecp)
        common /names/ phname,sname
        common /spinc/ spinod,spinph
	reslog = .false.
c	print '(a4,99f12.5)', 'lags',(lags(ispec),ispec=1,nspec)
c	write(31, '(a4,99f12.5)') 'lags',(lags(ispec),ispec=1,nspec)

C  Restore absent species in present phases
        do 81 lph=1,nph
         if (absent(lph)) go to 81
	 if (.not. allow(lph)) go to 81
         do 82 ispec=1,nspec
          if (f(lph,ispec) .eq. 0.) go to 82
	  if (spinod(ispec)) go to 82
	  if (resall .and. lags(ispec) .lt. lagsmax) go to 83
c	  if (resall) go to 83
	  if (lags(ispec) .gt. 0.) go to 82
83	  continue
          if (absents(ispec)) then
	   reslog = .true.
           absents(ispec) = .false.
           dn(ispec) = tradd
	   print '(1x,a17,i5,2x,a4,f12.5)', 'Restoring species',ispec,sname(ispec),lags(ispec)
	   write(31,'(1x,a17,i5,2x,a4,f12.5)') 'Restoring species',ispec,sname(ispec),lags(ispec)
          end if
82       continue
81      continue

        return
        end
