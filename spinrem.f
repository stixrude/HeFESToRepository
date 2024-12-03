        subroutine spinrem(f,nspec,nph,absent,absents,allow,spinlog)

        include 'P1'
	include 'lag.inc'
	include 'numpar.inc'

	integer nspec,nph,ispec,iph
        double precision f(nphasep,nspecp),n(nspecp)
        logical absents(nspecp),absent(nphasep)
        logical spinod(nspecp),spinph(nphasep),allow(nphasep),spinlog
        common /spinc/ spinod,spinph

	spinlog = .false.

C  Remove present species that have experienced spinodal instability
        do 81 iph=1,nph
         if (absent(iph)) go to 81
         do 82 ispec=1,nspec
          if (f(iph,ispec) .eq. 0.) go to 82
          if (.not. absents(ispec)) then
           if (spinod(ispec)) then
c	    write(31,*) 'WARNING: Retaining spinodally unstable species',ispec,n(ispec)
c	    go to 82
	    write(31,*) 'Removing spinodal instability',ispec
            absents(ispec) = .true.
	    n(ispec) = 0.
	    allow(iph) = .false.
	    spinlog = .true.
	   end if
          end if
82       continue
81      continue

C  Check for phases with all spinodally unstable species
        do 74 iph=1,nph
         spinph(iph) = .true.
         do 75 ispec=iphase(iph),iphase(iph)+mphase(iph)-1
          spinph(iph) = spinph(iph) .and. spinod(ispec)
75       continue
         if (spinph(iph)) write(31,'(a34,1x,i5)') 'WARNING: Spinodally unstable phase',iph
74      continue

        return
        end
