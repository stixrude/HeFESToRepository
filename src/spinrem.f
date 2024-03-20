        subroutine spinrem(f,nspec,nph,absent,absents,allow,spinlog)

        include 'P1'
	include 'numpar.inc'

	integer nspec,nph,ispec,lph
        double precision f(nphasep,nspecp),n(nspecp)
        logical absents(nspecp),absent(nphasep)
        logical spinod(nspecp),spinph(nphasep),allow(nphasep),spinlog
        common /spinc/ spinod,spinph
	spinlog = .false.

C  Remove present species that have experienced spinodal instability
        do 81 lph=1,nph
         if (absent(lph)) go to 81
         do 82 ispec=1,nspec
          if (f(lph,ispec) .eq. 0.) go to 82
          if (.not. absents(ispec)) then
           if (spinod(ispec)) then
	    write(31,*) 'Removing spinodal instability',ispec
            absents(ispec) = .true.
	    n(ispec) = 0.
	    allow(lph) = .false.
	    spinlog = .true.
	   end if
          end if
82       continue
81      continue

        return
        end
