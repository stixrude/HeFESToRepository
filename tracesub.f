        subroutine tracesub(n,iphase,mphase,nspec,nph,absent,absents,add)

        include 'P1'

	integer nspec,nph,i,iph,itrmin,itrsmin
        double precision n(nspecp)
        logical absents(nspecp),absent(nphasep)
        integer iphase(nphasep),mphase(nphasep)
        logical add
        add = .false.

        call tform(n,iphase,mphase,itrmin,itrsmin,nspec,nph,absent,absents)

C  Remove least abundant trace species
	if (itrsmin .ne. 0) then
         absents(itrsmin) = .true. 
         add = .true.
         n(itrsmin) = 0.
         print*, 'Subtracting species in tracesub',itrsmin
         write(31,*) 'Subtracting species in tracesub',itrsmin
	end if

C  Remove phases containing only absent species
	do 11 iph=1,nph
	 absent(iph) = .true.
	 do 12 i=iphase(iph),iphase(iph)+mphase(iph)-1
	  absent(iph) = absent(iph) .and. absents(i)
12	 continue
11	continue

        return
        end
