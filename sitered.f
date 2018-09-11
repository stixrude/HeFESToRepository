        subroutine sitered(r,wreg,f,nsite,iphase,mophase,iastate,nph,nspec,nco)

        include 'P1'

	integer nph,nspec,nco,i,j,j2,k,kbase,kk,kr,lph,nk,ispec,jspec
        logical other
	logical iastate(nspecp,nspecp,nsitep),iastatet(nspecp,nspecp,nsitep)
        double precision r(ncompp,nspecp,nsitep)
        double precision f(nphasep,nspecp)
	integer nsite(nphasep)
        integer iphase(nphasep),mophase(nphasep)
        double precision rt(ncompp,nspecp,nsitep)
        double precision wreg(nphasep,nsitep,nspecp,nspecp)
        double precision wregt(nphasep,nsitep,nspecp,nspecp)
        integer nred(nphasep),ired(nphasep,nsitep)

        do 11 lph=1,nph
         do 11 i=1,nspec
          if (f(lph,i) .eq. 0.0) go to 11
          do 13 k=1,nsite(lph)
           do 12 j=1,nco
            do 12 j2=1,nco
             rt(j,i,k) = r(j,i,k)
             wregt(lph,k,j,j2) = wreg(lph,k,j,j2)
12         continue
	   do 14 jspec=1,nspec
	    iastatet(i,jspec,k) = iastate(i,jspec,k)
14	   continue
13	  continue
11      continue

        do 1 lph=1,nph
         nred(lph) = 0
         do 2 k=1,nsite(lph)
          do 3 i=1,nspec
           if (f(lph,i) .eq. 0) go to 3
           other = .false.
           do 4 j=1,nco
c            if (r(j,i,k) .ne. r(j,iphase(lph),k)) go to 2
c            if (r(j,i,k) .ne. 0.0) then
            if (r(j,i,k) .ne. r(j,iphase(lph),k) .and. r(j,i,k) .gt. 0.0) go to 2
            if (r(j,i,k) .gt. 0.0) then
             if (other) go to 2
             other = .true.
            end if
4          continue
3         continue
	  do 41 i=1,nspec
	   do 41 j=1,nspec
	    if (f(lph,i) .eq. 0) go to 41
	    if (f(lph,j) .eq. 0) go to 41 
	    if (iastate(i,j,k)) then
	     write(31,*) 'sitered Found site to avoid reduction',i,j,k
	     go to 2
	    end if
41	  continue
c	  if (mophase(lph) .eq. 0) then
           nred(lph) = nred(lph) + 1
           ired(lph,nred(lph)) = k
c	  end if
2        continue
1       continue

        do 21 lph=1,nph
	 print*, 'in sitered',lph,nred(lph)
         do 22 kr=1,nred(lph)
          k = ired(lph,kr)
          kbase = k - kr
          nk = 0
          do 222 kk=k+1,nsite(lph)
           nk = nk + 1
           do 23 i=1,nspec
            if (f(lph,i) .eq. 0) go to 23
            do 24 j=1,nco
             do 24 j2=1,nco
              r(j,i,kbase+nk) = rt(j,i,kk)
              wreg(lph,kbase+nk,j,j2) = wregt(lph,kk,j,j2)
24          continue
	    do 25 jspec=1,nspec
	     iastate(i,jspec,kbase+nk) = iastatet(i,jspec,kk)
25	    continue
23         continue
222       continue
22       continue
21      continue

        do 26 lph=1,nph
         nsite(lph) = nsite(lph) - nred(lph)
26      continue

        do 27 lph=1,nph
         do 27 ispec=1,nspec
          if (f(lph,ispec) .eq. 0.0) go to 27
           do 28 k=nsite(lph)+1,nsitep
            do 28 j=1,nco
             do 28 j2=1,nco
              r(j,ispec,k) = 0.0
              wreg(lph,k,j,j2) = 0.0
28         continue
27      continue

        do 31 lph=1,nph
         do 32 k=1,nsite(lph)
          do 33 i=1,nspec
           if (f(lph,i) .eq. 0) go to 33
           do 34 j=1,nco
            write(7,*) lph,k,i,j,r(j,i,k)
34         continue
33        continue
32       continue
31      continue

        return
        end
