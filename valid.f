        logical function valid(vsum,ncp)

C  Determine validity of solution: amount of each component on each site in each phase is non-negative

        include 'P1'
        include 'chem.inc'
        include 'absent.inc'

	integer ic,iph,ispec,jsp,kst,nconstr,nsitecp,i
	double precision vsum
        double precision nikp,ncp(nspecp)
	double precision constr(nspecp),dconstr(nspecp,nspecp),dconstrtmp(nspecp)
        double precision, parameter :: vsmall = 1.e-15

	nconstr = 0
	if (nnull .lt. 0) then
	 valid = .false.
	 return
	end if
        call nform(ncp,n,n1,q2,nspec,nnull)
        valid = .true.
        vsum = 0.

C---> Constrain the amount of each present species to be greater than zero
c        do 21 ispec=1,nspec
c         if (absents(ispec)) go to 21
c         if (n(ispec) .lt. -vsmall) then
c          valid = .false.
c          vsum = vsum + 0.5*(abs(n(ispec)) - n(ispec))
c         end if
c21      continue
c	return
C<---

	do 10 iph=1,nph
c	 if (absent(iph)) go to 10
         nsitecp = nsite(iph)
         do 1 kst=1,nsitecp
          do 3 ic=1,nco
           nikp = 0.0
           do 4 jsp=1,nspec
            if (absents(jsp)) go to 4
            nikp = nikp + f(iph,jsp)*r(ic,jsp,kst)*n(jsp)
	    dconstrtmp(jsp) = f(iph,jsp)*r(ic,jsp,kst)
4          continue
	   if (nikp .ne. 0.) then
	    nconstr = nconstr + 1
	    constr(nconstr) = nikp
	    do 5 jsp=1,nspec
5	    dconstr(nconstr,jsp) = dconstrtmp(jsp)
c	    write(31,'(a6,3i5,99f12.5)') 'valid',iph,kst,ic,nikp
c	    write(31,'(a6,i5,99f12.5)') 'constr',nconstr,constr(nconstr),(dconstr(nconstr,jsp),jsp=1,nspec)
	   end if
           if (nikp .lt. -vsmall) then
            valid = .false.
	    vsum = vsum + 0.5*(abs(nikp) - nikp)
c           print*, 'Negative nikp',iph,kst,ic
           end if
3         continue
1        continue
	 do 11 ispec=1,nspec
          if (absents(ispec)) go to 11
	  if (f(iph,ispec) .eq. 0) go to 11
          if (nsitecp .eq. 0 .or. nsitsp(ispec) .gt. 0) then
           if (n(ispec) .lt. -vsmall) then
            valid = .false.
	    vsum = vsum + 0.5*(abs(n(ispec)) - n(ispec))
           end if
          end if
11       continue
10      continue

c       print '(a5,31e12.5,e12.5)', 'valid',(n(i),i=1,nspec),vsum
c       write(31,'(a5,31e12.5,e12.5)') 'valid',(n(i),i=1,nspec),vsum

c        print*, 'valid',vsum
c        print '(71f12.4)', (ncp(i),i=1,nnull)

        return
        end
