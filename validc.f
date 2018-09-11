        logical function validc(vsum,ncp,dconstr,nconstr)

C  Determine validity of solution: amount of each component on each site in each phase is non-negative

        include 'P1'
        include 'chem.inc'
        include 'absent.inc'

	integer nconstr,i,ic,iph,ispec,j,jsp,kst,nsitecp
	double precision vsum,constr
        double precision nikp,ncp(nspecp)
	double precision dconstr(nspecp,nspecp),dconstrtmp(nspecp)
        double precision, parameter :: vsmall = 1.e-15

	do 6 i=1,nspecp
	 dconstrtmp(i) = 0.
	 do 6 j=1,nspecp
	  dconstr(i,j) = 0.
6	continue
	nconstr = 0
        call nform(ncp,n,n1,q2,nspec,nnull)
c        write(31,'(a9,71f9.5)') 'validc',(n(i),i=1,nspec)
        validc = .true.
        vsum = 0.
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
	    constr = -nikp
	    do 5 jsp=1,nspec
5	    dconstr(jsp,nconstr) = -dconstrtmp(jsp)
c	    write(31,'(a6,3i5,99f12.5)') 'valid',iph,kst,ic,nikp
c	    write(31,'(a6,i5,99f12.5)') 'constr',nconstr,constr,(dconstr(jsp,nconstr),jsp=1,nspec)
	   end if
           if (nikp .lt. -vsmall) then
            validc = .false.
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
            validc = .false.
	    vsum = vsum + 0.5*(abs(n(ispec)) - n(ispec))
           end if
	   nconstr = nconstr + 1
	   constr = -n(ispec)
	   dconstr(ispec,nconstr) = -1.0
c	   write(31,'(a6,i5,99f12.5)') 'valids',ispec,n(ispec)
c	   write(31,'(a6,i5,99f12.5)') 'constr',nconstr,constr,(dconstr(jsp,nconstr),jsp=1,nspec)
          end if
11       continue
10      continue

c       print '(a5,31e12.5,e12.5)', 'valid',(n(i),i=1,nspec),vsum

c        print*, 'valid',vsum
c        print '(71f12.4)', (ncp(i),i=1,nnull)

        return
        end
