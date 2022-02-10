        subroutine petsub(nnew,itersum,iprint)

        include 'P1'
	include 'numpar.inc'
        include 'const.inc'
        include 'chem.inc'
        include 'absent.inc'
	include 'lag.inc'
	
	integer itersum,iprint,i,iadd,ic,icase,iph,ires,ispec,iter,loop,lphmin,lphtmp,nadd,ncall,nphprestmp,iextra
	double precision apar,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,feas,fret,fretmin,fretsav,ftol,phugo,Pi,prdiff
	double precision qual,qualmax,qualsav,vtarg,starg,Ti,tlast,val,vhugo,vsum,gspec,func,wmagg
	double precision start,finish,gphadd,gpetsub,pstart,pfinish
        logical add,reset,allow(nphasep),valid,resall,vtest
	logical chcalc,adcalc,hucalc,adcalcsav,chcalcsav
        double precision nnew(nspecp)
        double precision nnewsav(nspecp)
	double precision bold(ncompp)
        integer jphase(nspecp),nvet,nvep
        logical absentsav(nspecp),absentssav(nspecp)
        logical absenttmp(nspecp),absentstmp(nspecp)
	logical spinod(nspecp),spinph(nphasep),spinlog,reslog,succes
        integer, parameter :: iaddmax=10,loopmax=1
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /spinc/ spinod,spinph
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
	save ncall
        data ncall/0/
        data qualmax/-1.e15/
	data gphadd/0.0/,gpetsub/0.0/
	call cpu_time(pstart)

	adcalcsav = adcalc
	chcalcsav = chcalc
        fretsav= +1.e15
        qualsav = +1.e15
	ftol = 1.e-5
        ncall = ncall + 1
        itersum = 0
        iter = 0
        loop = 1
	if (ncall .eq. 1) then
	 do 291 ic=1,nco
291	 bold(ic) = 0.
	end if
	resall = .false.

        do 4 ispec=1,nspec
c4       gspeca(ispec) = gspec(ispec)
       gspeca(ispec) = gspec(ispec)
c	 write(31,*) ispec,gspeca(ispec)
4	continue
	write(31,*) 'spinodal instabilities:',Ti,(spinod(ispec),ispec=1,nspec)

cC  Check for phases with all spinodally unstable species
c        do 74 iph=1,nph
c         spinph(iph) = .true.
c         do 75 i=iphase(iph),iphase(iph)+mphase(iph)-1
c          spinph(iph) = spinph(iph) .and. spinod(i)
c75       continue
c         if (spinph(iph)) write(31,'(a34,1x,i5)') 'WARNING: Spinodally unstable phase',iph
c74      continue

C  Reset initial guess
	print*, 'Reset initial guess'

	reslog = .false.
        reset = .false.
	if (ncall .eq. 1) reset = .true.

10      continue
        do 14 i=1,nph
14      allow(i) = .true.
	do 292 ic=1,nco
	 if (b(ic) .ne. bold(ic)) reset = .true.
	 bold(ic) = b(ic)
292	continue
        if (reset) then
c         do 2 i=1,nspecp
         do 2 i=1,nnull
          nnew(i) = 0.
2	 continue
c	 if (chcalc) then
c	  nnew(nnull+nvep) = 10.0
c	  Pi = 10.0
c	 end if
        end if
	call nform(nnew,n,n1,q2,nspec,nnull)
        do 29 i=1,nspecp
	 dn(i) = 0.
29	continue

C  Find equilibrium assemblage

        print*, 'Initial guess',Pi,Ti,(absents(ispec),ispec=1,nspec)
        write(31,'(/,a,99f12.5)') 'Initial guess',Pi,Ti
	write(31,'(99f5.1)') (n(i),i=1,nspec)
	write(31,*) (absents(ispec),ispec=1,nspec)

	call spinrem(f,nspec,nph,absent,absents,allow,spinlog)
	if (reslog .or. spinlog) call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	if (.not. spinlog) go to 84
C  Check whether removal of spinodal instability generates infeasible solution and add phases as necessary
        if (.not. valid(vsum,nnew)) call nlfeas(n,nnew,nspec,feas,icase,ione)
	if (nnull .ne. nnulls .or. icase .ne. 0) then
         write(31,*) 'WARNING: sform failed after spinrem',nnull,nnulls,ires
         do 81 lph=1,nph
	  add = .false.
          if (.not. absent(lph)) go to 81
	  if (.not. allow(lph)) go to 81
	  if (spinph(lph)) go to 81
          write(31,*) 'Try adding phase',lph
	  add = .true.
	  absent(lph) = .false.
          do 82 ispec=1,nspec
           if (f(lph,ispec) .eq. 0.) go to 82
           if (.not. spinod(ispec)) absents(ispec) = .false.
82	  continue
	  call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	  if (nnull .ne. nnulls) then
	   write(31,*) 'Adding phase',lph,' failed in sform',nnull,nnulls
	   go to 81
	  end if
 	  if (.not. valid(vsum,nnew)) call nlfeas(n,nnew,nspec,feas,icase,ione)
          if (icase .ne. 0) then
           write(31,*) 'feasible failed for',lph,' skipping',icase
           go to 81
          end if
	  go to 83
81	 continue
	 if (.not. add) then
	  write(31,*) 'ERROR: Failed to find feasible solution'
	  stop
	 end if
	end if
83	continue
84	continue
        call gibmin(nnew,fret,iter)
        itersum = itersum + iter
	write(31,*) 'Calling ssave from petsub 1'
        call ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)
        call writeout(nnew,qual,itersum,-iprint)

        write(31,'(/,a,2f12.5)') 'Initial restore',Pi,Ti
	call lagscomp(s,n,lags,lagc,gspeca,lphase,nspec,nco,absent,absents)
        call restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)
	call spinrem(f,nspec,nph,absent,absents,allow,spinlog)
	if (spinlog .or. reslog) then
	 call newfrm3(s,dn,n,nspec,nc,ncs)
	 call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
         call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
         call gibmin(nnew,fret,iter)
         itersum = itersum + iter
	 write(31,*) 'Calling ssave from petsub 2'
         call ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)
         call writeout(nnew,qual,itersum,-iprint)
	end if

C  Skip phase addition or subtraction
c       go to 28
C

C  Check for phase addition or subtraction

        do 11 iadd=1,iaddmax
         if (adcalc) then
          adcalc = .false.
          val = func(nnew)
         end if
         if (chcalc) then
          chcalc = .false.
          val = func(nnew)
         end if
	 call cpu_time(start)
         call phaseadd(add,allow)
	 call cpu_time(finish)
	 gphadd = gphadd + finish - start
	 write(31,*) 'phaseadd time',gphadd
	 adcalc = adcalcsav
	 chcalc = chcalcsav
c	write(31,*) 'Back from phaseadd',Ti,nnull,nnew(nnull),nnew(nnull+nvet),nnewsav(nnull),nnewsav(nnull+nvet),add,Pi
         if (.not. add) go to 12
	 call newfrm3(s,dn,n,nspec,nc,ncs)
         call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
         call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
         if (nnulls .gt. nnull) print*, 'WARNING: Not Spanned',Pi,Ti
         nadd = nadd + 1
         write(31,'(/,a)') 'Phase Addition or Subtraction'
c         write(32,'(/,a)') 'Phase Addition or Subtraction'
106      nphpres = 0
         do 92 lph=1,nph
          absenttmp(lph) = absent(lph)
          if (absent(lph)) go to 92
          nphpres = nphpres + 1
          jphase(nphpres) = lph
92       continue
C   Extra degree of freedom
	iextra = 0
	if (adcalc) iextra = iextra + 1
	if (chcalc) iextra = iextra + 1
C  
         prdiff = nphpres - nco - iextra
	 if (prdiff .le. 0.) go to 101

C  Phase Rule Violated
C  Remove phases one by one and check for minimum energy solution
         print*, 'Phase Rule Violation',Pi,Ti,nphpres,nco
         write(31,'(/,a,2f12.3,2i5,1x,271l1,/)') 'Phase Rule Violation',Pi,Ti,nphpres,nco,(absents(ispec),ispec=1,nspec)
         if (prdiff .gt. 1) write(31,*) 'WARNING: Compound Violation',prdiff
         do 93 ispec=1,nspec
93       absentstmp(ispec) = absents(ispec)
         nphprestmp = nphpres
         fretmin = +1.e15
	 do 961 ispec=1,nspec
961	 dn(ispec) = 0.
	 print*, 'Test removing the following phases',(jphase(lph),lph=1,nphpres)
	 write(31,*) 'Test removing the following phases',(jphase(lph),lph=1,nphpres)
         do 94 lph=1,nphpres
          do 921 lphtmp=1,nph
921       absent(lphtmp) = absenttmp(lphtmp)
          do 931 ispec=1,nspec
931       absents(ispec) = absentstmp(ispec)
          absent(jphase(lph)) = .true.
          write (31,*) 'Test removing phase',jphase(lph)
          do 95 ispec=1,nspec
           if (f(jphase(lph),ispec) .eq. 0.) go to 95
           absents(ispec) = .true.
95        continue
          call restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)
	  call spinrem(f,nspec,nph,absent,absents,allow,spinlog)
          call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
	  if (nnull .ne. nnulls) then
	   write(31,*) 'WARNING: sform failed for',jphase(lph),' absent; skipping',nnull,nnulls
	   go to 102
	  end if
          do 952 ispec=1,nspecp
           nnew(ispec) = 0.
	   dn(ispec) = 0.
952       continue
	  icase = 0
          if (.not. valid(vsum,nnew)) call nlfeas(n,nnew,nspec,feas,icase,izero)
          if (icase .ne. 0) then
           write(31,*) 'WARNING: feasible failed for',jphase(lph),' absent; skipping',vsum
           go to 102
          end if
          call gibmin(nnew,fret,iter)
          itersum = itersum + iter
          if (.not. valid(vsum,nnew)) then
           write(31,*) 'WARNING: Solution for',jphase(lph),' absent is invalid',vsum
           print*, 'WARNING: Solution for',jphase(lph),' absent is invalid',vsum
	   if (vsum .gt. ssmall) then
	    print*, 'skipping'
	    write(31,*) 'skipping'
	    go to 102
	   end if
          end if
          call writeout(nnew,qual,iter,-iprint)
          if (fret .lt. fretmin) then
           fretmin = fret
           lphmin = lph
          end if
          print*, 'Energy when ',jphase(lph),' is absent = ',fret
          write(31,*) 'Energy when ',jphase(lph),' is absent = ',fret
102       absent(jphase(lph)) = .false.
          do 96 ispec=1,nspec
           if (f(jphase(lph),ispec) .eq. 0.) go to 96
           absents(ispec) = .false.
96        continue
94       continue
         do 941 lphtmp=1,nph
941      absent(lphtmp) = absenttmp(lphtmp)
         do 951 ispec=1,nspec
951      absents(ispec) = absentstmp(ispec)
         absent(jphase(lphmin)) = .true.
         do 105 ispec=1,nspec
          if (f(jphase(lphmin),ispec) .eq. 0.) go to 105
          absents(ispec) = .true.
105      continue
         print*, 'Lowest energy absent phase',jphase(lphmin)
         write(31,*) 'Lowest energy absent phase',jphase(lphmin)
         write(31,*) (absents(ispec),ispec=1,nspec)
         prdiff = prdiff - 1
         if (prdiff .gt. 0.) then
          print*, 'Need to remove another phase'
          write(31,*) 'Need to remove another phase'
          go to 106
         end if
         print*, 'Recovered from phase rule violation'
         write(31,*) 'Recovered from phase rule violation'
	 call spinrem(f,nspec,nph,absent,absents,allow,spinlog)
         call restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)
         call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
         do 953 ispec=1,nspecp
          nnew(ispec) = 0.
          dn(ispec) = 0.
953      continue
         if (.not. valid(vsum,nnew)) call nlfeas(n,nnew,nspec,feas,icase,ione)

C  Phase Rule Not Violated
101      continue
         call gibmin(nnew,fret,iter)
         itersum = itersum + iter
	 write(31,*) 'Calling ssave from petsub 3'
         call ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)
         call writeout(nnew,qual,itersum,-iprint)
         if (.not. valid(vsum,nnew)) then
          write(31,*) 'WARNING: Solution is invalid',vsum
          print*, 'WARNING: Solution is invalid',vsum
         end if
	 call lagscomp(s,n,lags,lagc,gspeca,lphase,nspec,nco,absent,absents)
         call restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)
	 call spinrem(f,nspec,nph,absent,absents,allow,spinlog)
	 if (spinlog .or. reslog) then
	  call newfrm3(s,dn,n,nspec,nc,ncs)
          call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
          call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
          call gibmin(nnew,fret,iter)
          itersum = itersum + iter
          write(31,*) 'Calling ssave from petsub 4'
          call ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)
          call writeout(nnew,qual,itersum,-iprint)
	 end if
         if (.not. valid(vsum,nnew)) then
          write(31,*) 'WARNING: Solution is invalid',vsum
          print*, 'WARNING: Solution is invalid',vsum
         end if
11      continue
12      continue
	
28      continue

C  Check the following two lines 21/1/15
c	write(31,*) 'After 12 continue',Ti,nnull,nnew(nnull),nnew(nnull+nvet),nnewsav(nnull),nnewsav(nnull+nvet),add,Pi
	fretsav = fret
	qualsav = qual
        do 34 i=1,nspec
         absents(i) = absentssav(i)
         nnew(i) = nnewsav(i)
34      continue
        do 37 i=1,nph
37      absent(i) = absentsav(i)
c	write(31,*) 'After 37 continue',Ti,nnull,nnew(nnull),nnew(nnull+nvet),nnewsav(nnull),nnewsav(nnull+nvet),add,Pi
        call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
c	write(31,*) 'After sform',Ti,nnull,nnew(nnull),nnew(nnull+nvet),nnewsav(nnull),nnewsav(nnull+nvet),add,Pi
	if (adcalc) Ti = nnewsav(nnull+nvet)
	if (chcalc) Pi = nnewsav(nnull+nvep)
c	write(31,*) 'Before writeout',Ti,nnull,nnew(nnull),nnew(nnull+nvet),nnewsav(nnull),nnewsav(nnull+nvet),add,Pi
        call writeout(nnew,qual,itersum,-iprint)
	if (nnull .gt. 0) write(31,*) 'After writeout',Ti,nnull,nnew(nnull),nnew(nnull+nvet),nnewsav(nnull),nnewsav(nnull+nvet),add,Pi

C  Check solution for consistency

        nphpres = 0
        do 91 iph=1,nph
         if (absent(iph)) go to 91
         nphpres = nphpres + 1
91      continue

        if (nphpres .gt. nco+iextra) then
         write(99,100) 'WARNING: Phase Rule Violation',Pi,Ti,nphpres,nco
         write(31,100) 'WARNING: Phase Rule Violation',Pi,Ti,nphpres,nco
        end if
        if (iadd .eq. iaddmax+1) then
         write(31,*) 'WARNING: Phase Addition Process May Not Be Converged'
         write(31,*) 'WARNING: Try Increasing iaddmax'
        end if
        if (.not. valid(vsum,nnew)) then
         write(31,*) 'WARNING: Solution is invalid',vsum
        end if

        if (loop .gt. loopmax .and. qual .gt. qtest) then
         write(31,*) 'WARNING: Loop inconsistency in petsub',fret,qual,qtest
        end if
c        if (loop .le. loopmax .and. qual .gt. qtest) then
        if (loop .le. loopmax) then
	 resall = .true.
         call restore(f,dn,lags,nspec,nph,absent,absents,allow,reslog,resall)
	 resall = .false.
	 if (reslog .or. qual .gt. qtest) then
	  reset = .true.
	  write(31,*) 'INFORMATION: Looping back in petsub with restore and reset'
          loop = loop + 1
	  go to 10
	 end if
	end if
	reset = .false.

        call writeout(nnew,qual,itersum,iprint)
        qualmax = max(qual,qualmax)
        write(97,*) Pi,Ti,itersum,nadd,qual,qualmax
	call cpu_time(pfinish)
	gpetsub = gpetsub + pfinish - pstart
	write(31,*) 'time in petsub',gpetsub

        return
100     format(a20,2f12.5,3i5)
        end
