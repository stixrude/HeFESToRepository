        subroutine gibmin(nnew,fret,itersum)

        include 'P1'
        include 'const.inc'
        include 'chem.inc'
        include 'absent.inc'
        include 'lag.inc'
        logical add,valid

	integer itersum,i,icase,ispec1,ispec2,iter,k,lph1,lph2,lphrep1,lphrep2,ndim
	integer nvet,nvep,ires
	double precision fret,apar,dhdpmol,dhdtmol,dvdpmol,dsdtmol,ehugo,feas,freths,fretnl,fretsav,ftol,phugo
	double precision Pi,qual,qualsav,vtarg,starg,Ti,tlast,vhugo,vsum,wmagg,Tfeas,Pfeas,func,yy,val
        double precision nnew(nspecp)
        double precision nnewsav(nspecp)
	double precision vec1(nspecp),vec2(nspecp),dnrm2,ddot,test
        double precision cpa(nspecp)
	double precision start,finish,ggibmin
        logical absentsav(nspecp),absentssav(nspecp),succes
	logical chcalc,adcalc,hucalc
        character*80 phname(nphasep),sname(nspecp)
        common /names/ phname,sname
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        common /chempot/ cpa
        double precision, parameter :: ftolcg=1.e-10, ftolqn=1.e-9, ntol=1.e-4, ftolhs=1.e-10, ftollm=1.e-8
        external func,dfunc,hesspot
	data ggibmin/0.0/
	call cpu_time(start)

	print*, 'in gibmin'

	fretsav = +1.e15
	qualsav = +1.e15
C  with large values of the ftol force all trace species to be removed.
c	ftol = 10.
	icase = 0

	write(31,*) 'enter gibmin',valid(vsum,nnew)
c        write(31,'(71f9.5)') (n(i),i=1,nspec)
	itersum = 0
	if (.not. valid(vsum,nnew)) call nlfeas(n,nnew,nspec,feas,icase,ione)
	if (icase .ne. 0) then
	 write(31,*) 'WARNING: gibmin starting with invalid solution'
c         write(31,'(71e12.5)') (n(i),i=1,nspec)
c         write(31,*) (nnew(i),i=1,nnull)
	end if
103     continue

	ndim = nnull
	if (adcalc) then
	 ndim = ndim + 1
	 nnew(ndim) = Ti
	end if
	if (chcalc) then
	 ndim = ndim + 1
	 nnew(ndim) = Pi
	end if
	if (adcalc) call Tlfeas(Tfeas,nnew)
	if (chcalc) call Plfeas(Pfeas,nnew)
	if (chcalc) write(31,*) 'Back from Plfeas = ',Pfeas
	if (adcalc) write(31,*) 'Back from Tlfeas = ',Tfeas
	call nlmin_L(nnew,ndim,fretnl,iter,ires)
	fret = fretnl
        itersum = itersum + iter
        call nform(nnew,n,n1,q2,nspec,nnull)
	if (adcalc) Ti = nnew(nnull+nvet)
	if (chcalc) Pi = nnew(nnull+nvep)
        write(31,'(a6,5i5,99f12.5)') 'gibmn',iter,nnull,nvep,nvet,ndim,Pi,Ti,(n(i),i=1,nspec)
c        write(31,'(a6,5i5,99f12.5)') 'vspec',iter,nnull,nvep,nvet,ndim,Pi,Ti,(vspeca(i),i=1,nspec)
	
	write(31,*) 'Calling ssave from gibmin 1'
	call ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)

        call tracesub(n,iphase,mphase,nspec,nph,absent,absents,add)
        if (add) then
         write(31,*) 'INFORMATION: Continuing in gibmin after trace removal',fret,qual,qtest
         call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
         call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
         go to 103
	end if

c        call hessmin(nnew,nnull,ftolhs,iter,freths)
c        fret = freths
c        itersum = itersum + iter

	ftol = 1.e-5

	write(31,*) 'Calling ssave from gibmin 2'
        call ssave(fret,qual,nnew,fretsav,qualsav,absentsav,absentssav,nnewsav,ftol,succes,ires)

	fret = fretsav
	qualsav = qual
        do 34 i=1,nspec
         absents(i) = absentssav(i)
         nnew(i) = nnewsav(i)
34      continue
        do 37 i=1,nph
37      absent(i) = absentsav(i)

        call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)

c	return

C  Check for artificial isochemical repetition of phases
        do 21 lph1=1,nph-1
         do 22 lph2=lph1+1,nph
          if (phname(lph1) .ne. phname(lph2)) go to 22
          if (absent(lph2)) go to 22
          if (absent(lph1)) then
	   write(31,*) 'Found second member of repeated phase artificially populated',lph1,lph2,phname(lph1),phname(lph2)
     &       ,absent(lph1),absent(lph2),(n(i),i=1,nspec)
	   go to 11
	  end if
	  write(31,*) 'Found repeated phase',lph1,lph2,phname(lph1),phname(lph2)
          k = 0
          do 23 ispec1=iphase(lph1),iphase(lph1)+mphase(lph1)-1
           ispec2 = iphase(lph2) + k
	   vec1(k+1) = n(ispec1)
	   vec2(k+1) = n(ispec2)
           k = k + 1
23        continue
	  test = ddot(mphase(lph1),vec1,ione,vec2,ione)/(dnrm2(mphase(lph1),vec1,ione)*dnrm2(mphase(lph1),vec2,ione)) - 1.
	  write(31,*) 'Phase compositions',(vec1(k),k=1,mphase(lph1)),(vec2(k),k=1,mphase(lph2)),test
	  if (abs(test) .gt. ntol) go to 22
C  Move phase 2 to phase 1.  set all phase 2 species to absent.
11	  write(31,*) 'Removing phase',lph2
	  if (absent(lph1)) write(31,*) 'Adding phase',lph1
          k = 0
          do 24 ispec1=iphase(lph1),iphase(lph1)+mphase(lph1)-1
           ispec2 = iphase(lph2) + k
	   absents(ispec1) = absents(ispec2)
	   absents(ispec2) = .true.
           n(ispec1) = n(ispec1) + n(ispec2)
           n(ispec2) = 0.
           k = k + 1
24        continue
	  absent(lph2) = .true.
	  absent(lph1) = .false.
c          call tracesub(n,iphase,mphase,nspec,nph,absent,absents,add)
          call sform(s,b,n1,q1,q2,nspec,nco,nc,ncs,nnull,nnulls,absents)
          call newfrm(q2,n,n1,nnew,nspec,nnull,absents)
	  if (adcalc) nnew(nnull+nvet) = Ti
	  if (chcalc) nnew(nnull+nvep) = Pi
c          go to 103
22	 continue
21	continue

	call cpu_time(finish)
	ggibmin = ggibmin + finish - start
	val = func(nnew)
c	write(31,*) 'time in gibmin',ggibmin,nnull,(nnew(i),i=1,nnull),fret,val,(n(i),i=1,nspec)
	write(31,*) 'time in gibmin',ggibmin,nnull

c	print*, 'end gibmin'

        return
        end
