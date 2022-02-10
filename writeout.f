        subroutine writeout(nnew,qual,iter,iprint)

        include 'P1'
        include 'chem.inc'
        include 'absent.inc'
        include 'lag.inc'
	include 'const.inc'

	integer iter,iprint,i,ispec,ispec1,ispec2,j,jspec,lph1,lph2,lphrep1,lphrep2,ncall,ndim,k,npres,nw,iextra,nvet,nvep
	integer ii,jj
	double precision qual,apar,dhdpmol,dhdtmol,diff11,diff12,dvdpmol,dsdtmol,ehugo,ftest,gibbs,gprint,phugo,Pi,pold,vtarg,starg
	double precision sumn,Ti,tlast,told,vhugo,func,depth,vsum,wmagg,nswap(nspecp),dnrm2
        logical newpt,valid
	logical chcalc,adcalc,hucalc
        character*80 phname(nphasep),sname(nspecp)
        double precision nnew(nspecp),nnewp(nspecp),hess(nspecp,nspecp),qh(nspecp,nspecp),np(nspecp),nnewpp(nspecp)
        double precision cpa(nspecp),xi(nspecp)
        double precision nprint(nspecp),nold1,nold2
	double precision fvec(nspecp),fjac(nspecp,nspecp),fvecp(nspecp),err(nspecp)
        common /names/ phname,sname
        common /state/ apar(nspecp,nparp),Ti,Pi
        common /chempot/ cpa
        common /tfindc/ tlast,vtarg,starg,phugo,vhugo,ehugo,dvdpmol,dsdtmol,dhdpmol,dhdtmol,wmagg,chcalc,adcalc,hucalc,nvet,nvep
        data ncall/0/,Pold/-1.e15/,Told/-1.e15/
        ncall = ncall + 1
        sname(nspecp) = ' '
        newpt = .false.
        if (Pold .ne. Pi .or. Told .ne. Ti) newpt = .true.
        Pold = Pi
        Told = Ti

C Analysis
        nw = max(nnull,nco)
        do 3 i=nnull+1,nspecp
3       xi(i) = 0.0
        call nform(nnew,n,n1,q2,nspec,nnull)
	ndim = nnull
	if (adcalc) then
	 ndim = nnull + 1
	 nnew(nnull+nvet) = Ti
	end if
        gibbs = func(nnew)
        call dfunc(nnew,xi)
c        write(31,*) 'in writeout',nnull,(nnew(i),i=1,nnull),gibbs,(n(i),i=1,nspec)
c        gibbs = sfunc(nnew)
c        call dsfunc(nnew,xi)
c	call physub(nnew,rho,wmagg,freeagg,2)
	go to 999
C--> Check projected derivatives
	fvec(1) = gibbs
	do 91 i=1,ndim
91	fjac(1,i) = xi(i)
	write(31,*) 'Check Derivatives',nnull,ndim,adcalc
	write(31,*) 'x',(nnew(i),i=1,ndim)
	write(31,'(a13,71f29.15)') 'n',(n(i),i=1,nspec)
	write(31,'(a13,71f29.15)') 'cp',(cpa(i),i=1,nspec)
	write(31,*) 'fjac',(fjac(1,i),i=1,ndim)
	write(31,*) 'f',fvec(1)
	call chkderlps(ione,ndim,nnew,fvec,fjac,nspecp,nnewp,fvecp,ione,err)
	write(31,*) 'xp',(nnewp(i),i=1,ndim)
        fvecp(1) = func(nnewp)
c        fvecp(1) = sfunc(nnewp)
	write(31,'(a13,71f19.15)') 'n',(n(i),i=1,nspec)
	write(31,'(a13,71f29.15)') 'cp',(cpa(i),i=1,nspec)
	write(31,*) 'fp',fvecp(1)
	call chkderlps(ione,ndim,nnew,fvec,fjac,nspecp,nnewp,fvecp,itwo,err)
	write(31,*) 'Error',Pi,Ti,err(1)
c  test
	do 112 j=1,nnull+1
	 do 111 i=1,nnull+1
111	 nnewpp(i) = nnew(i)
	 nnewpp(j) = nnewp(j)
         ftest = func(nnewpp)
c         ftest = sfunc(nnewpp)
	 write(31,*) 'ftest',ftest,nnewpp(j) - nnew(j),(ftest-fvec(1))/(nnewpp(j) - nnew(j))
112	continue
c  test
        gibbs = func(nnew)
        call dfunc(nnew,xi)
	write(31,*) 'Hessian'
	do 92 ispec=1,nspec
	 call hessian(ispec,n,hess)
92	continue
        call dgemm('Transpose Q2','Normal H',nnull,nspec,nspec,one,q2,nspecp,hess,nspecp,zero,qh,nspecp)
        call dgemm('Normal H','Normal Q2',nspec,nnull,nspec,one,qh,nspecp,q2,nspecp,zero,fjac,nspecp)
	do 95 i=1,nnull
	 fvec(i) = xi(i)
95	continue
	write(31,*) 'x',(nnew(i),i=1,nnull)
	write(31,'(a13,71f19.15)') 'n',(n(i),i=1,nspec)
	write(31,'(a13,71f29.15)') 'cp',(cpa(i),i=1,nspec)
	write(31,*) 'f',(fvec(i),i=1,nnull)
	write(31,*) 'fjac',((fjac(i,j),j=1,nnull),i=1,nnull)
	call chkderlps(nnull,nnull,nnew,fvec,fjac,nspecp,nnewp,fvecp,ione,err)
        gibbs = func(nnewp)
        call dfunc(nnewp,xi)
	do 94 ispec=1,nspec
	 fvecp(ispec) = xi(ispec)
94	continue
	write(31,*) 'xp',(nnewp(i),i=1,nnull)
	write(31,'(a13,71f19.15)') 'n',(n(i),i=1,nspec)
	write(31,'(a13,71f29.15)') 'cp',(cpa(i),i=1,nspec)
	write(31,*) 'fp',(fvecp(i),i=1,nnull)
	call chkderlps(nnull,nnull,nnew,fvec,fjac,nspecp,nnewp,fvecp,itwo,err)
	write(31,*) 'Hessian error',Pi,Ti,err(1)
C<--
        call nform(nnew,n,n1,q2,nspec,nnull)
        gibbs = func(nnew)
        call dfunc(nnew,xi)
999	continue
	go to 998
C--> Check straight derivatives
	fvec(1) = gibbs
	do 901 i=1,nspec
901	fjac(1,i) = cpa(i)
	write(31,*) 'Check Derivatives'
	write(31,*) 'x',(n(i),i=1,nspec)
	write(31,*) 'fjac',(fjac(1,i),i=1,nspec)
	write(31,*) 'f',fvec(1)
	call chkderlps(ione,nspec,n,fvec,fjac,nspecp,np,fvecp,ione,err)
	write(31,*) 'xp',(np(i),i=1,nspec)
	call newfrm(q2,np,n1,nnewp,nspec,nnull)
        fvecp(1) = func(nnewp)
	write(31,*) 'fp',fvecp(1)
	call chkderlps(ione,nspec,n,fvec,fjac,nspecp,np,fvecp,itwo,err)
	write(31,*) 'Error',Pi,Ti,err(1)
        gibbs = func(nnew)
	write(31,*) 'Hessian'
	do 902 ispec=1,nspec
	 npres = npres + 1
	 call hessian(ispec,n,hess)
	 do 902 jspec=1,nspec
	  fjac(ispec,jspec) = hess(ispec,jspec)
902	continue
	do 905 i=1,nspec
	 fvec(i) = cpa(i)
905	continue
	write(31,*) 'x',(n(i),i=1,nspec)
	write(31,*) 'f',(fvec(i),i=1,nspec)
	write(31,*) 'fjac',((fjac(i,j),j=1,nspec),i=1,nspec)
	call chkderlps(nspec,nspec,n,fvec,fjac,nspecp,np,fvecp,ione,err)
        call newfrm(q2,np,n1,nnewp,nspec,nnull)
        gibbs = func(nnewp)
	do 904 ispec=1,nspec
	 fvecp(ispec) = cpa(ispec)
904	continue
	write(31,*) 'xp',(np(i),i=1,nspec)
	write(31,*) 'fp',(fvecp(i),i=1,nspec)
	call chkderlps(nspec,nspec,n,fvec,fjac,nspecp,np,fvecp,itwo,err)
	write(31,*) 'Error',Pi,Ti,err(1)
C<--
        call nform(nnew,n,n1,q2,nspec,nnull)
        gibbs = func(nnew)
        call dfunc(nnew,xi)
998	continue

        sumn = 0.0
        do 4 i=1,nspec
4       sumn = sumn + n(i)
	
        call qcalc(nnew,qual)
        do 1 i=1,nph
         nphase(i) = 0.0
         do 1 j=1,nspec
          nphase(i) = nphase(i) + f(i,j)*n(j)
1       continue
        do 22 ispec=1,nspec
22      nprint(ispec) = n(ispec)

	go to 121
        lphrep1 = 0
        lphrep2 = 0
        do 21 lph1=1,nph-1
         do 21 lph2=lph1+1,nph
          if (phname(lph1) .eq. phname(lph2)) then
           lphrep1 = lph1
           lphrep2 = lph2
          end if
21      continue
        if (lphrep1 .ne. 0 .and. lphrep2 .ne. 0) then
         diff11 = abs(n(iphase(lphrep1))-nold1)
         diff12 = abs(n(iphase(lphrep1))-nold2)
         if (abs(n(iphase(lphrep1))-nold1) .gt. abs(n(iphase(lphrep1))-nold2)) then
          k = 0
          do 23 ispec1=iphase(lphrep1),iphase(lphrep1)+mphase(lphrep1)-1
           ispec2 = iphase(lphrep2) + k
           nprint(ispec1) = n(ispec2)
           nprint(ispec2) = n(ispec1)
           k = k + 1
23        continue
         end if
        end if
121	continue

        do 61 lph1=1,nph-1
         if (phname(lph1) .ne. phname(lph1+1)) go to 61
         if (n(iphase(lph1)) .ge. n(iphase(lph1)+mphase(lph1))) go to 61
	 if (absents(iphase(lph1)+mphase(lph1))) go to 61
         print*, 'Swapping phase',lph1,Pi,Ti,n(iphase(lph1)),n(iphase(lph1)+mphase(lph1))
         do 62 jj=1,mphase(lph1)
          nprint(iphase(lph1)+jj-1) = n(iphase(lph1)+mphase(lph1)+jj-1)
          nprint(iphase(lph1)+mphase(lph1)+jj-1) = n(iphase(lph1)+jj-1)
62       continue
61      continue
          
C Output
        if (ncall .eq. 1) then
         write(99,'(a7,2a9,105a12)') 'Pi','depth','Ti',(sname(ispec)(1:9),ispec=1,nspec)
        end if
c        write(7,*) 'Pressure, Temperature',Pi,Ti
c        write(7,*) 'Solution = ',(nnew(i),i=1,nnull)
c        write(7,*) 'Composition = ',(n(i),i=1,nspec)
c        write(7,*) 'Number of Iterations = ',iter
c        write(7,*) 'Gibbs Free Energy = ',gibbs/1000.
c       if (iprint .lt. 1) then
        if (iprint .eq. -1) then
c         if (newpt) then
c          write(31,'(//,a)') '------------------------- Pressure (GPa), Depth (km), Temperature (K) -------------------------'
c          depthp = depth(Pi)
cc          write(31,700)  Pi,depth(Pi),Ti
c          write(31,700)  Pi,depthp,Ti
c         end if
	 iextra = 0
	 if (adcalc) iextra = iextra + 1
	 if (chcalc) iextra = iextra + 1
         write(31,'(/,a23,i20)') 'Number of Iterations = ',iter
         write(31,*) 'Gibbs Free Energy = ',gibbs/1000.
         write(31,*) 'Quality = ',qual 
         write(31,'(a11,35f8.4)') 'Solution   ',(nnew(i),i=1,nnull+iextra)
         write(31,'(a11,35f8.4)') 'Derivatives',(xi(i),i=1,nnull+iextra)
         write(31,*) (absents(i),i=1,nspec)
         write(31,*) 'nc,nco,ncs,nspec =',nc,nco,ncs,nspec
c         write(7,200) Pi,Ti,(nprint(i),i=1,nspec),gibbs/1000.,qual
c         write(31,'(a11,105f8.4)') 'Composition',(nprint(i),i=1,nspec),
c     &    gibbs/1000.,qual
c         write(31,'(a11,105f8.4)') 'Chem. Pots.',(cpa(i)/1000.,i=1,nspec)
        end if

c        if (iprint .ge. 1) then
        if (iprint .eq. 1) then
         write(31,'(/,a)') 'Phase equilibria'
         write(31,'(a7,a9,a9,105a12)') 'Pi','depth','Ti',(sname(ispec)(1:8),ispec=1,nspec)
         write(31,700) Pi,depth(Pi),Ti,(nprint(i),i=1,nspec),
     &    gibbs/1000.,qual
	 do 31 lph1=1,nph-1
	  do 31 lph2=lph1+1,nph
	   if (phname(lph1) .ne. phname(lph2)) go to 31
	   if (absent(lph1)) go to 31
	   if (absent(lph2)) go to 31
	   write(31,*) 'Exsolved phase ', Pi,Ti,phname(lph1)
31	 continue
         gprint = gibbs/1000.
         if (abs(gprint) .gt. 10.) gprint = gibbs/10000.
         write(99,700) Pi,depth(Pi),Ti,(nprint(i),i=1,nspec),
     &    gprint,qual
         if (lphrep1 .ne. 0 .and. lphrep2 .ne. 0) then
          nold1 = nprint(iphase(lphrep1))
          nold2 = nprint(iphase(lphrep2))
         end if
         if (qual .gt. qtest) then
c          write (99,300) 'Poor Solution',Pi,Ti,qual
          write (31,300) 'Poor Solution',Pi,Ti,qual
          print 300, 'Poor Solution',Pi,Ti,qual
         end if
         if (.not. valid(vsum,nnew)) then
          write (99,300) 'Invalid Solution',Pi,Ti,vsum
          write (31,300) 'Invalid Solution',Pi,Ti,vsum
          print 300, 'Invalid Solution',Pi,Ti,qual
         end if
        end if

        return
200     format(f8.2,f8.1,35f8.4)
300     format(a13,3f12.5)
c700     format(f7.2,f9.2,f9.2,105f9.5)
700     format(f7.2,f9.2,f9.2,105f12.8)
        end
